"""
a3m and a2m input for BioPython

This registers new iterators that can be used for reading sequences and MSAs with Bio.SeqIO and Bio.AlignIO.

To use, simply import this module in your project and then use e.g. Bio.AlignIO.parse using any of these new formats:

    a3m         Regular A3M format
    a2m         Regular A2M format

    a3m-nogaps  A3M format, ignoring all insert states
    a2m-nogaps  A2M format, ignoring all insert states

See the demo at the bottom of the file for details
"""

import Bio.SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaIterator
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import collections
import re

def pad(x, l, pad_char="-"):
    """Pad x with pad_char characters until it has length l"""
    return x + pad_char * (l - len(x))


def insert_gaps(sequences):
    """Return a iterator over sequences that corresponds to the input sequences with gaps added"""

    re_match = re.compile(r'([A-Z0-9~-])')

    # split input sequences into columns of single-character match (upper case and -) states and
    # variable-length inserts (lower case)
    inserts = [
        re_match.split(seq)
        for seq in sequences
    ]

    # for each column, determine maximum length
    insert_max_lengths = [
        max(
            len(inserts[i][j])
            for i in range(len(inserts))
        )
        for j in range(len(inserts[0]))
    ]

    return (
        "".join(
            pad(insert, insert_len)
            for insert, insert_len in zip(seq, insert_max_lengths)
        )
        for seq in inserts
    )


def A2MA3MIterator(format='a3m', remove_inserts=False):
    """Create a SeqIO-style iterator from parameters

    Arguments:
        format          -- What input format to parse ('a3m' or 'a2m') (default: 'a3m')
        remove_inserts  -- Whether inserts with respect to the query sequence should be removed (default: False)

    """
    def inner_iterator(handle, alphabet=single_letter_alphabet, title2ids=None):
        """Generator function to iterate over a3m records (as SeqRecord objects).

        Arguments:
            handle      -- input file
            alphabet    -- optional alphabet
            title2ids   -- A function that, when given the title of the FASTA file (without the beginning >),
                           will return the id, name and description (in that order) for the record as a tuple
                           of strings. If this is not given, then the entire title line will be used as the
                           description, and the first word as the id and name.
        """

        re_insert = re.compile(r'[a-z.]')

        if title2ids is None:
            title2ids = lambda t: (t.split(None, 1)[0], t.split(None, 1)[0], t)


        # piggyback on the fasta parser for splitting file into title and sequence
        parsed = collections.OrderedDict(SimpleFastaParser(handle))

        titles = parsed.keys()
        sequences = parsed.values()

        if format == 'a3m' and not remove_inserts:
            sequences = insert_gaps(sequences)

        elif format == 'a2m' and not remove_inserts:
            sequences = (
                seq.replace(".", "-")
                for seq in sequences
            )

        elif remove_inserts:
            sequences = (
                re_insert.sub('', seq)
                for seq in sequences
            )

        else:
            raise ValueError("Unknown format: {0}".format(format))

        for title, seq in zip(titles, sequences):

            id, name, description = title2ids(title)
            yield SeqRecord(Seq(seq, alphabet), id=id, name=name, description=description)

    return inner_iterator


def monkey_patch(name, iterator):
    """Monkey patch the Bio.SeqIO iterator registry if no such iterator exists yet"""

    if name not in Bio.SeqIO._FormatToIterator:
        Bio.SeqIO._FormatToIterator[name] = iterator

monkey_patch('a3m', A2MA3MIterator('a3m', False))
monkey_patch('a3m-nogaps', A2MA3MIterator('a3m', True))
monkey_patch('a2m', A2MA3MIterator('a2m', False))
monkey_patch('a2m-nogaps', A2MA3MIterator('a2m', True))


if __name__ == '__main__':

    def assert_equal(iter_a, iter_b, msg):
        for rec_a, rec_b in zip(iter_a, iter_b):
            # assert that string representations of SeqRecords match - this is
            # quite hacky and should not be done :)
            assert repr(rec_a) == repr(rec_b), msg + ": " + repr((rec_a, rec_b))


    # Compare our parser against output of reformat.pl from the HH-suite

    a2m = list(Bio.SeqIO.parse("data/test.a2m", "a2m"))
    a3m = list(Bio.SeqIO.parse("data/test.a3m", "a3m"))
    fas = list(Bio.SeqIO.parse("data/test.fasta", "fasta"))

    assert_equal(a2m, a3m, "a2m vs a3m")
    assert_equal(a3m, fas, "a3m vs fas")
    assert_equal(a2m, fas, "a2m vs fas")

    a2m_ng = list(Bio.SeqIO.parse("data/test.a2m", "a3m-nogaps"))
    a3m_ng = list(Bio.SeqIO.parse("data/test.a3m", "a3m-nogaps"))
    fas_ng = list(Bio.SeqIO.parse("data/test.fasta_nogaps", "fasta"))

    assert_equal(a2m_ng, a3m_ng, "a2m_ng vs a3m_ng")
    assert_equal(a3m_ng, fas_ng, "a3m_ng vs fas_ng")
    assert_equal(a2m_ng, fas_ng, "a2m_ng vs fas_ng")


