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

def pad(x, l):
    """Pad x with '-' characters until it has length l"""
    return x + "-" * (l - len(x))


def A2MA3MIterator(format='a3m', remove_inserts=False):
    def inner_iterator(handle, alphabet=single_letter_alphabet, title2ids=None):
        """Generator function to iterate over a3m records (as SeqRecord objects).

        handle - input file
        alphabet - optional alphabet
        title2ids - A function that, when given the title of the FASTA file (without the beginning >),
                    will return the id, name and description (in that order) for the record as a tuple
                    of strings. If this is not given, then the entire title line will be used as the
                    description, and the first word as the id and name.
        """

        re_match = re.compile(r'([A-Z0-9~-])')
        re_insert = re.compile(r'[a-z.]')

        if title2ids is None:
            title2ids = lambda t: (t.split(None, 1)[0], t.split(None, 1)[0], t)

        if format == 'a3m' and not remove_inserts:

            # piggyback on the fasta parser for splitting file into title and sequence
            parsed = collections.OrderedDict(SimpleFastaParser(handle))

            # split input sequences into columns of single-character match (upper case and -) states and
            # variable-length inserts (lower case)
            inserts = [
                re_match.split(seq)
                for seq in parsed.values()
            ]

            # for each column, determine maximum length
            insert_max_lengths = [
                max(
                    len(inserts[i][j])
                    for i in range(len(inserts))
                )
                for j in range(len(inserts[0]))
            ]

            # pad so that every element of a column has the same length
            for title, seq in zip(parsed.keys(), inserts):

                seq_padded = "".join(
                    pad(insert, insert_len)
                    for insert, insert_len in zip(seq, insert_max_lengths)
                )

                id, name, description = title2ids(title)
                yield SeqRecord(Seq(seq_padded, alphabet), id=id, name=name, description=description)

        elif format == 'a2m' and not remove_inserts:

            for title, seq in SimpleFastaParser(handle):
                seq = seq.replace(".", "-")

                id, name, description = title2ids(title)
                yield SeqRecord(Seq(seq, alphabet), id=id, name=name, description=description)

        else:

            for title, seq in SimpleFastaParser(handle):
                seq_noinserts = re_insert.sub('', seq)

                id, name, description = title2ids(title)
                yield SeqRecord(Seq(seq_noinserts, alphabet), id=id, name=name, description=description)

    return inner_iterator

# Monkey patch the Bio.SeqIO iterator registry
Bio.SeqIO._FormatToIterator['a3m'] = A2MA3MIterator('a3m', False)
Bio.SeqIO._FormatToIterator['a3m-nogaps'] = A2MA3MIterator('a3m', True)
Bio.SeqIO._FormatToIterator['a2m'] = A2MA3MIterator('a2m', False)
Bio.SeqIO._FormatToIterator['a2m-nogaps'] = A2MA3MIterator('a2m', True)

if __name__ == '__main__':

    def assert_equal(iter_a, iter_b, msg):
        for rec_a, rec_b in zip(iter_a, iter_b):
            assert repr(rec_a) == repr(rec_b), msg + ": " + repr((rec_a, rec_b))


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


