# A3M / A2M input for BioPython

This module enables BioPython to parse the Soeding lab's custom HH-suite alignment formats A3M and A2M. The goal of the module is to replicate the results of the reformat.pl included in the HH-suite as faithfully as possible

## Usage

Simply `import A3MIO` in addition to your BioPython modules. This creates support for the following new formats in both `SeqIO` and `AlignIO` packages:

 * `a3m` - A3Ms with insertions
 * `a3m-nogaps` A3Ms that ignore insertions w.r.t. query sequence
 * `a2m` - A2Ms with insertions
 * `a2m-nogaps` A2Ms that ignore insertions w.r.t. query sequence


Example:

	import Bio.SeqIO
	import A3MIO

	for record in Bio.SeqIO.parse("example.a3m", "a3m"):
		print(record)

## Status
Not yet tested on larger sequence databases for accuracy or performance. Be sure to evaluate before you use it :)

## License
GNU GPL v3



