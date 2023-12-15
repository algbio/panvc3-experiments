# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import gzip
import sys


def parse_input(fp):
	contig_name = None
	sequence_lines = None
	is_first = True
	for line in fp:
		if line.startswith(">"):
			if is_first:
				is_first = False
			else:
				yield contig_name, sequence_lines

			line = line.rstrip("\n")
			line = line[1:]
			fields = line.split("\t")
			contig_name = fields[0]
			sequence_lines = []
		else:
			sequence_lines.append(line)
	if not is_first:
		yield contig_name, sequence_lines


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Convert vcf2multialign's non-founder seqeunce output for PanVC 2.")
	parser.add_argument('founder_count', type = int, help = "Number of founder sequences, effectively the number of copies of each sequence to be made.")
	parser.add_argument('output_prefix', type = str, help = "Output prefix")
	args = parser.parse_args()

	for contig_name, sequence_lines in parse_input(sys.stdin):
		with open(f"{args.output_prefix}{contig_name}", "x") as fp:
			fp.write(">REF\n")
			for line in sequence_lines:
				fp.write(line)
			for ii in range(args.founder_count):
				fp.write(f">{1 + ii}\n")
				for line in sequence_lines:
					fp.write(line)
