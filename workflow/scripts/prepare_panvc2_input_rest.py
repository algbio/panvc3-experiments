# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import gzip
import sys


def parse_input(fp, begin_contig, end_contig):
	is_first = True
	for line in fp:
		if line.startswith(">"):
			line = line.rstrip("\n")
			line = line[1:]
			fields = line.split("\t")
			contig_name = fields[0]

			if is_first:
				is_first = False
			else:
				end_contig()

			begin_contig(contig_name)
		else:
			yield line
	if not is_first:
		end_contig()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Convert vcf2multialign's non-founder seqeunce output for PanVC 2.")
	parser.add_argument('founder_count', type = int, help = "Number of founder sequences, effectively the number of copies of each sequence to be made.")
	parser.add_argument('output_prefix', type = str, help = "Output prefix")
	args = parser.parse_args()

	fp = None
	def do_open(name):
		fp = open(f"{args.output_prefix}{x}", "wx")

	for line in parse_input(sys.stdin, do_open, lambda: fp.close()):
		fp.write(line)
