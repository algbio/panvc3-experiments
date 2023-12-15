# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import itertools
import sys


def contig_order(fasta_index_path):
	def helper():
		with open(fasta_index_path, "r") as fp:
			for line in fp:
				line = line.rstrip("\n")
				fields = line.split("\t")
				yield fields[0]
	return {name: idx for idx, name in enumerate(helper())}


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Reorder the @SQ entries in the given SAM header as listed in the given FASTA index..")
	parser.add_argument('fasta_index', type = str, help = "The FASTA index path")
	args = parser.parse_args()

	order = contig_order(args.fasta_index)

	headers_1 = []
	headers_2 = []
	sq_headers = []

	seen_sq = False
	for line in sys.stdin:
		line = line.rstrip("\n")
		if line.startswith("@SQ\t"):
			seen_sq = True
			fields = line.split("\t")
			sn = None
			for ff in fields:
				if ff.startswith("SN:"):
					sn = ff[3:]
					break
			sq_headers.append((sn, line))
		elif seen_sq:
			headers_2.append(line)
		else:
			headers_1.append(line)
	
	sorted_sq_headers = [None] * len(sq_headers)
	for key, line in sq_headers:
		sorted_sq_headers[order[key]] = line

	for line in itertools.chain(headers_1, sorted_sq_headers, headers_2):
		sys.stdout.write(line)
		sys.stdout.write("\n")
