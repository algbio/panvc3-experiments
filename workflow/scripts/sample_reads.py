# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

from Bio import SeqIO
import argparse
import numpy.random as npr
import sys

if __name__ == "__main__":
	parser = argparse.ArgumentParser("Sample reads from a FASTQ file.")
	parser.add_argument('sampling-probability', type = float, help = "Sampling probability (using a uniform distribution).")
	parser.add_argument('seed', type = int, help = "Random seed.")
	args = parser.parse_args()

	npr.seed(args.seed)
	for rec in SeqIO.parse(sys.stdin, "fastq"):
		rr = npr.uniform()
		if rr < args.sampling_probability:
			SeqIO.write(rec, sys.stdout, "fastq")
