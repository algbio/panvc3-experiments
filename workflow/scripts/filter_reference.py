# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

from Bio import SeqIO, SeqRecord
import argparse
import sys


def process(contig_ids, contig_id_output_fh):

	def helper():
		for rec in SeqIO.parse(sys.stdin, "fasta"):
			if contig_id_output_fh:
				contig_id_output_fh.write(rec.id)
				contig_id_output_fh.write('\n')
			
			if rec.id not in contig_ids:
				yield SeqRecord.SeqRecord(rec.seq, id = rec.id)
	
	SeqIO.write(helper(), sys.stdout, "fasta-2line")


if __name__ == "__main__":
	parser = argparse.ArgumentParser("Output reference sequences as read from stdin with the given entries removed.")
	parser.add_argument('-c', '--contig', type = str, nargs = '+', help = "Contig identifier to be removed")
	parser.add_argument('-o', '--output-contig-ids', type = argparse.FileType('w'), help = "Output contig identifiers from the reference")
	args = parser.parse_args()
	
	process(frozenset(args.contig), args.output_contig_ids)
