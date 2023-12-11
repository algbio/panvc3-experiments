# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import itertools
import math
import pysam
import sys


"Custom exception for stopping merge."
class StopMerge(Exception):
	pass


"Take the next value from the given iterator or exhaust the other iterator and pass the values to missing_fn."
def take_next(it, other_it, missing_fn):
	try:
		return next(it)
	except StopIteration:
		for val in other_it:
			missing_fn(val)
		raise StopMerge


"Merge the values from the two given iterables in such a way that non-matching values are passed to the given callbacks."
def merge(lhs_it, rhs_it, do_cmp, lhs_missing, rhs_missing):
	lhs_it = iter(lhs_it)
	rhs_it = iter(rhs_it)
	lhs = None
	rhs = None

	try:
		lhs = take_next(lhs_it, rhs_it, lhs_missing)
		rhs = take_next(rhs_it, itertools.chain([lhs], lhs_it), rhs_missing)

		while True:
			res = do_cmp(lhs, rhs)
			if 0 == res:
				yield (lhs, rhs)
				lhs = take_next(lhs_it, rhs_it, lhs_missing)
				rhs = take_next(rhs_it, itertools.chain([lhs], lhs_it), rhs_missing)
			elif res < 0:
				rhs_missing(lhs)
				lhs = take_next(lhs_it, itertools.chain([rhs], rhs_it), lhs_missing)
			else:
				lhs_missing(rhs)
				rhs = take_next(rhs_it, itertools.chain([lhs], lhs_it), rhs_missing)
	except StopMerge:
		return


"Determine the opening mode for Pysam from the file name extension."
def open_mode(path):
	return "rb" if path.lower().endswith(".bam") else "r"


"Check that the given alignment file is sorted by QNAME."
def check_sort_order(alignment_file, path):
	header = alignment_file.header.to_dict()
	hd = header.get('HD')
	if hd is None:
		print(f"ERROR: no @HD header in {path}.", file = sys.stderr)
		sys.exit(1)
	so = hd.get('SO')
	if so is None:
		print(f"ERROR: no sort order specified in {path}.", file = sys.stderr)
		sys.exit(1)
	if "queryname" != so:
		print(f"ERROR: {path} is sorted by {so} instead of QNAME.", file = sys.stderr)
		sys.exit(1)


"cmp() replacement."
def cmp_str(lhs, rhs):
	if lhs < rhs:
		return -1
	elif rhs < lhs:
		return 1
	return 0


"Partition the alignments by the segment index"
def partition_by_aligned_segment(alns):
	first = []
	second = []
	other = []
	for aln in alns:
		if aln.is_read1:
			first.append(aln)
		elif aln.is_read2:
			second.append(aln)
		else:
			other.append(aln)
	return first, second, other


"Check that the count of the true alignments is sensible."
def check_true_alns(alns, qname, seg_idx):
	size = len(alns)
	if 0 == size:
		print(f"WARNING: No alignments in truth for read {qname} segment {seg_idx}.", file = sys.stderr)
		return False
	if 1 < size:
		print(f"WARNING: Multiple alignments in truth for read {qname} segment {seg_idx}.", file = sys.stderr)
		return False
	return True


"Check that the count of the tested alignments is sensible."
def check_tested_alns(alns, qname, seg_idx, should_report_missing):
	if 0 == len(alns):
		if should_report_missing:
			print(f"WARNING: No alignments in tested for read {qname} segment {seg_idx}.", file = sys.stderr)
		return False
	return True


def sanity_check_truth(qname, true_1, true_2, true_other):
	if not check_true_alns(true_1, qname, 1):
		return False
	if not check_true_alns(true_2, qname, 2):
		return False
	if 0 != len(true_other):
		print(f"WARNING: Additional aligned segments in truth for read {qname}.", file = sys.stderr)
		return False
	return True


def sanity_check_tested(qname, tested_1, tested_2, tested_other, should_report_missing):
	check = lambda tested, seg_idx: check_tested_alns(tested, qname, seg_idx, should_report_missing)
	if not check(tested_1, 1):
		return False
	if not check(tested_2, 2):
		return False
	if 0 != len(tested_other):
		print(f"WARNING: Additional aligned segments in test for read {qname}.", file = sys.stderr)
		# Continue.
	return True


"Aligned position of the given segment"
def alignment_position(aln):
	retval = aln.reference_start

	# Check for hard clipping since it affects the position.
	cigar = aln.cigartuples
	first_item = cigar[0]
	if pysam.CHARD_CLIP.value == first_item[0]:
		retval -= first_item[1]

	return retval


def missing_from_truth(aln_group):
	print(f"WARNING: Read {aln_group[0].query_name} is missing from the truth.", file = sys.stderr)


def fp_div(numerator, denominator):
	if 0 == denominator:
		return math.nan
	return numerator / denominator


def ref_name_mappings(path):
	if path is None:
		return {}
	
	def helper():
		with open(path, "r") as fp:
			for line in fp:
				if line.startswith("#"):
					continue
				line = line.rstrip("\n")
				fields = line.split("\t")
				yield fields[0], fields[1]
	return {k: v for k, v in helper()}


"""
Calculate the precision and the recall of the given alignments.

We mostly follow the method in Robert Lindner, Caroline C. Friedel. A Comprehensive Evaluation of Alignment Algorithms in the Context of RNA-Seq.
"""
def calculate_precision_and_recall(
	truth_path,
	tested_path,
	distance_threshold,
	truth_rname_mapping_path,
	should_report_missing_in_tested
):
	truth_rname_mappings = ref_name_mappings(truth_rname_mapping_path)

	total_reads = 0
	total_alns = 0
	true_positives = 0
	false_positives = 0
	false_negatives = 0

	def handle_missing_from_tested(truth_eqc):
		qname = truth_eqc[0].query_name
		true_1, true_2, true_other = partition_by_aligned_segment(truth_eqc)
		# Sanity checks.
		if not sanity_check_truth(qname, true_1, true_2, true_other):
			return

		# There should be one alignment in both true_1 and true_2.
		nonlocal total_reads
		nonlocal false_negatives
		total_reads += 2
		false_negatives += 2

	with pysam.AlignmentFile(truth_path, open_mode(truth_path)) as truth, pysam.AlignmentFile(tested_path, open_mode(tested_path)) as tested:
		check_sort_order(truth, truth_path)
		check_sort_order(tested, tested_path)
		
		# Map the reference names to indices.
		tested_ref_ids = {rname: idx for idx, rname in enumerate(tested.header.references)}

		# Make equivalence classes by QNAME lazily.
		do_map = lambda aln_file: map(lambda x: list(x[1]), itertools.groupby(iter(aln_file), lambda x: x.query_name))
		truth_alns	= do_map(truth)
		tested_alns	= do_map(tested)

		# Classify.
		seen_non_matching_reference_names = set()
		merged = merge(
			truth_alns,
			tested_alns,
			lambda lhs, rhs: cmp_str(lhs[0].query_name, rhs[0].query_name),
			missing_from_truth,
			handle_missing_from_tested
		)
		for truth_eqc, tested_eqc in merged:
			qname = truth_eqc[0].query_name
			true_1, true_2, true_other = partition_by_aligned_segment(truth_eqc)
			tested_1, tested_2, tested_other = partition_by_aligned_segment(tested_eqc)

			# Sanity checks.
			if not sanity_check_truth(qname, true_1, true_2, true_other):
				continue
			if not sanity_check_tested(qname, tested_1, tested_2, tested_other, should_report_missing_in_tested):
				continue

			for seg_idx, true, tested in [(1, true_1[0], tested_1), (2, true_2[0], tested_2)]:
				if not true.is_mapped:
					print(f"WARNING: Truth unmapped for read {qname} segment {seg_idx}.", file = sys.stderr)
					continue

				rname = true.reference_name
				rname = truth_rname_mappings.get(rname, rname)
				expected_ref_id = tested_ref_ids.get(rname)
				if expected_ref_id is None:
					if true.reference_name not in seen_non_matching_reference_names:
						print(f"WARNING: Reference name {rname} not found in tested set. (Further warnings will be suppressed.)", file = sys.stderr)
						seen_non_matching_reference_names.add(true.reference_name)

					continue

				total_reads += 1
				true_pos = alignment_position(true)

				found_true_positive = False
				found_false_negative = False
				for tested_aln in tested:
					total_alns += 1

					if not tested_aln.is_mapped:
						found_false_negative = True
						continue

					# Compare the positions.
					if true.reference_id != tested_aln.reference_id:
						false_positives += 1
						continue

					tested_pos = alignment_position(tested_aln)
					distance = abs(true_pos - tested_pos)
					if distance_threshold < distance:
						false_positives += 1
						continue

					# Found a match. We do not care if the alignment is secondary or supplementary.
					found_true_positive = True

				if found_true_positive:
					true_positives += 1
					false_positives += int(found_false_negative)
				else:
					# No true positives found; classify as false negative and also count the wrong alignments as false positives.
					false_negatives += 1
	
	precision = fp_div(float(true_positives), true_positives + false_positives)
	recall = fp_div(float(true_positives), true_positives + false_negatives)
	f1_score = fp_div(2.0 * true_positives, 2.0 * true_positives + false_positives + false_negatives)
	return total_reads, total_alns, precision, recall, f1_score


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Calculate the precision, recall, and F1 score from the given alignments by comparing the aligned positions with the given distance threshold.")
	parser.add_argument('-d', '--distance-threshold', type = int, default = 0, metavar = "DISTANCE", help = "Distance threshold")
	parser.add_argument('truth', type = str, help = "True alignments as SAM or BAM, sorted by QNAME")
	parser.add_argument('tested', type = str, help = "Tested alignments as SAM or BAM, sorted by QNAME")
	parser.add_argument('-o', '--omit-header', action = "store_true", help = "Omit header from output")
	parser.add_argument('--truth-rname-mapping', type = str, metavar = "PATH", help = "Map RNAMEs in truth as specified in the given TSV file")
	parser.add_argument('--report-missing-in-tested', action = 'store_true', help = "Report alignments missing from tested to stderr.")
	args = parser.parse_args()
	total_reads, total_alns, precision, recall, f1_score = calculate_precision_and_recall(args.truth, args.tested, args.distance_threshold, args.truth_rname_mapping, args.report_missing_in_tested)
	if not args.omit_header:
		print("TOTAL_READS\tTOTAL_ALIGNMENTS\tPRECISION\tRECALL\tF1_SCORE")
	print(f"{total_reads}\t{total_alignments}\t{precision}\t{recall}\t{f1_score}")
