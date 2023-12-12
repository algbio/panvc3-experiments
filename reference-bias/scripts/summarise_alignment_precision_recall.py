# 
# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

# Summarise alignment precision and recall, output as TSV.

import os
import re
import sys


# E.g. HG001.panvc3-bowtie2-f14-d25.max-mapq.d5.tsv
file_name_pattern = re.compile(
	r"^ (?: .* [/])? (?P<sample> [^.]+) [.] (?P<wf> .+) [.] d(?P<dist> \d+) [.]tsv $",
	re.VERBOSE
)


if __name__ == "__main__":
	print("SAMPLE\tWORKFLOW\tDISTANCE\tTOTAL_READS\tTOTAL_ALIGNMENTS\tPRECISION\tRECALL\tF1_SCORE")
	for ds in os.scandir("alignment-precision-recall"):
		mm = file_name_pattern.match(ds.path)
		if mm is None:
			print(f"Unable to match {ds.path}.\n", file = sys.stderr)
			sys.exit(1)
		sample = mm.group("sample")
		wf = mm.group("wf")
		dist = int(mm.group("dist"))

		with open(ds.path) as ff:
			is_first = True
			for line in ff:
				if is_first:
					is_first = False
					continue
				sys.stdout.write(f"{sample}\t{wf}\t{dist}\t")
				sys.stdout.write(line)
