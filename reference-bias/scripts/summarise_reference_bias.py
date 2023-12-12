# 
# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

# Summarise reference bias, output as TSV.

import gzip
import os
import re
import sys


# E.g. HG001.panvc3-bowtie2-f14-d25.mapq-recalculated.1.all.mc10.txt
file_name_pattern = re.compile(
	r"^ (?: .* [/])? (?P<sample> [^.]+) [.] (?P<wf> .+) [.] (?P<chromosome> \d+) [.] (?P<regions> [^.]+) [.] mc (?P<min_cov> \d+) [.]txt[.]gz $",
	re.VERBOSE
)


if __name__ == "__main__":
	print("SAMPLE\tWORKFLOW\tCHROMOSOME\tREGIONS\tMIN_COVERAGE\tBALANCE\tREF_LENGTH\tALT_LENGTH")
	for ds in os.scandir("reference-bias"):
		mm = file_name_pattern.match(ds.path)
		if mm is None:
			print(f"Unable to match {ds.path}.\n", file = sys.stderr)
			sys.exit(1)
		sample = mm.group("sample")
		wf = mm.group("wf")
		chromosome = mm.group("chromosome")
		regions = mm.group("regions")
		min_cov = int(mm.group("min_cov"))

		with gzip.open(ds.path, "rt") as ff:
			is_first = True
			for line in ff:
				if is_first:
					is_first = False
					continue
				sys.stdout.write(f"{sample}\t{wf}\t{chromosome}\t{regions}\t{min_cov}\t")
				sys.stdout.write(line)
