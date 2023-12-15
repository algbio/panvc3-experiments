# 
# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

# Summarise Truvari’s results, output as TSV.

import json
import os
import re
import sys

file_name_pattern = re.compile(
    # {sample_id}.{wf}.{variant_caller}.{vc_regions}.{eval_regions}.{filter}
	# E.g. NA24385.panvc3-bowtie2-f25-d50.max-mapq.manta.all.chr1_selected.all
	r"^ (?P<sample> [^.]+) [.] (?P<wf> .+) [.] (?P<vc> [^.]+) [.] (?P<vc_regions> [^.]+) [.] (?P<eval_regions> [^.]+) [.] (?P<filters> [^.]+) $",
	re.VERBOSE
)


if __name__ == "__main__":
	print("SAMPLE\tWORKFLOW\tVARIANT_CALLER\tVC_REGIONS\tEVALUATED_REGIONS\tFILTERS\tCOUNT\tPRECISION\tRECALL")
	for ds in os.scandir("truvari"):
		mm = file_name_pattern.match(ds.name)
		if mm is None:
			print(f"Unable to match {ds.name}.\n", file = sys.stderr)
			sys.exit(1)
		sample = mm.group("sample")
		workflow = mm.group("wf")
		variant_caller = mm.group("vc")
		vc_regions = mm.group("vc_regions")
		eval_regions = mm.group("eval_regions")
		filters = mm.group("filters")
		with open(f"{ds.path}/summary.json") as ff:
			summary = json.load(ff)
			print(f"{sample}\t{workflow}\t{variant_caller}\t{vc_regions}\t{eval_regions}\t{filters}\t{summary['comp cnt']}\t{summary['precision']}\t{summary['recall']}")
