# 
# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

# Summarise Truvari’s results, output as TSV.

import json
import os
import sys


if __name__ == "__main__":
	print("ALIGNER\tWORKFLOW\tVARIANT_CALLER\tVC_REGIONS\tEVALUATED_REGIONS\tFILTERS\tCOUNT\tPRECISION\tRECALL")
	for ds in os.scandir("truvari"):
		aligner, workflow, variant_caller, vc_regions, evaluated_regions, filters = ds.name.split(".")
		with open(f"{ds.path}/summary.json") as ff:
			summary = json.load(ff)
			print(f"{aligner}\t{workflow}\t{variant_caller}\t{vc_regions}\t{evaluated_regions}\t{filters}\t{summary['comp cnt']}\t{summary['precision']}\t{summary['recall']}")
