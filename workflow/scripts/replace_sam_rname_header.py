# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import sys


def parse_replacements(path):
	def helper():
		with open(path, "r") as fp:
			for line in fp:
				if line.startswith("#"):
					continue
				line = line.rstrip("\n")
				yield line.split("\t")
	return {k: v for k, v in helper()}


def modify(field, replacements):
	if not field.startswith("SN:"):
		return field
	
	rname = field.lstrip("SN:")
	new_rname = replacements.get(rname, rname)
	return f"SN:{new_rname}"


if __name__ == "__main__":
	parser = argparse.ArgumentParser("Rewrite the reference names in SAM header")
	parser.add_argument('replacements', type = str, help = "Path to replacements as TSV")
	args = parser.parse_args()

	replacements = parse_replacements(args.replacements)

	for line in sys.stdin:
		if not line.startswith("@SQ\t"):
			sys.stdout.write(line)
			continue

		line = line.rstrip("\n")
		fields = line.split("\t")
		new_fields = [modify(field, replacements) for field in fields]
		sys.stdout.write("\t".join(new_fields))
		sys.stdout.write("\n")
