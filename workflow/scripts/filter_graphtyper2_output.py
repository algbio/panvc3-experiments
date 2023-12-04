# Copyright (c) 2021-2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import sys
import vcf.filters

class AAScoreFilter(vcf.filters.Base):
	'Filter sites by AAScore'

	name = 'AAScore'

	@classmethod
	def customize_parser(self, parser):
		parser.add_argument('--aa-score', type = float, default = 0.5, help = "Filter sites below this value")

	def __init__(self, args):
		self.threshold = args.aa_score

	def __call__(self, record):
		aa_scores = record.INFO["AAScore"]
		if 1 < len(aa_scores):
			sys.stderr.write("Got multiple ALT values:\n")
			sys.stderr.write(f"{record}\n")
			return 0

		score = aa_scores[0]
		if score < self.threshold:
			return score

		return None
