#!/usr/bin/env	python

import sys

regulon = {}
for line in sys.stdin:
	tf, score, target = line.rstrip().split("\t")

	if tf not in regulon:
		regulon[tf] = []

	regulon[tf].append( (target, score) )


for tf in regulon:

	printstr = tf
	for (target, score) in regulon[tf]:
		printstr += target+'\t'+score

	print printstr	
