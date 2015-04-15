#!/usr/bin/env	python

import sys

edges = {}
for line in sys.stdin:
	a, i, b = line.rstrip().split("\t")	
	if (a,b) in edges or (b,a) in edges:
		continue
	print line.rstrip()
	edges[(a,b)] = 1
