#!/usr/bin/env python

"""
Take a VIPER-inferred matrix of protein activities, median center on DMSO/Control samples
,normalize and print out just the non-control samples.
"""

import sys
from numpy import *

perturbation_idx = []
control_idx = []

header = False
for line in sys.stdin:
	
	parts = line.rstrip().split('\t')
	if not header:
		header = parts
		for i in range(1, len(parts)):
			if parts[i].startswith('DMSO'):
				control_idx.append(i)
			else:
				perturbation_idx.append(i)
		print 'Key\t'+'\t'.join([header[i] for i in perturbation_idx])
		continue


	printstr = parts[0]

	control_values = []
	for i in control_idx:
		control_values.append(float(parts[i]))

	for i in perturbation_idx:
		normalized_viper = float(float(parts[i]) - float(median(control_values))) / float(std(control_values))
		printstr += '\t'+str(normalized_viper)

	print printstr
