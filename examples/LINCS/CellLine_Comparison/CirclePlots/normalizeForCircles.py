#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()
import sys, os

import math


def parseMatrix(file, restrict_samples=None, binary_threshold=0.0, transpose=False):
	''' 
		Sample IDS should be the header line. Gene ids are the row names
		
		Input:
			binary_threshold: 'include data values only if they fall above this range (abs val)
			tf_parents: 

		Options:
			transpose: index by rows, then columns, instead of the default column/row spec
			
	'''


	# indexed by sample then by gene	
	data = {}
	 
	first = True
	sampleIDS = None
	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		row_id = parts[0]
		vals = parts[1:]
		if first:
			first = False
			column_ids = vals
			continue

		for i in range(0,len(vals)):
			val = None
			try:
				val = float(vals[i])
			except:
				continue
			column_id = column_ids[i]		

			if restrict_samples and column_id not in restrict_samples:
				continue
			if abs(val) < binary_threshold:
				continue	

			###
			### Get the gene expression, indexed by samples
			###
			if not transpose:
				if column_id not in data:
					data[column_id] = {}
				data[column_id][row_id] = val
			else:
				if row_id not in data:
					data[row_id] = {}
				data[row_id][column_id] = val

	return data

CONST = 1.0
data = parseMatrix(args[0], None, 0.0)

# find max value over all 
max_val = 0.0
genes = None
for sample in data.keys():
	genes = data[sample].keys()
	for key in data[sample]:
		v = abs(float(data[sample][key]))
		if v > max_val:
			max_val = v

sample_order = data.keys()
print 'Key\t'+'\t'.join(sample_order)
for gene in genes:
	printstr = gene
	for sample in sample_order:
		v = float(data[sample][gene])
		normalized = v/max_val
		normalized *= CONST
		normalized = math.pow(normalized, 2)

		printstr += '\t'+str(normalized)

	print printstr
