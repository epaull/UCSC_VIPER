#!/usr/bin/env	python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()

import sys
import operator
from collections import defaultdict

# Parse drug metadata and index by drug name
fields = set(['DrugName', 'TargetName', 'Role', 'Pathway', 'DrugClass'])
field_map = {}
drug_data = {}
header = False
for line in open(args[0], 'r'):
	
	parts = line.rstrip().split("\t")
	if not header:
		header = True
		# save the index for each data field
		for i in range(0, len(parts)):
			key = parts[i]
			if key in fields:
				field_map[i] = key
		continue

	attributes = {}
	for i in range(0, len(parts)):
		if i in field_map:
			key = field_map[i]
		else:
			continue

		val = parts[i]
		attributes[key] = val

	drug_data[attributes['DrugName']] = attributes	


# parse a data matrix, the first line is the header with drug names
header = False

# build an index order
sorted_index_order = []

drug_data['DMSO'] = {'DrugClass':'Control'}

drug_class_order = ['Control', 'ERBB', 'MAPK', 'MEK', 'NRTK', 'RTK', 'PI3K/MTOR', 'AKT']

for line in sys.stdin:
	
	parts = line.rstrip().split('\t')
	if not header:
		header = parts
		index_toKey_map = {}
		key_toIndex_map = defaultdict(set)
		# find the drug class each drug specified in the header belongs to (don't include if it doesn't have one)
		for i in range(1, len(parts)):
			drug_name = parts[i].split('_')[0]
			if drug_name in drug_data:
				metadata = drug_data[drug_name]
				key = metadata['DrugClass']
				key_toIndex_map[key].add(i)

		# now build the sorted index using the pre-defined order of drug classes
		for drug_class in drug_class_order:
			if drug_class not in key_toIndex_map:
				continue
			for index in key_toIndex_map[drug_class]:
				sorted_index_order.append(index)

		print 'Name\tDesc\t'+'\t'.join(parts[i] for i in sorted_index_order)
		continue

	# CDT format
	printstr = parts[0]+'\t'+parts[0]
	for index in sorted_index_order:
		printstr += '\t'+parts[index]
	print printstr
	
