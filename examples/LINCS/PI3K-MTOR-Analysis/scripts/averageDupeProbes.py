#!/usr/bin/env	python


import sys
from collections import defaultdict

# indexed by sample, then by gene
data = {}
all_genes = set()
header = None
for line in sys.stdin:
	parts = line.rstrip().split('\t')
	if not header:
		header = parts
		for i in range(1, len(header)):
			data[header[i]] = defaultdict(set)
		continue


	# add all values for any gene matching this label
	gene = parts[0]
	for i in range(1, len(parts)):
		data[header[i]][gene].add(float(parts[i]))

	all_genes.add(gene)

print '\t'.join(header)
for gene in all_genes:

	printstr = gene

	# average over all values with this gene-label 
	for sample in header[1:len(header)]:
		count = 0.0
		sum = 0.0
		for val in data[sample][gene]:
			sum += val	
			count += 1.0

		mean = sum/count
		printstr += '\t'+str(mean)

	print printstr
