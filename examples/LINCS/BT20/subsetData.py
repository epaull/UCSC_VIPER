#!/usr/bin/env	python

import sys

gene_set = set()
for line in open("Cluster_BT20_T24_INFchdir.tsv", 'r'):
	parts = line.rstrip().split("\t")
	gene = parts[0]
	gene_set.add(gene)

	

for line in open("multinet.adj", 'r'):
	parts = line.rstrip().split('\t')
	tf = parts[0]
	if tf not in gene_set:
		continue
	printstr = tf
	for i in range(1, len(parts), 2):
		gene, mi = (parts[i], parts[i+1])
		if gene in gene_set:
			printstr += '\t'+gene+'\t'+mi 

	print printstr
	
