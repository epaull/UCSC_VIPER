#!/usr/bin/env python

###
### mapSamples: Script to map individual drug data onto a network, yeilding a sample-specific
### network
###
### Authors:
###
###		Evan Paull (epaull@soe.ucsc.edu)
###
###
### Minimum Data Inputs: 
###		
###		- drug_matrix: a gene expression matrix with normalized, median centered values per sample (column), with
###		gene names on the rows. Typically the values should be centered by the expression of a group of normal-adjacent 
###		drugs.
###		- a search pathway in .sif format (geneA <interaction> geneB). Most likely a TieDIE global solution network
###
###

import traceback
import os, sys, re
from optparse import OptionParser
import numpy as np
from collections import defaultdict
parser = OptionParser()
parser.add_option("-n","--network",dest="network",action="store",default=None,help="Full network .sif file for directionality")
parser.add_option("-c","--classes",dest="eq_classes",action="store",default=None)
(opts, args) = parser.parse_args()

def parseEdges(file):

	"""
	Get equivalence classes mapping class name to set of nodes
	"""
	edges = {}
	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		s, i, t = parts[1], parts[2], parts[3]
		
		edges[(s,i,t)] = (parts[0], parts[4], parts[5], parts[6])

	return edges

def parseClasses(file):

	"""
	Get equivalence classes mapping class name to set of nodes
	"""
	classes = defaultdict(set)
	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		classes[parts[0]] = parts[1:len(parts)]

	return classes


edges = parseEdges(opts.network)
# 'PI3K':set(PIK3CA, PIK3RA...)
classes = parseClasses(opts.eq_classes)


processed_edges =  {}
for (s, i, t) in edges: 

	data = edges[(s,i,t)]
	source_identifier = s
	target_identifier = t
	# replace node label with the class name if it's in any...
	for eq_class in classes:
		if s in classes[eq_class]:
			source_identifier = eq_class
		if t in classes[eq_class]:
			target_identifier = eq_class

	processed_edges[(source_identifier, i, target_identifier)] = data

for (s, i, t) in processed_edges:
	data = processed_edges[(s,i,t)]	
	print '\t'.join([data[0], s, i, t, data[1], data[2], data[3]])

