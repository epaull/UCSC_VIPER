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
from scipy import stats
from collections import defaultdict
parser = OptionParser()
parser.add_option("-x","--subtypes", dest="subtypes",action="store",default=None,help="Subtype definitions")
parser.add_option("-s","--drug_events",dest="events",action="store",default=None,help="VIPER scores for phospho regulators")
parser.add_option("-t","--drug_activity",dest="activities",action="store",default=None,help="VIPER scores for expression regulators")
parser.add_option("-o","--output",dest="output",action="store",default=None,help="Full network .sif file for inferring TF regulons")
(opts, args) = parser.parse_args()

from tiedie_util import *
from collections import defaultdict
from pathway import Pathway, BasicPathValidator 
import operator

def parseSubtypes(file, setA, setB):
	"""
	setA: types of classes i.e. MAPK or AKT_PI3K 
	"""
	samples_A = set()
	samples_B = set()

	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		drug_name = parts[0]
		drug_class = parts[10]
	
		if drug_class in setA:
			for conc in ['0.04', '0.12', '0.37', '1.11', '3.33', '10']:
				samples_A.add(drug_name+'_'+conc)

		if drug_class in setB:
			for conc in ['0.04', '0.12', '0.37', '1.11', '3.33', '10']:
				samples_B.add(drug_name+'_'+conc)

	return (samples_A, samples_B)

def refineEdges(edges, sources, targets, max_depth):
	## do a djikstra search, setting edge weights to the differential score

	max_k = 10
	g = nx.DiGraph()
	for (s,i,t) in edges:
		edge_cost = 1 - abs(edges[(s,i,t)])
		g.add_edge(s,t,cost=edge_cost,i=i)
		# hack semi-directional graph with extra direction for undirected edges
		if i == 'undirected' or i == 'PPI>':
			g.add_edge(t,s,cost=edge_cost,i=i)

	# this is one way to do it: the other is to weight each path by the 
	# least-specific edge it contains. 
	#paths = nx.all_pairs_dijkstra_path(g, weight='cost')
	# for all source/target pairs, get shortest paths
	# weighted by differential score
	filtered_edges = set()
	for source in sources:
		for target in targets:
			try:
				#for path in nx.all_shortest_paths(g, source, target, weight='cost'):
				scored_paths = {}
				for path in nx.all_simple_paths(g, source, target, cutoff=max_depth):
					max_cost = 0.0
					for i in range(0, len(path)-1):
						this_source = path[i]
						this_target = path[i+1]
						type = g[this_source][this_target]['i']
						c = g[this_source][this_target]['cost']
						
						if float(c) > max_cost:
							max_cost = float(c)

					scored_paths['_'.join(path)] = max_cost		

				k = 1
				for (path, path_cost) in sorted(scored_paths.items(), key=operator.itemgetter(1)):
					if k > max_k:
						break

					path = path.split('_')
					for i in range(0, len(path)-1):
						this_source = path[i]
						this_target = path[i+1]
						type = g[this_source][this_target]['i']
						filtered_edges.add( (this_source, type, this_target) )

					k += 1

			except:
				#traceback.print_exc(file=sys.stdout)
				continue

	return filtered_edges

def doTtest(data, list1, list2):

	vec_1 = []
	vec_2 = []

	for sample in data:
		if sample in list1:
			vec_1.append(data[sample])
		if sample in list2:
			vec_2.append(data[sample])


	t, prob = stats.ttest_ind(vec_1, vec_2)
	return (t, prob)

# data 
events = parseMatrix(opts.events, None, 0.0000001)
activities = parseMatrix(opts.activities, None, 0, transpose=True)

## filter VIPER scores by quantile
#filtered_activities = {}
#for gene in activities:
#	l = []
#	for sample in activities[sample]:

# split experiment samples/networks in two categories
samples_A, samples_B = parseSubtypes(opts.subtypes, set(['AKT_PI3K', 'MTOR']), set(['MAPK']))

for viper_tf in activities.keys():
	t, p_val = doTtest(activities[viper_tf], samples_A, samples_B)
	print viper_tf+'\t'+str(p_val)
