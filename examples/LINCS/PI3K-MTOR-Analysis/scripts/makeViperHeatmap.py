#!/usr/bin/env python

import traceback
import os, sys, re
from optparse import OptionParser
import numpy as np
from scipy import stats
from collections import defaultdict
parser = OptionParser()
parser.add_option("-x","--subtypes", dest="subtypes",action="store",default=None,help="Subtype definitions")
parser.add_option("-t","--drug_activity",dest="activities",action="store",default=None,help="VIPER scores for expression regulators")
(opts, args) = parser.parse_args()

from tiedie_util import *
from collections import defaultdict
from pathway import Pathway, BasicPathValidator 
import operator

def parseSubtypes(file, setA, setB, covered_drugs):
	"""
	setA: types of classes i.e. MAPK or AKT_PI3K 
	"""
	drugs_A = {}
	drug_classs_B = {}

	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		drug_name = parts[0]
		protein_target = parts[4]
		drug_class = parts[10]
	
		if drug_class in setA:
			if protein_target not in drugs_A:
				drugs_A[protein_target] = set()
			for conc in ['0.04', '0.12', '0.37', '1.11', '3.33', '10']:
				key = drug_name+'_'+conc
				if key not in covered_drugs:
					continue
				drugs_A[protein_target].add( key )

		if drug_class in setB:
			if protein_target not in drug_classs_B:
				drug_classs_B[protein_target] = set()
			for conc in ['0.04', '0.12', '0.37', '1.11', '3.33', '10']:
				key = drug_name+'_'+conc
				if key not in covered_drugs:
					continue
				drug_classs_B[protein_target].add( key )

	return (drugs_A, drug_classs_B)

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

def doTtest(vec1, vec2):
	t, prob = stats.ttest_ind(vec1, vec2)
	return t

# data indexed by gene, then by drug 
activities = parseMatrix(opts.activities, None, 0, transpose=True)
# all drugs with activity data
drug_conc_covered = set()
for gene in activities:
	for drug in activities[gene]:
		drug_conc_covered.add(drug)
	# it's a square matrix, just use the first row...
	break

# split experiment drug_classs/networks in two categories
drugs_A, drugs_B = parseSubtypes(opts.subtypes, set(['AKT_PI3K', 'MTOR']), set(['MAPK']), drug_conc_covered)

# calculate differential scores between genes/regulators in these subclasses
diff_scores = {}
for gene in activities:

	class_A_scores = []
	class_B_scores = []

	for drug_class in drugs_A:
		for drug in drugs_A[drug_class]:
			if drug in activities[gene]:
				class_A_scores.append(float(activities[gene][drug]))

	for drug_class in drugs_B:
		for drug in drugs_B[drug_class]:
			if drug in activities[gene]:
				class_B_scores.append(float(activities[gene][drug]))

	t = doTtest(class_A_scores, class_B_scores)
	#sys.stderr.write(gene+'\t'+str(t)+'\n')
	diff_scores[gene] = t

# sort by t-statistic between the 2 major drug classes
gene_order = [gene for (gene, score) in sorted(diff_scores.items(), key=operator.itemgetter(1), reverse=True)]


# order all perturbation inferences by drug class
drugs_order = []
# corresponding drug classes for these perturbations
drug_class_labels = []

for drug_class in drugs_A:
	if drug_class.startswith("PIK3"):
		for drug in drugs_A[drug_class]:
			drugs_order.append(drug)
			drug_class_labels.append('PIK3')
for drug_class in drugs_A:
	if drug_class.startswith("AKT"):
		for drug in drugs_A[drug_class]:
			drugs_order.append(drug)
			drug_class_labels.append('AKT')
for drug_class in drugs_A:
	if drug_class.startswith("MTOR"):
		for drug in drugs_A[drug_class]:
			drugs_order.append(drug)
			drug_class_labels.append('MTOR')
for drug_class in drugs_B:
	for drug in drugs_B[drug_class]:
		drugs_order.append(drug)
		drug_class_labels.append('MAPK')

# print out matrix in this order

# set 
print 'Key\tType\t'+'\t'.join(drug_class_labels)
print 'Class\tType_Color\t'+'\t'.join(drugs_order)
for gene in gene_order:
	printstr = gene+'\t'+gene
	for drug in drugs_order:
		printstr += '\t'+str(activities[gene][drug])
	print printstr
