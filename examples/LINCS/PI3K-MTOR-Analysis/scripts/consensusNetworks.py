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
parser.add_option("-x","--subtypes", dest="subtypes",action="store",default=None,help="Subtype definitions")
parser.add_option("-s","--drug_events",dest="events",action="store",default=None,help="VIPER scores for phospho regulators")
parser.add_option("-t","--drug_activity",dest="activities",action="store",default=None,help="VIPER scores for expression regulators")
parser.add_option("-n","--network",dest="network",action="store",default=None,help="Full network .sif file for directionality")
parser.add_option("-d","--directory",dest="directory",action="store",default=None,help="Directory with individual network .sif files")
parser.add_option("-o","--output",dest="output",action="store",default=None,help="Full network .sif file for inferring TF regulons")
parser.add_option("-z","--depth",dest="depth",action="store",default=3)
parser.add_option("-c","--classes",dest="eq_classes",action="store",default=None)
(opts, args) = parser.parse_args()

from tiedie_util import *
from collections import defaultdict
from pathway import Pathway, BasicPathValidator 
import operator

def parseNodeClasses(file):

	"""
	Get equivalence classes mapping class name to set of nodes
	"""
	# map each node to it's node-class/abstract (often the same node name)
	classes = {}
	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		for i in range(1, len(parts)):
			classes[parts[i]] = parts[0]

	return classes

# build an index, source to targets fro the directed graph
def parseNet(network, node_classes, header=False):
	"""
	Build a directed network from a .sif file. 
	
	Inputs:
		A network in .sif format, tab-separated (<source> <interaction> <target>)

	Returns
		A network in hash key format, i.e. convert two lines of a file:
			<source>	<interaction1>	<target1>
			<source>	<interaction2>	<target2>
		To:	
			{'source': set( (interaction, target1), (interaction, target2) )
	"""
	net = {}
	lineno = 0
	for line in open(network, 'r'):

		parts = line.rstrip().split("\t")
		source = parts[0]
		interaction = parts[1]
		target = parts[2]

		# map to node class abstractions, if available
		if source in node_classes:
			source = node_classes[source]
		if target in node_classes:
			target = node_classes[target]

		if header and lineno == 0:
			lineno += 1
			continue

		if source not in net:
			net[source] = set()

		net[source].add((interaction, target))

	return net


def getEdges(net):

	edgelist = {}
	for s in net:
		for (i,t) in net[s]:
			# ignore self-links
			if s == t:
				continue

			edgelist[(s,t)] = i

	return edgelist
	
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

# data 
events = parseMatrix(opts.events, None, 0.0000001)
activities = parseMatrix(opts.activities, None, 0, transpose=True)

## filter VIPER scores by quantile
#filtered_activities = {}
#for gene in activities:
#	l = []
#	for sample in activities[sample]:

# get node mappings to node-classes to summarize the network
node_classes = parseNodeClasses(opts.eq_classes)
scaffold_network = parseNet(opts.network, node_classes)
scaffold_edges = getEdges(scaffold_network)

# split experiment samples/networks in two categories
samples_A, samples_B = parseSubtypes(opts.subtypes, set(['AKT_PI3K', 'MTOR']), set(['MAPK']))

directed_drug_networks = {}
edge_counts = {}
edge_counts['A'] = {}
edge_counts['B'] = {}
edge_counts['total'] = defaultdict(int)

drug_networks_nodes = {}
# these are the samples to evaluate
test_drugs = set(events.keys()).intersection(activities.keys())
samples_A = samples_A.intersection(test_drugs)
samples_B = samples_B.intersection(test_drugs)

# file names will correspond to sample names: the drug plus the concentration...
# these actually have networks
active_networks = set()
for network in os.listdir(opts.directory):
	if network.endswith('activities.txt'):
		continue

	drug_network = parseNet(opts.directory+'/'+network, node_classes, header=True)

	all_nodes = set()

	drug_name = network.rstrip('.txt')
	if drug_name not in test_drugs:
		continue

	drug_class = None
	if drug_name in samples_A:
		drug_class = 'A'
	elif drug_name in samples_B:
		drug_class = 'B'
	else:
		continue


	active_networks.add(drug_name)

	# FIXME: merge/add directionality data
	directed_net = {}
	for s in drug_network:
		for (i, t) in drug_network[s]:
			directed_interaction = None
			if (s,t) in scaffold_edges:
				if s not in directed_net:
					directed_net[s] = set()
				directed_net[s].add( (scaffold_edges[(s,t)],t) )
			elif (t,s) in scaffold_edges:
				if t not in directed_net:
					directed_net[t] = set()
				directed_net[t].add( (scaffold_edges[(t,s)],s) )
			else:
				# otherwise add the HPRD interaction
				if s not in directed_net:
					directed_net[s] = set()
				directed_net[s].add( (i,t) )


	
	# save for this drug
	directed_drug_networks[drug_name] = directed_net
	for source in directed_drug_networks[drug_name]:
		all_nodes.add(source)
		for (i, t) in directed_drug_networks[drug_name][source]:
			all_nodes.add(t)
			edge = (source, i, t)
			if edge not in edge_counts[drug_class]:
				edge_counts[drug_class][edge] = 0
			edge_counts[drug_class][edge] += 1
			edge_counts['total'][edge] += 1

	
	drug_networks_nodes[drug_name] = all_nodes

##
## Remove any samples without networks...

samples_A = samples_A.intersection(active_networks)
samples_B = samples_B.intersection(active_networks)

# indexed by edges: store statistics for
edge_score = {}
# indexed by node
node_score = {}

# tally edge counts for each 
for edge in edge_counts['total']:
	# the overall frequency of this edge, over both subtypes
	edge_score[edge] = {}
	# initialize node scores
	node_score[edge[0]] = defaultdict(int)
	node_score[edge[2]] = defaultdict(int)

	edge_score[edge]['total'] = edge_counts['total'][edge]/float(len(samples_A)+len(samples_B))
	# 
	if edge in edge_counts['A']:
		# add edge scores
		edge_score[edge]['A'] = edge_counts['A'][edge]/float(len(samples_A))
	else:
		edge_score[edge]['A'] = 0.0
		
	if edge in edge_counts['B']:
		# edges...
		edge_score[edge]['B'] = edge_counts['B'][edge]/float(len(samples_B))
	else:
		edge_score[edge]['B'] = 0.0

	edge_score[edge]['average'] = (edge_score[edge]['A'] + edge_score[edge]['B'])/2.0
	edge_score[edge]['diff'] = edge_score[edge]['A'] - edge_score[edge]['B']
	

# get node statistics
for node in node_score:
	count_A = 0
	for sample in samples_A:
		if node in drug_networks_nodes[sample]:
			count_A += 1
	count_B = 0
	for sample in samples_B:
		if node in drug_networks_nodes[sample]:
			count_B += 1

	node_score[node]['A'] = count_A/float(len(samples_A))
	node_score[node]['B'] = count_B/float(len(samples_B))
	node_score[node]['average'] = (node_score[node]['A'] + node_score[node]['B'])/2.0

fh = open(opts.output+'/consensus-networks.txt', 'w')
for edge in edge_score:
	fh.write('\t'.join(edge)+'\t'+str(edge_score[edge]['diff'])+'\t'+str(edge_score[edge]['average'])+'\n')
fh.close()

fh = open(opts.output+'/consensus-networks.nodes.txt', 'w')
for node in node_score:
	fh.write(node+'\t'+str(node_score[node]['average'])+'\n')
fh.close()

#for node in node_counts_egfr:
#	print node+'\t'+str(node_counts_egfr[node])
#FIXME: get these from the input matrix files, not the summary files
#mek_nodes = set()
#for line in open('../DATA/mek.upstream.txt', 'r'):
#	mek_nodes.add(line.split('\t')[0])
#
#pi3k_nodes = set()
#for line in open('../DATA/pi3k.upstream.txt', 'r'):
#	pi3k_nodes.add(line.split('\t')[0])
#
#tf_nodes = set()
#for line in open('tiedie/input/viperScores.30TF_whitelist.txt', 'r'):
#	tf_nodes.add(line.split('\t')[0])
#
### do a djikstra search, setting edge weights to the differential score
#mek_summary = refineEdges(mek_edges, mek_nodes, tf_nodes, int(opts.depth))
#for edge in mek_summary:
#	print 'MEK\t'+'\t'.join(edge)
#
#pi3k_summary = refineEdges(pi3k_edges, pi3k_nodes, tf_nodes, int(opts.depth))
#for edge in pi3k_summary:
#	print 'PI3K\t'+'\t'.join(edge)
