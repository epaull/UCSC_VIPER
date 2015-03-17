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
(opts, args) = parser.parse_args()

from tiedie_util import *
from collections import defaultdict
from pathway import Pathway, BasicPathValidator 

def getEdges(net):

	edgelist = {}
	for s in net:
		for (i,t) in net[s]:
			# ignore self-links
			if s == t:
				continue

			edgelist[(s,t)] = i

	return edgelist
	
def parseSubtypes(file):
	data = defaultdict(set)
	for line in open(file, 'r'):
		drug, subtype = line.rstrip().split("\t")
		data[subtype].add(drug)

	return data

# data 
events = parseMatrix(opts.events, None, 0.0000001)
activities = parseMatrix(opts.activities, None, 1.5)
scaffold_network = parseNet(opts.network)
scaffold_edges = getEdges(scaffold_network)

subtypes = parseSubtypes(opts.subtypes)

directed_drug_networks = {}
edge_counts = {}
# parse networks and save for each drug
for subtype in subtypes:
	edge_counts[subtype] = defaultdict(int)

for network in os.listdir(opts.directory):
	if network.endswith('activities.txt'):
		continue

	drug_network = parseNet(opts.directory+'/'+network, header=True)

	drug_name = network.split('.')[0]
	subtype = None
	for s in subtypes:
		if drug_name in subtypes[s]:
			subtype = s
			break

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

	# save for this drug
	directed_drug_networks[drug_name] = directed_net

test_drugs = set(events.keys()).intersection(activities.keys())

# filter drug networks
all_edges = {}
drug_networks = {}

# score by subtype
edge_counts = {}
for subtype in subtypes:
	edge_counts[subtype] = defaultdict(int)

for drug in directed_drug_networks:

	# get upstream nodes
	#print "events:\t"+'\t'.join(events[drug].keys())
	#print "activities:\t"+'\t'.join(activities[drug].keys())

	upstream_nodes = set(events[drug].keys())
	downstream_nodes = activities[drug].keys()
	input_sets = {}
	input_sets['source'] = upstream_nodes
	input_sets['target'] = downstream_nodes

	if len(directed_drug_networks[drug]) == 0:
		print "No network for "+drug
		continue

	validator = BasicPathValidator(input_sets)
	pathway = Pathway(directed_drug_networks[drug], validator=validator, opts={'undirected_edges':set(['PPI>'])})
	edges = pathway.allPaths(upstream_nodes, downstream_nodes, 3)
	drug_networks[drug] = edges


	this_subtype = None
	for subtype in subtypes:
		if drug in subtypes[subtype]:
			this_subtype = subtype
			break	

	for edge in edges:
		edge_counts[this_subtype][edge] += 1

#
# 
#
score = {}
for subtype in subtypes:
	for edge in edge_counts[subtype]:
		if edge not in score:
			score[edge] = defaultdict(float)
		score[edge][subtype] = edge_counts[subtype][edge]/float(len(subtypes[subtype]))

all_subtypes = subtypes.keys()
for edge in score:
	diff = score[edge][all_subtypes[0]] - score[edge][all_subtypes[1]]
	print '\t'.join(all_subtypes)+'\t'+'\t'.join([edge[0], edge[2]])+'\t'+str(diff)

