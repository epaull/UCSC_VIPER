#!/usr/bin/env python

###
### mapSamples: Script to map individual sample data onto a network, yeilding a sample-specific
### network
###
### Authors:
###
###		Evan Paull (epaull@soe.ucsc.edu)
###
###
### Minimum Data Inputs: 
###		
###		- sample_matrix: a gene expression matrix with normalized, median centered values per sample (column), with
###		gene names on the rows. Typically the values should be centered by the expression of a group of normal-adjacent 
###		samples.
###		- a search pathway in .sif format (geneA <interaction> geneB). Most likely a TieDIE global solution network
###
###

import os, sys, re
from optparse import OptionParser
import numpy as np
from collections import defaultdict
parser = OptionParser()
parser.add_option("-x","--subtypes", dest="subtypes",action="store",default=None,help="Subtype definitions")
parser.add_option("-s","--sample_events",dest="events",action="store",default=None,help="VIPER scores for phospho regulators")
parser.add_option("-t","--sample_activity",dest="activities",action="store",default=None,help="VIPER scores for expression regulators")
parser.add_option("-n","--network",dest="network",action="store",default=None,help="Full network .sif file for directionality")
parser.add_option("-d","--directory",dest="directory",action="store",default=None,help="Directory with individual network .sif files")
parser.add_option("-o","--output",dest="output",action="store",default=None,help="Full network .sif file for inferring TF regulons")
(opts, args) = parser.parse_args()

from tiedie_util import *
from collections import defaultdict
from pathway import Pathway, BasicPathValidator 

# data 
events = parseMatrix(opts.events, None, 0.0)
activities = parseMatrix(opts.activities, None, 0.0)
scaffold_network = parseNet(opts.network)

def parseSubtypes(file):
	data = defaultdict(set)
	for line in open(file, 'r'):
		drug, subtype = line.rstrip().split("\t")
		data[subtype].add(drug)

	return data

subtypes = parseSubtypes(options.subtypes)


networks = {}
# parse networks and save for each sample
for subtype in subtypes:
	edge_counts[subtype] = defaultdict(int)

for network in os.listdir(options.directory):
	if network.endswith('activities.txt'):
		continue
	sample_network = parseNet(options.directory+'/'+network)

	drug_name = network.split('.')[0]
	subtype = None
	for s in subtypes:
		if drug_name in subtypes[s]:
			subtype = s
			break

	# FIXME: merge/add directionality data

	# save for this sample
	networks[network] = netObj

test_samples = set(events.keys()).intersection(activities.keys())

# filter sample networks
all_edges = defaultdict(float)
sample_networks = {}

# score by subtype
edge_counts = {}
for sample in test_samples:

	sample_networks[sample] = set()
	print sample+'\tevents\t:'+','.join(events[sample])
	print sample+'\tactivities\t:'+','.join(activities[sample])

	# get upstream nodes
	upstream_nodes = set(events[sample].keys())
	downstream_nodes = activities[sample].keys()

	validator = BasicPathValidator()
	pathway = Pathway(networks[sample], validator=validator, opts={'undirected_edges':set(['PPI>'])})
	edges = pathway.allPaths(upstream_nodes, downstream_nodes, 3)



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
	print '\t'.join(all_subtypes)+'\t'+'\t'.join([edge[0], all_edges[edge], edge[1]])+'\t'+str(diff)




