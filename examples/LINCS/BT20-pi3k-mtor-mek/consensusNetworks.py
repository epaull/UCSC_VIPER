#!/usr/bin/env	python

#
# Take the directory with networks for each, also parse the subtypes. Generate a consensus network for each
# by counting the fractional occurance of each edge. Diff the final consensus networks: 
#
# Also diff viper scores, above the 1.5 threshold. Look at differential regulators at the bottom. 

from optparse import OptionParser
import os
from collections import defaultdict

parser = OptionParser()
parser.add_option("-d","--directory",dest="directory",action="store",default=None)
parser.add_option("-s","--subtypes",dest="subtypes",action="store",default=None)
(options, args) = parser.parse_args()

def parseNet(network):

	net = {}
	for line in open(network, 'r'):
		parts = line.rstrip().split("\t")
		source = None
		target = None
		interaction = None
		if len(parts) > 2:
			source = parts[0]
			interaction = parts[1]
			target = parts[2]
		else:
			source = parts[0]
			target = parts[1]
		
		if source == 'source':
			continue		
		net[(source, target)] = interaction

	return net

def parseSubtypes(file):
	data = defaultdict(set)
	for line in open(file, 'r'):
		drug, subtype = line.rstrip().split("\t")
		data[subtype].add(drug)

	return data

networks = {}
subtypes = parseSubtypes(options.subtypes)

edge_counts = {}
all_edges = {}

for subtype in subtypes:
	edge_counts[subtype] = defaultdict(int)

for network in os.listdir(options.directory):
	if network.endswith('activities.txt'):
		continue
	netObj = parseNet(options.directory+'/'+network)
	
	drug_name = network.split('.')[0]
	subtype = None
	for s in subtypes:
		if drug_name in subtypes[s]:
			subtype = s
			break

	networks[network] = netObj
	for edge in netObj:
		edge_counts[subtype][edge] += 1	
		# save the interaction as well
		all_edges[edge] = netObj[edge]


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



