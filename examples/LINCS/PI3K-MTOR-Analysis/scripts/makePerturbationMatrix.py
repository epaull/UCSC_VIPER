#!/usr/bin/env	python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-t","--ljp_targets",dest="ljp_targets",action="store",type="string")
(opts, args) = parser.parse_args()

import sys

header = True
descriptions = {}
for line in open(opts.ljp_targets, 'r'):
	if header:
		header = False
		continue
	name, pathway_role, pathway, drug_class = line.rstrip().split("\t")
	#descriptions[name] = (pathway_role, pathway, drug_class)
	descriptions[name] = pathway_role

## parse gene-expansions: map pathway_role to gene sets of interest
expansion_map = {}
for arg in args:

	map_from = set()
	for line in open(arg, 'r'):
		parts = line.rstrip().split('\t')
		if len(parts) == 1:
			map_from.add(parts[0])
			expansion_map[parts[0]] = set()
		else:
			for gene in map_from:
				expansion_map[gene].add(parts[0])
		
##
## Take in a list of drugs + concentrations, separated by a ':', and assign the 
## gene targets
##

all_targets = set()
# for each drug/conc id save the gene target
perturb_data = {}
for line in sys.stdin:
	drug, concentration = line.rstrip().split('_')
	target = descriptions[drug]
	if target in expansion_map:
		perturb_data[line.rstrip()] = set()
		for t in expansion_map[target]:
			all_targets.add(t)
			perturb_data[line.rstrip()].add(t)
	else:
		all_targets.add(target)
		perturb_data[line.rstrip()] = set([target])

targets = list(all_targets)	
print 'Drug\t'+'\t'.join(targets)
for id in perturb_data:
	printstr = id
	for gene in targets:
		if gene in perturb_data[id]:
			printstr += '\t1'
		else:
			printstr += '\t0'

	print printstr


