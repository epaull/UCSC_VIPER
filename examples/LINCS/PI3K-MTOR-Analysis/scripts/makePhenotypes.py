#!/usr/bin/env	python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()


header = True
descriptions = {}
for line in open(args[0], 'r'):
	if header:
		header = False
		continue
	name, pathway_role, pathway, drug_class = line.rstrip().split("\t")
	descriptions[name] = (pathway_role, pathway, drug_class)




index = 1
print '\t'.join(['Index', 'Drug', 'Dose', 'Role', 'Pathway', 'Drug Class'])
for line in open(args[1], 'r'):
	parts = line.rstrip().split('\t')
	for i in range(1, len(parts)):
		drug, dose = parts[i].split(':')

		role, pathway, drug_class = descriptions[drug]
		if drug_class == 'PI3K/MTOR' or drug_class == 'MEK':
			print '\t'.join([str(index), drug, dose, role, pathway, drug_class])
			index += 1

	break	
