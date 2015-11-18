#!/usr/bin/env	python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()


header = False
fields = set(['DrugName', 'TargetName', 'Role', 'Pathway', 'DrugClass'])
field_map = {}
drug_data = {}
for line in open(args[0], 'r'):
	
	parts = line.rstrip().split("\t")
	if not header:
		header = True
		# save the index for each data field
		for i in range(0, len(parts)):
			key = parts[i]
			if key in fields:
				field_map[i] = key
		continue

	attributes = {}
	for i in range(0, len(parts)):
		if i in field_map:
			key = field_map[i]
		else:
			continue

		val = parts[i]
		attributes[key] = val

	drug_data[attributes['DrugName']] = attributes	


index = 1
# parse the header line
print '\t'.join(['Index', 'Drug', 'Dose', 'Role', 'Pathway', 'Drug Class'])
for line in open(args[1], 'r'):
	parts = line.rstrip().split('\t')
	for i in range(1, len(parts)):
		drug, dose = parts[i].split(':')

		if drug not in drug_data:
			continue

		role = drug_data[drug]['Role']
		drug_class = drug_data[drug]['DrugClass']
		pathway = drug_data[drug]['Pathway']
		## 	PIK3/MTOR, MEK AND ERBB are the only drug classes we look at ...
		if drug_class == 'PI3K/MTOR' or drug_class == 'MEK' or drug_class == 'ERBB':
			print '\t'.join([str(index), drug, dose, role, pathway, drug_class])
			index += 1

	break	
