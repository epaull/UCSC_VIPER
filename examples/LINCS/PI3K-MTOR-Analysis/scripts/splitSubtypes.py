#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()

from tiedie_util import *

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

samples_A, samples_B = parseSubtypes(args[0], set(['AKT_PI3K', 'MTOR']), set(['MAPK']))

activities = parseMatrix(args[1], None, 0)

samples_A = list(samples_A.intersection(set(activities.keys())))
samples_B = list(samples_B.intersection(set(activities.keys())))
regulators = None
for sample in activities:
	regulators = activities[sample].keys()	
	

fh = open('samples_AktPI3K.txt', 'w')
fh.write('Key'+'\t'+'\t'.join(samples_A)+'\n')
for r in regulators:
	printstr = r
	for sample in samples_A:
		printstr += '\t'+str(activities[sample][r])
	fh.write(printstr+'\n')
fh.close()

fh = open('samples_MAPK.txt', 'w')
fh.write('Key'+'\t'+'\t'.join(samples_B)+'\n')
for r in regulators:
	printstr = r
	for sample in samples_B:
		printstr += '\t'+str(activities[sample][r])
	fh.write(printstr+'\n')
fh.close()

