#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()

from tiedie_util import *
from distributions import EmpiricalDist

genomic_events = parseMatrix(args[0])
# indexed by sample, then by gene
tf_activity = parseMatrix(args[1])

sample_names = list(set(genomic_events.keys()).intersection(set(tf_activity.keys())))
print '\t'.join(['id', 'Type'])+'\t'+'\t'.join(sample_names)

genes = None
for sample in genomic_events:
	genes = genomic_events[sample].keys()
	break

for gene in genes:
	printstr = gene+'\tGenomic_Event'
	for sample in sample_names:
		printstr += '\t'+str(int(genomic_events[sample][gene]))
	print printstr

# add TF activity
##
## VIPER inferences: 
##
# z-normalize by gene
viper_z_scores = EmpiricalDist(tf_activity).getZscores()

genes = None
for sample in tf_activity:
	genes = tf_activity[sample].keys()
	break

for gene in genes:
	printstr = gene+'\tTF_Activity'
	for sample in sample_names:

		score = 0.0
		# retain only Z > 1
		z_score = viper_z_scores[sample][gene]
		if z_score < -1.0:
			score = z_score	
		printstr += '\t'+str(score)
	print printstr
