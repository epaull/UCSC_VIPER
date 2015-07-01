#!/usr/bin/env  python

import os, sys

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()

all_edges = {}

if len(args) > 0:
        files = args
else:
        files = ['flat_directed_pathways/'+file for file in os.listdir('flat_directed_pathways')]


for file in files:
        fh = open(file, 'r')
        for line in fh:
                src, dest, direction, type = line.rstrip().split('\t')
                if src == 'src':
                        continue

                # just merge and set the edge to 'conflicted'
                if (direction == 'directed' and (dest, src) in all_edges):
                        #print "Direction Conflict: "+src+'\t'+dest
                        #sys.exit(1)
                        direction = 'conflicted'

                all_edges[(src, dest)] = (direction, type)

        fh.close()

# final pass for order indepedence: 
# update
filtered_edges = {}
for edge in all_edges:
        src = edge[0]
        dest = edge[1]
        if (dest, src) in all_edges:
                # update the other direction to 'conflicted' if this is the case
                filtered_edges[(dest, src)] = ('conflicted', all_edges[(dest, src)][1])
        else:
                filtered_edges[(src, dest)] = all_edges[(src, dest)]

for edge in filtered_edges:
        print '\t'.join(edge)+'\t'+'\t'.join(filtered_edges[edge])

