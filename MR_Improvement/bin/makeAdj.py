#!/usr/bin/env python2.7

# find the downstream TF set of each
# matrix: i.e. "gene event <samples>"

import sys, os
from optparse import OptionParser
from collections import defaultdict
parser = OptionParser()
(opts, args) = parser.parse_args()

def parseEDGES(file):

    regs = {}
    for line in open(file, 'r'):
        regulator, target, MI = line.rstrip().split("\t")
        if regulator.startswith("Transcriptional"):
            continue
        if regulator not in regs:
            regs[regulator] = set()
        regs[regulator].add( (target, MI) )
    return regs

regulators = parseEDGES(args[0])

for regulator in regulators:

    printstr = regulator
    for (target, MI) in regulators[regulator]:
        printstr += '\t'+target+'\t'+MI
    print printstr
