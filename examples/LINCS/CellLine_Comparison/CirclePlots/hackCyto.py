#!/usr/bin/python

import sys

num = 259
for line in sys.stdin:
	gene = line.rstrip()
	printstr = "nodeCustomGraphics1.default-Node\ Custom\ Graphics\ 1-Discrete\ Mapper.mapping.map."+gene+\
	"=cytoscape.visual.customgraphic.impl.bitmap.URLImageCustomGraphics,"+str(num)+\
	",/Users/evanpaull/Desktop/KIRC/HotLink/CytoScape_SuperGene/img/"+gene+\
	".png,bitmap image"

	num += 1

	print printstr
