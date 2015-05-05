#!/usr/bin/env	python

import sys, os

img_dir = os.getcwd()+'/img/'

output_file = open('imageTableCytoscape.txt', 'w')
output_file.write('name\tImageLoc\n')
for img in os.listdir(img_dir):
	img_loc = img_dir+img
	gene_name = img.rstrip('.png')
	output_file.write(gene_name+'\t'+'file://'+img_loc+'\n')
output_file.close()
