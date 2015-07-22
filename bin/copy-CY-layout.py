#!/usr/bin/env	python

import sys, os, subprocess
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--from_key",dest="from_key",action="store",type="string",default=None)
parser.add_option("--to_key",dest="to_key",action="store",type="string",default=None)
(opts, args) = parser.parse_args()

copy_layout_from = args[0]
copy_layout_to = args[1]

if not copy_layout_from.endswith('.cys') or not copy_layout_to.endswith('.cys'):
	print 'Usage: <layout copy from>.cys <copy to>.cys --from_key <network view (FROM)> --to_key <network view (TO)>'
	sys.exit(1)

def getCoordsXGMLL(file):

	node_coordinates = {}
	# 0 = search for <node id, 1 = get coordinates from line, 2 = line after graphics initial state
	state = 0
	# does is end here '/>' or are there attributes?
	current_label = None
	for line in open(file, 'r'):
		line = line.rstrip().lstrip()
    	# <graphics y="359.7348954718036" x="1154.7450509691116" z="0.0">
		if state == 1:
		
			x = None	
			y = None	
			z = None	
			for axis in line.split(' ')[1:4]:
				coord, val = axis.split('=')
				val = val.lstrip('"')
				val = val.rstrip('>')
				val = val.rstrip('/')
				val = val.rstrip('"')
				if coord == "x":
					x = val
				elif coord == "y":
					y = val
				elif coord == "z":
					z = val

			node_coordinates[current_label] = (x,y,z)	
			state = 0
			current_label = None
			continue

  		# <node id="1376" label="CDK2" cy:nodeId="386">
		if line.startswith('<node id='):
			for attr in line.split(' '):
				if attr.startswith('label'):
					label = attr.split("=")[1].lstrip('"').rstrip('"')
			state = 1
			current_label = label
			continue

	return node_coordinates

def editCoords(xgmll_file, new_coordinates, output_file):

	fh = open(output_file, 'w')
	current_label = None
	# 0 = search for <node id
	# 1 = edit coordinates with current label key in new_coordinates
	state = 0
	for line in open(xgmll_file, 'r'):

		if state == 1:
			# edit the line, reset state and continue
			escape_seq = ''
			if line.rstrip().rstrip('>').rstrip('/') != line.rstrip().rstrip('>'):
				escape_seq = '/'
			
			fh.write('    <graphics x="'+new_coordinates[current_label][0]+'" y="'+new_coordinates[current_label][1]+'"'+' z="'+new_coordinates[current_label][2]+'"'+escape_seq+'>'+'\n')
			state = 0
			current_label = None
			continue
		elif line.lstrip().startswith('<node id'):
			node_label = None
			for attr in line.split(' '):
				if attr.startswith('label'):
					node_label = attr.split("=")[1].lstrip('"').rstrip('"')
			if node_label in new_coordinates:
				current_label = node_label
				state = 1
		fh.write(line.rstrip()+'\n')

	fh.close()
			

# extract source with layout
if not os.path.exists('tmp'):
	os.mkdir('tmp')

os.popen('unzip -d tmp/from '+copy_layout_from)
# extract target
os.popen('unzip -d tmp/to '+copy_layout_to)

from_xgmll = os.popen('ls tmp/from/CytoscapeSession*/views/*'+opts.from_key+'*').read().rstrip()
to_xgmll = os.popen('ls tmp/to/CytoscapeSession*/views/*'+opts.to_key+'*').read().rstrip()

if not os.path.isfile(to_xgmll):
	raise Exception("Error: cannot find destination view!")
if not os.path.isfile(from_xgmll):
	raise Exception("Error: cannot find source view!")

from_coords = getCoordsXGMLL(from_xgmll)

editCoords(to_xgmll, from_coords, 'tmp/edited_layout.xgmll')
# copy it to the new directory
os.popen('cp -R tmp/to tmp/new_layout ')
new_xgmll = 'tmp/new_layout'
for dir in to_xgmll.split('/')[2:]:
	new_xgmll += '/'+dir	
os.popen('cp tmp/edited_layout.xgmll '+new_xgmll)
os.popen('cd tmp/new_layout && zip -r session.new_CY_layout.cys CytoscapeSession* && mv session.new_CY_layout.cys ../../')
print "New session file created! session.new_CY_layout.cys"
os.popen('rm -rf tmp')

