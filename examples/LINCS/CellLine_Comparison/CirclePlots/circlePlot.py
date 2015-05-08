#!/usr/bin/env python
"""circlePlot.py: 

Usage:
  circlePlot.py [options] outputDir inputFile [inputFile_1:colorSpecFile_1 ...]

Options:
  -s str		list file containing samples to include
  -f str		list file containing features to include
  -c str		file to use as center colors
  -c str		file to use as center colors
  -l			print the feature identifier in the circle or not (default: FALSE)
  -q			run quietly
"""
## Written By: Steve Benz and Zack Sanborn
## Modified By: Sam Ng
## Last Updated: 10/10/2011
import getopt, math, os, sys, re
from matplotlib import *
use('Agg')
from pylab import *
from random import random
import mData

verbose = True
tstep = 0.01

class rgb:
	def __init__(self,r,g,b):
		self.r = int(round(r))
		self.g = int(round(g))
		self.b = int(round(b))
		
		if self.r > 255:
			self.r = 255
		elif self.r < 0:
			self.r = 0
		if self.g > 255:
			self.g = 255
		elif self.g < 0:
			self.g = 0
		if self.b > 255:
			self.b = 255
		elif self.b < 0:
			self.b = 0
		
	def tohex(self):
		r = self.r
		g = self.g
		b = self.b
		hexchars = "0123456789ABCDEF"
		return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]

def usage(code = 0):
	print __doc__
	if code != None: sys.exit(code)

def log(msg, die = False):
	if verbose:
		sys.stderr.write(msg)
	if die:
		sys.exit(1)

def syscmd(cmd):
	log("running:\n\t"+cmd+"\n")
	exitstatus = os.system(cmd)
	if exitstatus != 0:
		print "Failed with exit status %i" % exitstatus
		sys.exit(10)
	log("... done\n")

def scmp(a, b, feature, dataList):
	dataFeature = feature
	if (a not in dataList[0]) & (b in dataList[0]):
		return(1)
	elif (a in dataList[0]) & (b not in dataList[0]):
		return(-1)
	elif (b not in dataList[0]) & (a not in dataList[0]):
		return(0)
	if dataFeature not in dataList[0][a]:
		if "*" in dataList[0][a]:
			dataFeature = "*"
		else:
			return(0)
	val = cmp(dataList[0][a][dataFeature], dataList[0][b][dataFeature])
	if val == 0:
		if len(dataList) > 1:
			val = scmp(a, b, feature, dataList[1:])
		else:
			return(0)
	return(val)

def polar(r, val):
	theta = -2.0 * math.pi * val + math.pi/2.0
	x = r * math.cos(theta)
	y = r * math.sin(theta)
	return x, y

def getColor_Range(val, maxColor = rgb(0, 0, 0), minColor = rgb(0, 0, 0)):
	"""
		Gets the intermediate blended RGB color corresponding to two 'endpoint'
		colors, and a value between 0 and 1 that represents the relative weight.

		Input: 

		val: floating point value between 0 and 1
		minColor: the fade-in color to use for the start of the range
		maxColor: the fade-out color at the end of the range

		Returns:

		An 'rgb' color object
	"""
	fval = float(val)
	if fval < 0 or fval > 1:
		raise Exception("Error: value outside of range")

	zeroColor = rgb(0,0,0)

	r = fval * float(maxColor.r) + (1-fval)*minColor.r
	g = fval * float(maxColor.g) + (1-fval)*minColor.g
	b = fval * float(maxColor.b) + (1-fval)*minColor.b

	col = rgb(r, g, b)

	return col.tohex()


def mapValue_ColorRange(val, color_scheme_def):
	"""
	Get a positive value in interval [ , )

	"""

	color = None
	upper_bound = None
	lower_bound = None
	for (start, end) in color_scheme_def:
		if val >= start and val < end:
			normalized_val = (val - start)/float(end-start)
			rgbStart, rgbEnd = color_scheme_def[(start,end)]			
			return (normalized_val, rgbStart, rgbEnd)

		# largest 
		if val > end:
			upper_bound = (start, end)
		elif val < start:
			lower_bound = (start, end)
			


	# if this value is above the highest range, or below the lowest, assign
	# the max color and return
	rgbStart = None
	rgbEnd = None
	normalized_val = None
	if upper_bound:
		rgbStart, rgbEnd = color_scheme_def[upper_bound]
		normalized_val = 1.0
	elif lower_bound:
		rgbStart, rgbEnd = color_scheme_def[lower_bound]
		normalized_val = 0.0
	return (normalized_val, rgbStart, rgbEnd)

def getColor(val, color_scheme_def):

	try:
		fval = float(val)
		if fval != fval:
			raise ValueError
	except ValueError:
		col = rgb(200,200,200)
		return col.tohex()

	# zero is always white
	zeroColor = rgb(255,255,255)

	# get the color based on where the value falls in the range, according to the 
	# user-provided input spec, and normalize the value:
	# 1 is exactly this color, 0 is just white. 
	relativeVal, rgbStart, rgbEnd = mapValue_ColorRange(val, color_scheme_def)

	return getColor_Range(relativeVal, rgbEnd, rgbStart)
	## convert value to a point on the gradient
	#r = fval * float(rgbEnd.r - rgbStart.r) + rgbStart.r
	#g = fval * float(rgbEnd.g - rgbStart.g) + rgbStart.g
	#b = fval * float(rgbEnd.b - rgbStart.b) + rgbStart.b
#
#	
#	try:
#		col = rgb(r,g,b)
#	except ValueError:
#		col = rgb(200,200,200)
#	return col.tohex()

def plotScale(imgFile, minVal, maxVal):
	imgSize = (2, 4)
	fig = plt.figure(figsize=imgSize, dpi=100, frameon=True, facecolor='w')
	for i in range(10):
		val = minVal+i*(maxVal-minVal)/10
		col = getColor(val, minVal, maxVal)
		X = [float(i)/10, float(i+1)/10, float(i+1)/ 10, float(i)/10, float(i)/10]
		Y = [1, 1, 0, 0, 1]
		fill(X, Y, col, lw = 1, ec = col)
	savefig(imgFile)
	close()

def plotCircle(imgFile, label = "", centerCol = rgb(255, 255, 255).tohex(), circleCols = [[rgb(200, 200, 200).tohex()]], innerRadTotal=0.2, outerRadTotal=0.5, width = 5):
	## image settings
	imgSize = (width, width)
	fig = plt.figure(figsize=imgSize, dpi=100, frameon=True, facecolor='w')
	axes([0, 0, 1, 1], frameon=True, axisbg='w')
	axis('off')
	circleWid = (outerRadTotal-innerRadTotal)/float(len(circleCols))
	
	## color center
	outerRad = innerRadTotal
	outerRad -= .01
	X = []
	Y = []
	x, y = polar(outerRad, 0)
	X.append(x)
	Y.append(y)
	ti = 0
	while ti < 1:
		x, y = polar(outerRad, ti)
		X.append(x)
		Y.append(y)
		ti += tstep
		if ti > 1:
			break
	x, y = polar(outerRad, 1)
	X.append(x)
	Y.append(y)
	fill(X, Y, centerCol, lw = 1, ec = centerCol)
	
	## color rings
	for i in range(len(circleCols)):
		innerRad = (i*circleWid)+innerRadTotal
		outerRad = ((i+1)*circleWid)+innerRadTotal-.01
		for j in range(len(circleCols[i])):
			t0 = float(j)/len(circleCols[i])
			t1 = float(j+1)/len(circleCols[i])
			X = []
			Y = []
			x, y = polar(innerRad, t0)
			X.append(x)
			Y.append(y)
			ti = t0
			while ti < t1:
				x, y = polar(outerRad, ti)
				X.append(x)
				Y.append(y)
				ti += tstep
				if ti > t1:
					break
			x, y = polar(outerRad, t1)
			X.append(x)
			Y.append(y)
			ti = t1
			while ti > t0:
				x, y = polar(innerRad, ti)
				X.append(x)
				Y.append(y)
				ti -= tstep
				if ti < t0:
					break
			x, y = polar(innerRad, t0)
			X.append(x)
			Y.append(y)
			fill(X, Y, circleCols[i][j], lw = 1, ec = circleCols[i][j])
	
	## save image
	text(0, 0, label, ha='center', va='center')
	xlim(-0.5, 0.5)
	ylim(-0.5, 0.5)
	savefig(imgFile)
	close()

def parseColorScheme(file):

	map = {}	
	fh = open(file, 'r')
	for line in fh:
		parts = line.rstrip().split('\t')
		col_range = tuple([float(v) for v in parts[0].split(':')])
		rgbA, rgbB = parts[1].split(':')
		tmpX, tmpY, tmpZ = [int(v) for v in rgbA.rstrip(']').lstrip('[').split(',')]
		rgbA = rgb(tmpX, tmpY, tmpZ)
		tmpX, tmpY, tmpZ = [int(v) for v in rgbB.rstrip(']').lstrip('[').split(',')]
		rgbB = rgb(tmpX, tmpY, tmpZ)
		map[col_range] = (rgbA, rgbB)

	return map

def main(args):
	## parse arguments
	try:
		opts, args = getopt.getopt(args, "s:f:o:c:lq")
	except getopt.GetoptError, err:
		print str(err)
		usage(2)
	if len(args) < 2:
		usage(2)
	
	outputDir = args[0].rstrip("/")
	circleFiles = args[1:]
	
	sampleFile = None
	featureFile = None
	centerFile = None
	printLabel = False
	global verbose
	for o, a in opts:
		if o == "-s":
			sampleFile = a
		elif o == "-f":
			featureFile = a
		elif o == "-o":
			sa = re.split(";", a)
			if len(sa) == 1:
				orderFeature = sa[0]
				orderFiles = []
			else:
				orderFeature = sa[0]
				orderFiles = re.split(",", sa[1])
		elif o == "-c":
			centerFile = a
		elif o == "-l":
			printLabel = True
		elif o == "-q":
			verbose = False
	
	## execute
	samples = []
	features = []
	if sampleFile != None:
		samples = mData.rList(sampleFile)
	if featureFile != None:
		features = mData.rList(featureFile)

	print features
	## read circleFiles
	circleData = []
	circleColors = []
	##
	## record file types for each, effects the color scheme 
	## use the input index for each
	color_scheme_map = {}
	for i in range(len(circleFiles)):
		circleFile, colorScheme = circleFiles[i].split(':')
		color_scheme_map[i] = parseColorScheme(colorScheme)
		(data, cols, rows) = mData.rCRSData(circleFile, retFeatures = True)
		circleData.append(data)
		#circleColors.append( (minCol, zerCol, maxCol) )
		#if sampleFile == None:
		#	samples = list(set(cols) | set(samples))
		#if featureFile == None:
		#	features = list(set(rows) | set(features))
	
	## read centerFile
	centerData = None
	if centerFile != None:
		centerData = mData.r2Col(centerFile, header = True)
		
	## plot images
	for feature in features:
		log("Drawing %s\n" % (feature))
		imgName = re.sub("[/:]", "_", feature)
		if len(imgName) > 100:
			imgName = imgName[:100]
		imgFile = "%s/%s.png" % (outputDir, imgName)
		label = ""
		if printLabel:
			label = feature
		centerCol = rgb(255, 255, 255).tohex()
		if centerData != None:
			if feature in centerData:
				minVal = min([-0.01]+mData.floatList(centerData.values()))
				maxVal = max([0.01]+mData.floatList(centerData.values()))
				centerCol = getColor(centerData[feature], minVal, maxVal)
				log("\t%s,%s,%s,%s\n" % (centerData[feature],minVal,maxVal,centerCol))
		circleCols = []
		for i in range(len(circleData)):
			ringCols = []
			ringVals = []
			for sample in samples:
				if sample in circleData[i]:
					if feature in circleData[i][sample]:
						ringVals.append(circleData[i][sample][feature])
					elif "*" in circleData[i][sample]:
						ringVals.append(circleData[i][sample]["*"])
			minVal = min([-0.01]+mData.floatList(ringVals))
			maxVal = max([0.01]+mData.floatList(ringVals))
			for sample in samples:
				if sample in circleData[i]:
					if feature in circleData[i][sample]:
						ringCols.append(getColor(circleData[i][sample][feature], color_scheme_map[i]))
					elif "*" in circleData[i][sample]:
						ringCols.append(getColor(circleData[i][sample]["*"], color_scheme_map[i]))
					else:
						ringCols.append(rgb(200, 200, 200).tohex())
				else:
					ringCols.append(rgb(200, 200, 200).tohex())
			circleCols.append(ringCols)
		plotCircle(imgFile, label = label, centerCol = centerCol, circleCols = circleCols, innerRadTotal=0.2, outerRadTotal=0.5, width = 5)

if __name__ == "__main__":
	main(sys.argv[1:])
