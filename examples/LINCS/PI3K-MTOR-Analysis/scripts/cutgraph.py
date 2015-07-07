#!/usr/bin/env	python2.7

from tiedie_util import *
from linkers import *
import networkx as nx
import math

class cutGraph:
	"""
	Take a list of scores, rank them, and return the subnetworks
	found at that threshold.
	"""

	# how many nodes to step over when searching
	STEP_SIZE=50
	# maximum search range
	MAX_SEARCH=1000

	def __init__(self, network, directed=True, verbose=False, max_search_size=1000, step_size=50, penalty_multiplier=(1,5)):
		"""
			Converts the hash format network into a networkX graph. 
			
			scores : indexed by gene
			base_network : 
		"""

		# convert to a networkX object. 
		self.G = cutGraph.getNXgraph(network, directed)

		self.VERBOSE = verbose
		self.MAX_SEARCH = max_search_size
		self.STEP_SIZE = step_size
		self.MIN_PENALTY_MULTIPLIER = float(penalty_multiplier[0])
		self.MAX_PENALTY_MULTIPLIER = float(penalty_multiplier[1])

	def setScores(self, gene_scores, input_set):
		"""
		Rank the set of scores here

		input_set:
			the set of all genes, either upstream or downstream input nodes
			to be considered when constructing networks (i.e. the second
			category of nodes--linkers are the other). 
		"""
		self.input_set = input_set
		# convert scores to absolute value
		node_universe = self.G.nodes()
		for gene in gene_scores:
			# remove any nodes not in the network from the list
			if gene not in node_universe:
				continue
			gene_scores[gene] = abs(gene_scores[gene])

		self.ranked_genes = []
		# store the rank of 
		self.gene_ranks = {}
		i = 1
		for (gene, val) in sorted(gene_scores.items(), key=operator.itemgetter(1), reverse=True):
			self.ranked_genes.append(gene)
			self.gene_ranks[gene] = i
			i += 1

	def supervisedSearch(self, size_factor):
		"""
		Use TieDIE 1.0 logic to get the network based on a user-supplied size threshold
		"""
		# remove the input genes, rank the rest
		threshold_idx = int(math.floor(len(self.input_set)*size_factor))
		return self.getSubgraphs(self.G, threshold_idx)

	def kSubnetsSearch(self, k_min=3, k_max=15):
		self.k_min = k_min
		self.k_max = k_max
		threshold, subgraph = self.kSubnetsSearch_Internal(self.G, self.STEP_SIZE, self.MAX_SEARCH+1, self.STEP_SIZE)
		return threshold, subgraph

	def cutGraph_SW(self, g, k=25):
		"""
		Define edge weights based on the delta between node weights. 
		Use stoer_wagner to cut the graph into subgraphs of max node size (k)
		"""

		cut_subgraph = nx.Graph()
		for edge in g.edges():
		
			rank_A = self.gene_ranks[edge[0]]
			rank_B = self.gene_ranks[edge[1]]
			# set the weight as the relative delta in rankings 
			w = 1.0 - abs(rank_A-rank_B)/float(len(self.ranked_genes))
			g[edge[0]][edge[1]]['rank_delta'] = w

		cut_value, partitions = nx.stoer_wagner(g, weight='rank_delta')
		for partition in partitions:
			gc = None
			if len(partition) > k:
				# recurse and cut again until it meets the threshold
				gc = self.cutGraph_SW(g.subgraph(partition), k)	
			else:
				gc = g.subgraph(partition)
			# get the edges connecting this partition and add them to the new graph
			cut_subgraph.add_edges_from(gc.edges())
			# now add edge data
			for edge in gc.edges():
				cut_subgraph[edge[0]][edge[1]]['i'] = g[edge[0]][edge[1]]['i']

		return cut_subgraph

	def scoreExpressionRange(self, start, stop, step_size, source_set, downstream_expression):
		"""
		Search over cutoffs to find the threshold that gives the maximal 
		amount of explanatory power for an expression dataset, relative to
		network size
		"""

		# maintain a score for each threshold based on number of graphs of size K or above
		scores = {}

		# for the range function
		if stop - start < 0:
			tmp = start
			start = stop
			stop = tmp

		# linear search
		fractions = {}
		for threshold in range(start, stop, step_size):
			# FIMME: use has_path, threshold on the length if necessary. 
			# use node_connected_component(G,n) for each additional 
			# node and look at just those sources, targets
			covered_genes, coverage_fraction, covered_tfs, num_linkers = self.scoreThreshold(source_set, threshold, downstream_expression)
			fractions[threshold] = (coverage_fraction, covered_tfs, num_linkers)

		return fractions	

	def searchOptScore(self, start, stop, step_size, paths_to_targets, downstream_expression):
		"""
		Search over cutoffs to find the threshold that gives the maximal 
		amount of explanatory power for an expression dataset, relative to
		network size
		"""

		# maintain a score for each threshold based on number of graphs of size K or above
		scores = {}

		# for the range function
		if stop - start < 0:
			tmp = start
			start = stop
			stop = tmp

		# linear search
		scores = {}
		covered_genes = {}
		for threshold in range(start, stop, step_size):
			covered_genes, coverage_fraction, score = self.scoreThreshold(paths_to_targets, threshold, downstream_expression)
			scores[threshold] = score
			covered_genes[threshold] = covered_targets

		sorted_tuples = sorted(scores.items(), key=operator.itemgetter(1), reverse=True)

		threshold = None
		subgraph = None
		covered_targets = None

		if step_size == 1:
			threshold = sorted_tuples[0][0]
			subgraph = self.getSubgraphs(self.G, threshold)
			covered_targets = covered_genes[threshold]
		else:		
			# otherwise recurse once and search over a smaller range. Use the entire graph again
			threshold, subgraph, covered_targets = self.searchOptScore(sorted_tuples[0][0], sorted_tuples[1][0], int(math.ceil(step_size/10.0)), paths_to_targets, downstream_expression)
		
		return (threshold, subgraph, covered_targets)	

	def kSubnetsSearch_Internal(self, g, start, stop, step_size, recurse=True):
		"""
		Search over cutoffs to find the threshold that gives the maximal number of 
		subnetworks of at least size K.
		"""

		# maintain a score for each threshold based on number of graphs of size K or above
		scores = {}

		# for the range function
		if stop - start < 0:
			tmp = start
			start = stop
			stop = tmp

		for threshold in range(start, stop, step_size):
			subgraph = self.getSubgraphs(g, threshold)
			num_gt_k = 0
			for cc in nx.connected_components(subgraph):
				if len(cc) >= self.k_min:
					num_gt_k += 1
			scores[threshold] = num_gt_k
			if self.VERBOSE:
				print "Threshold - Num CC: "+str(threshold)+' '+str(num_gt_k)
		

		sorted_tuples = sorted(scores.items(), key=operator.itemgetter(1), reverse=True)
		# now do a binary search between the top 2 ranges
		threshold = None
		ccm = None
		if step_size == 1:
			# return the best solution
			threshold = sorted_tuples[0][0]
			subgraph_raw = self.getSubgraphs(g, threshold)
			subgraph = nx.Graph()
			for cc in nx.connected_components(subgraph_raw):
				cc_subgraph = g.subgraph(cc)
				cut_graph = None
				if self.k_max and len(cc) > self.k_max:
					# cut with the flow algorithm if necessary
					cut_graph = self.cutGraph_SW(cc_subgraph, self.k_max)
				else:
					cut_graph = cc_subgraph
	
				# now add the edges of this modified or unmodified connected component	
				subgraph.add_edges_from(cut_graph.edges())
				for edge in cut_graph.edges():
					subgraph[edge[0]][edge[1]]['i'] = g[edge[0]][edge[1]]['i']
					

			# foreach subgraph exceeding the max size self.k_max, do another stepwise search, get the 
			# broken up components of this subgraph
			# create a new graph object, build it up with new edges
			#minced_graph = nx.Graph()
			#for cc in nx.connected_components(subgraph):
			#	cc_subgraph = subgraph.subgraph(cc)
			#	if len(cc) > self.k_max:
			#		print len(cc)
			#		# recurse and get smaller subgraphs
			#		t, ccM = self.kSubnetsSearch_Internal(cc_subgraph, 1, len(cc), 1, False)
			#		minced_graph.add_edges_from(ccM.edges())
			#		for edge in ccM.edges():
			#			minced_graph[edge[0]][edge[1]]['i'] = ccM[edge[0]][edge[1]]['i']
			#	else:
			#		minced_graph.add_edges_from(cc_subgraph.edges())
			#		for edge in cc_subgraph.edges():
			#			minced_graph[edge[0]][edge[1]]['i'] = cc_subgraph[edge[0]][edge[1]]['i']
			#
			#subgraph = minced_graph
		else:		
			# otherwise recurse once and search over a smaller range. Use the entire graph again
			threshold, subgraph = self.kSubnetsSearch_Internal(g, sorted_tuples[0][0], sorted_tuples[1][0], int(math.ceil(step_size/10.0)), False)

		# return a 
		return (threshold, subgraph)

	@staticmethod
	def isDirectedEdge(i_str):

		directed = True
		if i_str == 'PPI>':
			directed = False
		elif i_str == 'INTERACTS>':	
			directed = False

		return directed

	@staticmethod
	def getNXgraph(network, directed=True):
	
		"""
			Convert a hash-style network object into a directed nx graph object
		"""	
		# use networkx to find the largest connected sub graph
		G = None
		if directed:
			G = nx.DiGraph()

		else:
			G = nx.Graph()
	
		for s in network:
			for (i,t) in network[s]:
				# ignore self-links
				if s == t:
					continue
				G.add_edge(s,t)
				G[s][t]['i'] = i

				# with a directed graph, add a duplicate edge in the other 
				# direction to create an undirected edge if the edge type matches
				if directed and not cutGraph.isDirectedEdge(i):
					G.add_edge(t,s)
					G[t][s]['i'] = i
	
		return G


	def getSubgraphs(self, g, threshold_idx):
		"""
		Returns: a networkX graph with edges between just this set of nodes and the 
		top X genes
		"""	
		gene_set = set()
		for gene in self.ranked_genes[0:threshold_idx]:
			gene_set.add(gene)
		if self.input_set:
			for gene in self.input_set:
				gene_set.add(gene)

		subgraph = g.subgraph(gene_set)		
		
		return subgraph

	def scoreThreshold(self, source_set, linker_cutoff_rank, downstream_expression):
		"""
		Use a heuristic: assuming no overlap the minimal possible number of linker genes
		to connect the set is X = (k - 1)/2 * (size of input set), where k is the minimum number of edges. 
		The other part of the heuristic is the maximum penalty to ramp up; call it 5 by default:
		X*5
		"""

		# includes the union of input set genes and these linkers...	
		g = self.getSubgraphs(self.G, linker_cutoff_rank)

		expression_gene_universe = set()
		# all genes covered by an active TF in this set
		tfs_covered = set()
		expression_covered = set()
		for target in downstream_expression.keys():
			if target not in g.nodes():
				continue
			# if you find a single path from any source, add it and continue
			for source in source_set:	
				if source not in g.nodes():
					continue
				if nx.has_path(g, source, target):
					expression_covered = expression_covered.union(downstream_expression[target])
					tfs_covered = tfs_covered.union(target)
					break

		for target in downstream_expression:
			expression_gene_universe = expression_gene_universe.union(downstream_expression[target])

		# the maximization term
		if len(expression_gene_universe) == 0:
			return None
		coverage_fraction = float(len(expression_covered))/len(expression_gene_universe)
		tfs_covered_fraction = len(tfs_covered)/float(len(downstream_expression.keys()))

		# count the number of linker nodes
		num_linker_nodes = len(set(g.nodes()).difference(self.input_set))

		# return the score along with the genes actually covered 
		return (expression_covered, coverage_fraction, tfs_covered_fraction, num_linker_nodes)
		#return (expression_covered, coverage_fraction, float(coverage_fraction - min_term))

	def computeConnectedTargets(self, paths_to_targets, linker_cutoff_rank):
		"""
		For each target, evalute all it's paths to see if a linker is below
		the cutoff. One 'broken' linker in each path will invalidate the path.
		If all paths are invalidated, the target it disconnected. 

		Returns:
			A set of targets that are still connected. 
		"""
		connected_targets = set()
		for target in paths_to_targets:

			# assume disconnected until we find one path that still connects
			# at this threshold
			disconnected = True
			for path in paths_to_targets[target]:

				connected_path = True
				for i in range(1, len(path)-1):
					linker_gene = path[i]
					if self.gene_ranks[linker_gene] > linker_cutoff_rank:
						connected_path = False
						break
				if connected_path:
					disconnected = False
					break

			if not disconnected:
				connected_targets.add(target)

		return connected_targets	
						

	def getLinkingPaths(self, source_set, target_set):
		"""
		Build a hash of all source -> target paths, indexed by target
		"""
		# build an index of paths for each target
		target_paths = defaultdict(set)
		for source in source_set:
			if source not in self.G.nodes():
				continue

			for target in target_set:
				if target not in self.G.nodes():
					continue
				for path in nx.all_simple_paths(self.G, source, target, cutoff=int(self.search_depth)):
					target_paths[target].add(tuple(path))
	
		return target_paths	

	def tiediePSNsearch(self, k, source_set, downstream_expression):
		"""
		Search all paths from each gene in the source set to any 
		gene in the target set up to length k. 

		Search linker thresholds, score each based on the number of linker
		genes included and the number of genes explained by each network

		target_scores:
			hash, target genes, each points to a (possibly overlapping)
			set of targets that are dis-regulated in this patient

		downstream_expression:
			a hash of active TFs in this patient, each pointing to a set of 
			genes that each regulates	
			
		"""
		# Use a heuristic: assuming no overlap the minimal possible number of linker genes
		# to connect the set is X = (k - 1)/2 * (size of input set), where k is the 
		# minimum number of edges. 
		# The other part of the heuristic is the maximum penalty to ramp up; call it 5 by default:
		self.search_depth = k
		self.c_linker_min = float(len(self.input_set))*((self.search_depth-1)/2.0)
		self.c_linker_max = int(math.ceil(self.c_linker_min*self.MAX_PENALTY_MULTIPLIER))
		# search all paths, source to target
		for edge in self.G.edges():
			rank_A = self.gene_ranks[edge[0]]
			rank_B = self.gene_ranks[edge[1]]
			# set the weight as the relative delta in rankings 
			w = abs(rank_A-rank_B)/float(len(self.ranked_genes))
			self.G[edge[0]][edge[1]]['rank_delta'] = w

		# build an index of paths for each target
		#target_paths = self.getLinkingPaths(source_set, set(downstream_expression.keys()))
		#if len(target_paths) == 0:
		#	return None

		scores = self.scoreExpressionRange(1, self.c_linker_max, 1, source_set, downstream_expression)
		return scores

class tiedieNetwork:
	"""
	Encapsulates a tiedie heat vector, gene inputs and a networkX graph subnetwork
	related to each
	"""

	def __init__(self, nxGraph, input_events, tiedie_vector):


		self.g = nxGraph
		# indexed by gene -- each contains a set of event types
		self.input_events = input_events
		# indexed by gene
		self.tiedie_vector = tiedie_vector

	def getLinkerGenes(self, threshold=1.5):
		"""
		Return linkers that are not in the input_event set
		"""

	@staticmethod	
	def parseLinkers(name):
		# indexed by sample, then by gene
		sample_vectors = {}
		header = False
		for line in open(name, 'r'):
			parts = line.rstrip().split("\t")
			if not header:
				header = parts
				continue
	
			sample = parts[0]
			sample_vectors[sample] = {}
			for i in range(1, len(header)):
				sample_vectors[sample][header[i]] = float(parts[i])	
	
		return sample_vectors

	@staticmethod	
	def parseEventsFile(name):
		# indexed by sample, then by gene: points to a set containing all event types that occur
		sample_vectors = {}
		header = False
		for line in open(name, 'r'):
			parts = line.rstrip().split("\t")
			if not header:
				header = parts
				continue
	
			gene = parts[0]
			type = parts[1]
			for i in range(2, len(header)):
				sample = header[i]
				if sample not in sample_vectors:
					sample_vectors[sample] = {}
				if gene not in sample_vectors[sample]:
					sample_vectors[sample][gene] = set()
				sample_vectors[sample][gene].add(type)
				
		return sample_vectors
	
	
	@staticmethod
	def fromFile(network, k, linkerFile, eventsFile):
		"""
		network: a hash-valued network object
		k: the minimum size subnet to optimize
		linkerFile: linker heats output by the initial tiedie
		run
		eventsFile: the events/inputs file printed out by tiedie.PSN
		"""
		tiedie_vectors = tiedieNetwork.parseLinkersFile(linkerFile)
		input_events = tiedieNetwork.parseEventsFile(eventsFile)

		network = parseNet(NETWORK)
		graph = cutGraph(network, False)

		# return a bunch of network objects indexed by sample id
		networkObjects = {}

		for sample in sorted(tiedie_vectors.keys()):
			graph.setScores(tiedie_vectors[sample], None)
			# look at the top 100 genes
			#subgraphs = graph.getSubgraphs(100)
			#print subgraphs.edges()
			threshold, g = graph.kSubnetsSearch(k)
			netObj = tiedieNetwork(g, input_events[sample], tiedie_vectors[sample])
			networkObjects[sample] = netObj
			
		return networkObjects

if __name__ == '__main__':
	unittest.main()

			
	
	
	
	
