# This Makefile computes the diffusion kernel of a network.
# The graph file should be a list of tab-delimited pairs, one
# pair on each line of the file. The main target is an n-by-n
# kernel matrix, where n is the number of unique nodes in the
# graph.

in          = ../signal-transduction.sif
matlab_code = $(MAPDIR)/Templates/Matlab/network_diffusion_kernel.m

targets  = \
	   data.tab.gz \

include $(MAPDIR)/Templates/Make/quick.mak

startup.m: ${in}
	cut -f 1,3 $< \
	| order_keys.pl \
	| sort -k 1,2 -t '	' -u \
	> in.tab; \
	echo 'exit;' \
	| cat $(matlab_code) - \
	> startup.m; 

kernel.mat: data.tab.gz

# make kernel_TIME.tab
kernel.tab: in.tab startup.m
	matlab -nodesktop -nojvm -nosplash; \
