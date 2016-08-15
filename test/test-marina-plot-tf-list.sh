
# Run MARINa and  plot a specific list of MR TFs instead of the default top "N" TFs. 
# In addition to normal input, must pass in a list of desired TFs to be plotted, 
# one per line, no header (e.g. data/TFs-to-plot.lst below).
# Any TFs (or other characters--e.g. "TEST" here) in this list but not in the MARINa MR 
# output/results (usually masterRegulators.txt) will be silently skipped during plotting.
mkdir test-marina-plot-tf-list
../bin/run-marina.R -e data/test.data.tab -p data/phenotypes.tab -n data/multinet.adj -t Tumor -r Normal --permutations 100 --output test-marina-plot-tf-list --tfs_to_plot data/TFs-to-plot.lst

