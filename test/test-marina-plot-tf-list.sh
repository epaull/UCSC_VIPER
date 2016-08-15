
mkdir test-marina-plot-tf-list
../bin/run-marina.R -e data/test.data.tab -p data/phenotypes.tab -n data/multinet.adj -t Tumor -r Normal --permutations 100 --output test-marina-plot-tf-list --tfs_to_plot data/TFs-to-plot.lst

