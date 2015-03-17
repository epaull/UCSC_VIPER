

upstream-signaling-summary.txt:
	./scripts/getUpstreamGenes.py  < LJP56_BT20_T24_p0.05.tsv > $@

viper-output:
	mkdir -p $@
	../../../bin/run-viper.R -o viper-output -n ../../../pathways/brca-tf-regulon.rda -e data.entrez.tab --unsupervised --phenotypes phenotypes.tab --regulon_minsize 5


consensus-network-diff.txt:
	./scripts/consensusNetworks.py \
	-x subtypes.tab \
	-d cytoscape \
	-s input/tiedie/upstream_mek_mtor_pi3k.tab \
	-t input/tiedie/viperScores.txt \
	-n ../pathway-priors/superpathway-2.1.sif 
