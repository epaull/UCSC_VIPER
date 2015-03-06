

upstream-signaling-summary.txt:
	./scripts/getUpstreamGenes.py  < LJP56_BT20_T24_p0.05.tsv > $@

viper-output:
	mkdir -p $@
	../../../bin/run-viper.R -o viper-output -n ../../../pathways/brca-tf-regulon.rda -e data.entrez.tab --unsupervised --phenotypes phenotypes.tab --regulon_minsize 5
