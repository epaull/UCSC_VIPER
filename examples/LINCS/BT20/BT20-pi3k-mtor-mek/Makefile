

upstream-signaling-summary.txt:
	./scripts/getUpstreamGenes.py  < LJP56_BT20_T24_p0.05.tsv > $@


viper-inferences-unsupervised:
	mkdir -p $@
	../../bin/run-viper-unsupervised.R -o $@ -n ../../pathways/brca-tf-regulon.rda -e data/expression.entrez.tab


concensus-networks-kegg.txt:
	./scripts/consensusNetworks.py \
	--subtypes subtypes.tab \
	-d cytoscape-kegg \
	-s input/tiedie/upstream_mek_mtor_pi3k.tab \
	-t input/tiedie/viperScores.txt \
	-n ../pathway-priors/pi3k-mek.sif \
	-z 4 > $@


analysis/net_sig_match_egfr.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_egfr.txt --output analysis/net_sig_match_egfr.REPORT

analysis/net_sig_match_erbb.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_erbb.txt --output analysis/net_sig_match_erbb.REPORT
