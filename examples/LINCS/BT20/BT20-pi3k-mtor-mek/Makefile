
perturbation_data="../DATA/perturbation.matrix"
viper_inferences="../viper-inferences/viperScores.30TF_whitelist.txt"

upstream-signaling-summary.txt:
	./scripts/getUpstreamGenes.py  < ../../DATA/LJP56_BT20_T24_p0.05.tsv > $@

analysis/concensus-networks-kegg.txt:
	./scripts/consensusNetworks.py \
	--subtypes subtypes.tab \
	-d cytoscape-kegg \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperScores.30TF_whitelist.txt \
	-n ../../pathway-priors/pi3k-mek.sif \
	-z 4 > $@

tiedie/input:
	mkdir -p $@
	cp ${perturbation_data} $@	
	cp ${viper_inferences} $@	

analysis/net_sig_match_egfr.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_egfr.txt --output analysis/net_sig_match_egfr.REPORT

analysis/net_sig_match_erbb.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_erbb.txt --output analysis/net_sig_match_erbb.REPORT
