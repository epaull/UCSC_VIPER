
perturbation_data="../DATA/perturbation.matrix"
viper_inferences="../viper-inferences/viperScores.30TF_whitelist.txt"

tiedie/input:
	mkdir -p $@
	cp ${perturbation_data} $@	
	cp ${viper_inferences} $@	
	transpose.pl ${perturbation_data} > 1.tmp
	transpose.pl ${viper_inferences} > 2.tmp

analysis/consensus-networks-kegg-small.txt:
	./scripts/consensusNetworks.py \
	--subtypes drug_classes.tab \
	-d tiedie/cytoscape-kegg-small \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperScores.30TF_whitelist.txt \
	-n ../../pathway-priors/pi3k-mek.sif \
	-z 4 > $@

analysis/consensus-networks-kegg-large.txt:
	./scripts/consensusNetworks.py \
	--subtypes drug_classes.tab \
	-d tiedie/cytoscape-kegg-large \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperScores.30TF_whitelist.txt \
	-n ../../pathway-priors/pi3k-mek.sif \
	-z 4 > $@

analysis/net_sig_match_test.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_test.txt --output $@ 
