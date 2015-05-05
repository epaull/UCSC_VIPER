
perturbation_data="../DATA/perturbation.matrix"
viper_inferences="../viper-inferences/viperScores.txt"
pathway_genes="../../pathway-priors/signal-transduction.lst"
pathway_net="../../pathway-priors/signal-transduction.sif"


tiedie/input:
	mkdir -p $@
	cp ${perturbation_data} $@	
	join.pl ${pathway_genes} ${viper_inferences} > $@/viperInferences.txt

analysis/consensus-networks-kegg.txt:
	./scripts/consensusNetworks.py \
	--subtypes drug_classes.tab \
	-d tiedie/cytoscape-kegg-final \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperInferences.txt \
	-n ../../pathway-priors/signal-transduction.sif \
	-c node_equivalence_classes/defs.txt \
	-z 6 > $@

analysis/net_sig_match_test.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_test.txt --output $@ 

