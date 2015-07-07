

# Cell Line Specific Data Here
# EDIT THIS
DATA=input/BT20_T24_INFchdir.tsv


# static
targets=../shared_data/LJP56_drugs_targets.tsv
concentrations=../shared_data/drugs.conc.lst
pathway_maps=../shared_data/pathway_def_maps
node_eq=../shared_data/node_equivalence_classes/defs.txt
hugo2entrez=../shared_data/hugo2entrez.txt
entrez2hugo=../shared_data/entrez2hugo.txt
drug_classes=../shared_data/drug_classes.tab
brca_regul=../../../../pathways/brca-tf-regulon.rda

# static
pathway_genes=../../pathway-priors/signal-transduction.lst
pathway_net=../../pathway-priors/signal-transduction.sif

all: phenotypes.tab 

clean:
	rm -f phenotypes.tab
	rm -f perturbation.matrix
	rm -rf viper-inferences
	rm -rf analysis
	rm -rf tiedie/input

phenotypes.tab:
	../scripts/makePhenotypes.py ${targets} ${DATA} > $@

preprocessed/expression.entrez.tab:
	join.pl ${hugo2entrez} ${DATA} | cut -f 2- | grep -v '^	' | ../scripts/averageDupeProbes.py  > $@

viper-inferences: preprocessed/expression.entrez.tab
	mkdir -p $@
	../../../../bin/run-viper-unsupervised.R -o $@ -n ${brca_regul}  -e $<

viper-inferences/viperInferences.txt: viper-inferences/viperScores.txt
	sed -e 's/:/_/' viper-inferences/viperScores.txt | transpose_fast.pl > 1.tmp
	head -n 1 1.tmp > $@
	join.pl ${entrez2hugo} 1.tmp | cut -f 2- >> $@
	rm -f 1.tmp

tiedie/input: viper-inferences/viperInferences.txt
	mkdir -p $@
	mkdir -p tiedie/output
	cp ../shared_preprocessing/perturbation.matrix $@
	head -n 1 $< > $@/viperInferences.txt 
	join.pl ${pathway_genes} $< >> $@/viperInferences.txt
	head -n 1 $< | sed -e 's/	/	Type	/' > tiedie/input/sample.data
	grep -v  '^Drug' tiedie/input/perturbation.matrix | sed -e 's/	/	Genomic_Event	/' >> tiedie/input/sample.data
	grep -v '^id' tiedie/input/viperInferences.txt | sed -e 's/	/	TF_Activity	/' >> tiedie/input/sample.data

analysis/consensus-networks-kegg.txt: tiedie/input
	mkdir -p analysis 
	../scripts/consensusNetworks.py \
	--subtypes ${drug_classes} \
	-d tiedie/cytoscape-kegg \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperInferences.txt \
	-n ${pathway_net} \
	-c ${node_eq} \
	-z 6 > 1.tmp
	grep '^MTOR-MEK' 1.tmp | cut -f 2- > $@
	rm -f 1.tmp

analysis/net_sig_match_test.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_test.txt --output $@ 


networks-kegg: tiedie/output/linker_scores.txt tiedie/input/sample.data
	mkdir -p tiedie/output/networks-kegg
	../scripts/cutGraph \
	-n ${pathway_net} \
	-e $< \
	-d input/sample.data \
	-v --output tiedie/output/networks-kegg --step 1 --max_nodes 100 --simple 3.0 
	../scripts/sherpa.PSN \
	-n tiedie/output/networks-kegg \
	--undirected_types 'undirected:conflicted' \
	-d input/sample.data \
	--output $@ \
	-x 6 > out.err &



