

# Cell Line Specific Data Here
# EDIT THIS
DATA=HS578T_T24_INFchdir.tsv
name=HS578T


# static
targets=../shared_data/LJP56_drugs_targets.tsv
concentrations=../shared_data/drugs.conc.lst
pathway_maps=../shared_data/pathway_def_maps
node_eq=../shared_data/node_equivalence_classes/defs.txt
hugo2entrez=../shared_data/hugo2entrez.txt
entrez2hugo=../shared_data/entrez2hugo.txt
drug_classes=../shared_data/LJP56_Drugs_annotated.tsv
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

##
## Pre-TieDIE analysis
##
phenotypes.tab:
	../scripts/makePhenotypes.py ${targets} ${DATA} > $@

preprocessed/expression.entrez.tab:
	mkdir -p preprocessed
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
	sed -e 's/:/_/g' ${DATA} > tiedie/input/expression.tab
	../scripts/mergeSampleData.py tiedie/input/perturbation.matrix tiedie/input/viperInferences.txt > tiedie/input/sample.data

##
## Post-TieDIE analysis
##

cytoscape-kegg: tiedie/output/linker_scores.txt tiedie/input/sample.data
	mkdir -p tiedie/output/networks-kegg
	../scripts/cutGraph \
	-n ${pathway_net} \
	-e $< \
	-d tiedie/input/sample.data \
	-v --output tiedie/output/networks-kegg --step 1 --max_nodes 100 --simple 3.0 
	../scripts/sherpa.PSN \
	-n tiedie/output/networks-kegg \
	--undirected_types 'undirected:conflicted' \
	-d tiedie/input/sample.data \
	--output $@ \
	-x 6 > out.err &

analysis: cytoscape-kegg
	mkdir -p analysis 
	../scripts/consensusNetworks.py \
	--subtypes ${drug_classes} \
	-d $< \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperInferences.txt \
	-n ${pathway_net} \
	-c ${node_eq} \
	-o analysis

tfList.txt:
	echo 'name	TF' > $@
	cut -f 1 tiedie/input/viperInferences.txt | sed -e 's/$$/	1/' >> $@

analysis/net_sig_match_test.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_test.txt --output $@ 

##
## Get cytoscape view/layout 
## Import the network into 'session.cys' first!

layout:
	../../../../bin/copy-CY-layout.py ../shared_visualization/master_layout.cys --from_key consensus.union session.cys --to_key ${name}
	mv session.new_CY_layout.cys ${name}.cys
