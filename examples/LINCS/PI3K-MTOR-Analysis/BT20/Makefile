

# Cell Line Specific Data Here
# EDIT THIS
DATA_T24_DMSO=input/BT20_allINF_T24.tsv
DATA_T3_DMSO=input/BT20_allINF_T3.tsv
DATA_T3=input/BT20_RTK_MAPK_AKT_T3.tsv
DATA_T24=input/BT20_RTK_MAPK_AKT_T24.tsv
name=BT20

# static
targets=../shared_data/LJP56_Drugs_annotated.tsv
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
phenotypes.t3.tab:
	../scripts/makePhenotypes.py ${targets} ${DATA_T3} > $@

phenotypes.t24.tab:
	../scripts/makePhenotypes.py ${targets} ${DATA_T24} > $@

preprocessed/expression.entrez.t3.tab:
	join.pl ${hugo2entrez} ${DATA_T3} | cut -f 2- | grep -v '^	' | ../scripts/averageDupeProbes.py  > $@

preprocessed/expression.entrez.t3.DMSO.tab:
	join.pl ${hugo2entrez} ${DATA_T3_DMSO} | cut -f 2- | grep -v '^	' | ../scripts/averageDupeProbes.py  > $@
	head -n 1 preprocessed/expression.entrez.t3.DMSO.tab | transpose.pl | ../scripts/makeFeatureInput.py \
	| sed -e 's/Key.*/Sample	description/' > preprocessed/features.t3.tab

preprocessed/expression.entrez.t24.DMSO.tab:
	join.pl ${hugo2entrez} ${DATA_T24_DMSO} | cut -f 2- | grep -v '^	' | ../scripts/averageDupeProbes.py  > $@
	head -n 1 preprocessed/expression.entrez.t24.DMSO.tab | transpose.pl | ../scripts/makeFeatureInput.py \
	| sed -e 's/Key.*/Sample	description/' > preprocessed/features.t24.tab

#viper-inferences-t3-control: preprocessed/expression.entrez.t3.DMSO.tab preprocessed/features.t3.tab
#	mkdir -p $@
#	../../../../bin/run-viper-supervised.R -o $@ -n ${brca_regul}  -e $< \
#	-t Test -r DMSO -p preprocessed/features.t3.tab
#
preprocessed/expression.entrez.t24.tab:
	join.pl ${hugo2entrez} ${DATA_T24} | cut -f 2- | grep -v '^	' | ../scripts/averageDupeProbes.py  > $@

viper-inferences-t3: preprocessed/expression.entrez.t3.DMSO.tab
	mkdir -p $@
	../../../../bin/run-viper-unsupervised.R -o $@ -n ${brca_regul}  -e $<

viper-inferences-t24: preprocessed/expression.entrez.t24.DMSO.tab
	mkdir -p $@
	../../../../bin/run-viper-unsupervised.R -o $@ -n ${brca_regul}  -e $<

viper-inferences-t3/viperInferences.txt: viper-inferences-t3/viperScores.txt
	sed -e 's/:/_/' $< | transpose_fast.pl > 1.tmp
	head -n 1 1.tmp > data.tmp
	join.pl ${entrez2hugo} 1.tmp | cut -f 2- >> data.tmp
	../scripts/centerSamples.py < data.tmp > $@
	rm -f 1.tmp data.tmp

viper-inferences-t24/viperInferences.txt: viper-inferences-t24/viperScores.txt
	sed -e 's/:/_/' $< | transpose_fast.pl > 1.tmp
	head -n 1 1.tmp > data.tmp
	join.pl ${entrez2hugo} 1.tmp | cut -f 2- >> data.tmp
	../scripts/centerSamples.py < data.tmp > $@
	rm -f 1.tmp data.tmp

# sort inferences for the heatmap: order by drug class
sort_inferences:
	../scripts/sortMatrix_DrugClass.py ../shared_data/LJP56_Drugs_annotated.tsv	 < viper-inferences-t24/viperInferences.txt \
	> viper-inferences-t24/viperInferences.byDrugClass.cdt
	../scripts/sortMatrix_DrugClass.py ../shared_data/LJP56_Drugs_annotated.tsv	 < viper-inferences-t3/viperInferences.txt \
	> viper-inferences-t3/viperInferences.byDrugClass.cdt

tiedie/input: viper-inferences/viperInferences.txt
	mkdir -p $@
	mkdir -p tiedie/output
	cp ../shared_preprocessing/perturbation.matrix $@
	head -n 1 $< > $@/viperInferences.txt 
	join.pl ${pathway_genes} $< >> $@/viperInferences.txt
	../scripts/mergeSampleData.py tiedie/input/perturbation.matrix tiedie/input/viperInferences.txt > tiedie/input/sample.data

##
## Post-TieDIE analysis
##

networks-kegg: tiedie/output/linker_scores.txt tiedie/input/sample.data
	mkdir -p tiedie/output/networks-kegg
	../scripts/cutGraph \
	-n ${pathway_net} \
	-e $< \
	-d input/sample.data \
	-v --output tiedie/output/networks-kegg --step 1 --max_nodes 100 --simple 3.0 

cytoscape-kegg:
	../scripts/sherpa.PSN \
	-n tiedie/output/networks-kegg \
	--undirected_types 'undirected:conflicted' \
	-d tiedie/input/sample.data \
	--output $@ \
	-x 3 > out.err &

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

differential-viper-analysis: 
	mkdir -p $@ 
	../scripts/diffViperScores.py \
	--subtypes ${drug_classes} \
	-s tiedie/input/perturbation.matrix \
	-t tiedie/input/viperInferences.txt \
	--output $@ 

analysis/net_sig_match_test.REPORT:
	./scripts/net-signature.R --scores analysis/net_sig_match_test.txt --output $@ 


layout:
	../../../../bin/copy-CY-layout.py ../shared_visualization/master_layout.cys --from_key consensus.union session.cys --to_key ${name}
	mv session.new_CY_layout.cys ${name}.cys

heatmap.cdt:
	../scripts/makeViperHeatmap.py \
	--subtypes ${drug_classes} \
	-t tiedie/input/viperInferences.txt \
	> $@

