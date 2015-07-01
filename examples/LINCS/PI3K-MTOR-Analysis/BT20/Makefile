

# Cell Line Specific Data Here
# EDIT THIS
DATA=input/BT20_T24_INFchdir.tsv


# static
targets=../shared_data/LJP56_drugs_targets.tsv
concentrations=../shared_data/drugs.conc.lst
pathway_maps=../shared_data/pathway_def_maps
hugo2entrez=../shared_data/hugo2entrez.txt
entrez2hugo=../shared_data/entrez2hugo.txt
brca_regul=../../../../pathways/brca-tf-regulon.rda

# static
pathway_genes=../../pathway-priors/signal-transduction.lst
pathway_net=../../pathway-priors/signal-transduction.sif

all: phenotypes.tab 

clean:
	rm -f phenotypes.tab
	rm -f perturbation.matrix
	rm -rf viper-inferences

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

tiedie/input:
	mkdir -p $@
	cp ../shared_preprocessing/perturbation.matrix $@
	join.pl ${pathway_genes} viper-inferences/viperScores.txt > $@/viperInferences.txt


