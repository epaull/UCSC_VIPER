

# Cell Line Specific Data Here
# EDIT THIS
DATA=input/BT20_T24_INFchdir.tsv


# static
targets=../shared_data/LJP56_drugs_targets.tsv
concentrations=../shared_data/drugs.conc.lst
pathway_maps=../shared_data/pathway_def_maps
hugo2entrez=../shared_data/hugo2entrez.txt
brca_regul=../../pathway-priors/brca-sig-regulon.rda

all: phenotypes.tab 

clean:
	rm -f phenotypes.tab
	rm -f perturbation.matrix

phenotypes.tab:
	../scripts/makePhenotypes.py ${targets} ${DATA} > $@


preprocessed/expression.entrez.tab:
	join.pl ${hugo2entrez} ${DATA} | cut -f 2- | grep -v '^	' | ../scripts/averageDupeProbes.py  > $@

viper-inferences: preprocessed/expression.entrez.tab
	mkdir -p $@
	../../../../bin/run-viper-unsupervised.R -o $@ -n ${brca_regul}  -e $<


