

# Cell Line Specific Data Here
DATA=input/BT20_T24_INFchdir.tsv


# static
targets=../shared_data/LJP56_drugs_targets.tsv
concentrations=../shared_data/drugs.conc.lst
pathway_maps=../shared_data/pathway_def_maps

all: phenotypes.tab perturbation.matrix

clean:
	rm -f phenotypes.tab
	rm -f perturbation.matrix

phenotypes.tab:
	../scripts/makePhenotypes.py ${targets} ${DATA} > $@

