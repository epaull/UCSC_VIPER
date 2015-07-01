

# static
concentrations=../shared_preprocessing/drugs.conc.lst
pathway_maps=../shared_preprocessing/pathway_def_maps
targets=../shared_preprocessing/LJP56_drugs_targets.tsv

clean:
	rm -f perturbation.matrix

perturbation.matrix:
	 ../scripts/makePerturbationMatrix.py < ${concentrations} --ljp_targets ${targets} ${pathway_maps}/* | transpose.pl > $@

