
phenotypes.tab:
	./makePhenotypes.py > $@

# BT20_T24_INFchdir.tsv - expression set
marina:
	mkdir -p $@
	 ../../../../bin/run-marina.R -n ~/Documents/UCSC_VIPER/pathways/brca-tf-regulon.rda \
	-e expression.tab -p phenotypes.tab -t MEK	-r PI3K-MTOR --permutations 100 --output $@


perturbation.matrix:
	 ./makePerturbationMatrix.py < drugs.conc.lst ../input/mek.upstream.txt ../input/pi3k.upstream.txt > $@
