

viper-output:
	mkdir viper-output
	../bin/run-viper.R -o viper-output -n multinet.adj -e data.tab -p phenotypes.tab -t Tumor -r Normal --num_results 25

viper-output-rda:
	mkdir $@
	../bin/run-viper.R -o $@ -n brca-sig-regulon.rda -e data.tab -p phenotypes.tab -t Tumor -r Normal --num_results 25 --permutations 100
