
mkdir test-viper
../bin/run-viper.R -e data/test.data.tab -p data/phenotypes.tab -n data/multinet.adj -t Tumor -r Normal --permutations 100 --output test-viper

