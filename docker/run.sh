

docker run --rm \
	-v /Users/evanpaull/Documents/GenINST/UCSC_VIPER/data:/data \
	-v /Users/evanpaull/Documents/GenINST/UCSC_VIPER/output:/output \
	ucsc-viper \
	-e /data/test.data.tab \
	-p /data/phenotypes.tab \
	-n /data/multinet.adj \
	-t Tumor -r Normal \
	--permutations 100 \
	--output /output
