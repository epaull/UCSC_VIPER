
mkdir -p output
docker run --rm \
	-v `pwd`/data:/data \
	-v `pwd`/output:/output \
	ucsc-viper-unsupervised \
	-e /data/test.data.tab \
	-n /data/multinet.adj \
	--output /output
