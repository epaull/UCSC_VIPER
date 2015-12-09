

# preprocessing
# might need to do this:  ../../../bin/findDuplicateRows < data.tab
data.tab: ../RNA-Seq/TGCT_IlluminaHiSeq_RNASeqV2.tab
	transpose_fast.pl $< > 1.tmp
	join.pl phenotypes.tab 1.tmp | cut -f 1,3- > 2.tmp
	transpose_fast.pl 2.tmp > $@
	rm *.tmp

# filter low variance genes out
data.filtered.tab: data.tab
	head -n 1 $< > $@
	 ../../../bin/get-variance.R  -m data.tab -o variances.txt
	select.pl -gte 10 -k 3 variances.txt | cut -f 2 > 2.tmp
	join.pl 2.tmp data.tab >> $@

marina: data.filtered.tab phenotypes.tab
	mkdir $@
	../../../bin/run-marina.R -e $< -p phenotypes.tab \
	-n ../../../pathways/multinet.adj \
	-t High \
	-r Low \
	--permutations 1000 \
	--output $@

marina-aracne-plot: data.filtered.tab phenotypes.tab
	mkdir $@
	../../../bin/run-marina.R -e $< -p phenotypes.tab \
	-n data.DPI0.1.adj \
	-t High \
	-r Low \
	--num_results 85 \
	--permutations 1000 \
	--output $@
