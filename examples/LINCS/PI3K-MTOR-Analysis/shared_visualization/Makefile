




consensus.union.txt:
	cut -f 1-3 ../BT20/analysis/consensus-networks-kegg.txt >> 1.tmp
	cut -f 1-3 ../HS578T/analysis/consensus-networks-kegg.txt >> 1.tmp
	cut -f 1-3 ../MCF7/analysis/consensus-networks-kegg.txt >> 1.tmp
	cut -f 1-3 ../MCF10A/analysis/consensus-networks-kegg.txt >> 1.tmp
	cut -f 1-3 ../MDAMB231/analysis/consensus-networks-kegg.txt >> 1.tmp
	cut -f 1-3 ../SKBR3/analysis/consensus-networks-kegg.txt >> 1.tmp
	sort -u 1.tmp > $@

