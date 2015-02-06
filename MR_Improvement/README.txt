MR_Improvement
==============
Charles Markello 1/30/2015

Documented here is the pipeline that was used to integrate DNase hypersensitivity 
data with Hi-C assay gene domain data to filter out genes in a gene network that 
lay within closed chromatin domain structures. The resultant gene network file is 
then run with RNAseq expression data in VIPER to compute lists of significant 
master regulators. Also includes steps to evaluate MR results and compare them 
with true positive gene set datasets via precision-recall analysis 
(compare_MR_results) or via hypergeometric tests (compute_hypergeoProb_MR_results).

1) Split Hi-C gene domain file into tab-delimited format:
	./split_HiC_file  lists_t_dh.tab > hicdomains.bed

2) Find and sort coverage of DNase domains (e.g. 
GSM753973_UW.Breast_vHMEC.ChromatinAccessibility.RM035.DS18406.bed) within Hi-C 
domains (e.g. hicdomains.bed):
	bedtools coverage -a GSM753973_UW.Breast_vHMEC.ChromatinAccessibility.RM035.DS18406.bed \
					  -b hicdomains.bed > dnaseCoverage.bed

	bedtools sort -i dnaseCoverage.bed > dnaseCoverage.sorted.bed

3) Filter gene network file (e.g. brca_rnaseq851_multinet.adj) by DNase coverage 
level threshold (e.g. 0.2) of Hi-C domain file and write to output file (e.g. 
brca_rnaseq851.filtered.open.0.2.adj:
	./calculateOpenDomains.R -i dnaseCoverage.sorted.bed \
							 -n brca_rnaseq851.adj \
							 -o brca_rnaseq851.filtered.open.0.2.adj \
							 -s 0.2

4) Make gene expression profile file on experimental dataset (e.g. directory 
containing breast tumor TCGA RNAseq gene quantification files) against normal 
control data (e.g. directory containing breast normal TCGA RNAseq gene 
quantification files). Write output expression file (e.g. exprData.tab) and 
phenotype file (e.g. phenotypes.tab) labelling experimental samples (e.g. tumor) and 
normal control samples (e.g. normal). Parse out only RAW count data (e.g. -n RAW) :
	./makeExpressionFile -i ALL_BRCA_RNAseq -j Normal_Breast_RNAseq \
						 -e exprData.tab -p phenotypes.tab \
						 -t tumor -r normal \
						 -n RAW
						 
5) Run RNAseq DESeq normalization on expression file and write to file (e.g. exprData.deseq.varNormalized.tab). Run normalization (e.g. based on DESeq variant normalization -n variance [recomended]):
	./normalize_expression_file.R -e exprData.tab \
								  -p phenotypes.tab \
								  -o exprData.deseq.varNormalized.tab \
								  -d deseq.out.tab \
								  -n variance
								  
6) Run MARINa algorithm using VIPER:
	mkdir viper-output_open
	./run-viper.R -y 180 -o viper-output_open.0.2 \
				  -n brca_rnaseq851.open.0.2.adj \
				  -e exprData.tab \
				  -p phenotypes.tab \
				  -t tumor \
				  -r normal
				  
7) Convert true positive gene list gene ids to gene symbol:
	./convert_reference_MR_genelist_names.R -i breast_cancer_reference_MR_genelist.txt \
											-o breast_cancer_reference_MR_genelist.geneSymbol.txt
											
8) Compare MR results against true positive data:
	./compare_MR_results -i breast_cancer_reference_MR_genelist.genesymbol.txt \
						 -j viper-output_open.0.2/masterRegulators.txt
						 
9) Run hypergeometric test on MR results against true positive data:
	./compute_hypergeoProb_MR_results -i breast_cancer_reference_MR_genelist.geneSymbol.txt \
									  -j viper-output_open.0.2/masterRegulators.txt
								  
				