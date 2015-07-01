#!/usr/bin/env  Rscript
library(getopt)
library(graphite)

opt = getopt(matrix(c(
    'outdir', 'd', 1, "character"
    ),ncol=4,byrow=TRUE));

for(i in names(kegg)){
    edge.i <- edges(convertIdentifiers(kegg[[i]],type="symbol"))
    name.i <- gsub(" ","_",i)
    name.i <- gsub("/_","",name.i)
    name.i <- gsub(",","",name.i)
    write.table(edge.i, file=paste0(opt$outdir, "/", name.i),row.names=FALSE, sep="\t", quote=F)

