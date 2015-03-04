
# convert the names
library(GenomicFeatures)
library(mygene)


GenEntToGenSym <- function (GeneList) {
  GenSym <- queryMany(GeneList, fields="symbol", species="human", size=1)  
  GenSym <- GenSym[,"symbol"] 
}

regulon_name <- 'brca-tf-regulon.rda'
load(regulon_name)

reg_names <- names(regul)
reg_genes <- GenEntToGenSym(reg_names)


net <- c()
for (i in 1:length(regul)) {
  target <- names(regul[[i]]$tfmode)
  vals <- regul[[i]]$tfmode
  for (j in 1:length(regul[i])) {
    net = rbind(net, c(reg_genes[i], vals[j], target[j]))
}
}

net[,3] <- GenEntToGenSym(net[,3])



write.table(net, 
            file=paste(unlist(strsplit(regulon_name, split='.', fixed=TRUE))[1], '.tsv', sep=''), 
            col.names = FALSE, row.names = FALSE, sep="\t", quote=F)
