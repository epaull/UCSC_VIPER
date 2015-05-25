#!/usr/bin/env  Rscript

# Get command line options:
library('getopt')

# the real network score --
opt = getopt(matrix(c(
    'dm', 'm', 1, "character",
    'output', 'o', 1, "character"
    ),ncol=4, byrow=TRUE));


parse.tab <- function (file) {
  m <- as.matrix(read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE))
  return (m)
}


data <- parse.tab(opt$dm)

variances <- matrix(nrow=dim(data)[1], ncol=2)

for (i in 1:dim(data)[1]) {

        name = rownames(data)[i]
        var <- sd(as.numeric(data[i,]))^2
        variances[i,] <- c(name, var)
}

write.table(variances, file=opt$output,  sep="\t", quote=FALSE, col.names=NA)

