#!/usr/bin/env	RScript

library(ggplot2)
library(RColorBrewer)

library('getopt')

# the real network score --
opt = getopt(matrix(c(
    'data', 'd', 1, "character",
    'regulator', 'r', 1, "character",
    'output', 'o', 1, "character"
    ),ncol=4, byrow=TRUE));

png("viper-unsupervised/REST.png")
# sample ratio pval ratio
data = read.table('viper-unsupervised/REST.tab' , header=TRUE, sep="\t")
# track ratio scores, segregate by sample...
data$type <- factor(data$type, levels = data$type[order(data$viper)])
p <- ggplot(data, aes(factor(type), viper))
p + geom_violin(aes(fill=factor(type))) + 
geom_jitter(width = 0.1, size = 1.5, aes(colour = type)) +
scale_fill_brewer(palette='Set3') +
theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20))

