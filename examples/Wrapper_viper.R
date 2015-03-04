
library(Biobase)
# viper libraries
library(mixtools)
library(viper)
# convert the names
library(GenomicFeatures)
library(mygene)
# Get command line options:
library(getopt)


GenSymToGenEnt <- function (GeneList) {
  GenEnt <- queryMany(GeneList, scopes="symbol", fields="entrezgenes", species="human", size=1)  
  GenEnt <- GenEnt[,"_id"]
}

GenEntToGenSym <- function (GeneList) {
  GenSym <- queryMany(GeneList, fields="symbol", species="human", size=1)  
  GenSym <- GenSym[,"symbol"] 
}

parse.tab <- function (file) {
  # The first two columns are the IDS and names, respectively. Note, the gene names 
  # might not be unique in the training data...
  m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
  return (as.matrix(m))
}

parse.phenotypes <- function (file) {
  m <- read.table(file, sep="\t", row.names=1, header=TRUE, quote="", check.names=FALSE)
  # hack to fix sample IDs and presevere the 'names' attribute. This needs to be fixed on
  # their side, most likely. 
  sampleID = m$sampleID
  attr(sampleID, 'names') <- rownames(m)
  Tissue = m$Tissue
  attr(Tissue, 'names') <- rownames(m)
  description = m$description
  attr(description, 'names') <- rownames(m)
  m = data.frame(I(sampleID), I(Tissue), I(description))
  metadata <- data.frame(labelDescription = c("TCGA Sample ID", "Tissue Type", "Sample Type"), row.names = c("sampleID", "Tissue", "description"))
  phenoData <- new("AnnotatedDataFrame", data = m, varMetadata = metadata)
  return (phenoData)
}

makeExpressionSet <- function(data.matrix, phenoData) { 
  exprSet <- new("ExpressionSet", exprs = data.matrix, phenoData = phenoData, 
                 annotation = "viper_input")
  return (exprSet)
}

# get inputs -------------------------------


# the real network score 
opt = getopt(matrix(c(
  'expression', 'e', 1, "character",
  'phenotypes', 'p', 1, "character",
  'output', 'o', 1, "character",
  'regulon', 'r', 1, "character"
),ncol=4,byrow=TRUE));

print(paste(opt$expression, " -> ", opt$output, sep=""))


# Load data ---------------------------------------------------------------

expdata <- parse.tab(opt$expression)

geall <- GenSymToGenEnt(row.names(expdata))
idx <- which(!is.na.data.frame(geall))

expdata <- expdata[idx,]
row.names(expdata) <- geall[idx]

pdata = parse.phenotypes(opt$phenotypes)
pdata@data <- pdata@data[1:ncol(expdata),]

expset.obj = makeExpressionSet(expdata, pdata)


# Run VIPER ---------------------------------------------------------------


load(opt$regulon)

vpres <- viper(expset.obj, regul)
print (vpres)

vpres.matrix= exprs(vpres);

row.names(vpres.matrix) = GenEntToGenSym(row.names(vpres.matrix))

write.table(format(vpres.matrix, digits=2), file=opt$output, col.names = NA, sep="\t", quote=F)

