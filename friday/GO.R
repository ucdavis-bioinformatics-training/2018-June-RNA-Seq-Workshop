# install topGO and org.At.tair.db if not already installed
# source("https://bioconductor.org/biocLite.R")
# biocLite("topGO")
# biocLite("org.At.tair.db")

setwd("C:/users/bpdurbin/desktop/biocore/2017 course") # set to location of DE results

library(topGO)

infile <- "I5_v_C_time6.txt"
outfile <- "I5_v_C_time6_GO.txt"

tmp <- read.delim(infile)
geneList <- tmp$P.Value
names(geneList) <- tmp$Gene
# Create topGOData object
GOdata <- new("topGOdata",
	ontology = "BP",
	allGenes = geneList,
	geneSelectionFun = function(x)x,
	annot = annFUN.org, mapping = "org.At.tair.db")
# Kolmogorov-Smirnov testing
resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
tab <- GenTable(GOdata, KS = resultKS, topNodes = length(resultKS@score), numChar = 120)
write.table(tab, file = outfile, quote = F, sep = "\t", row.names = F)

##############################################################################
### Fisher's Exact Test example, with GO terms from Biomart
# install biomaRt if not already installed
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

setwd("C:/users/bpdurbin/desktop/biocore/2017 course") # set to location of DE results

library(topGO)
library(biomaRt)

infile <- "ensembl_example_input.txt"
outfile <- "ensembl_example_GO.txt"

# Pull GO terms out of biomart
mart <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=mart)
# get all genes used in analysis
mygenes <- read.delim(infile, stringsAsFactors = F)$Gene
# take off part of ensembl ID after dot
mygenes2 <- unlist(lapply(strsplit(mygenes, split = ".", fixed = T), function(x)x[1])) 
gene.go <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"), 
			filter = "ensembl_gene_id", values = mygenes2, mart = ensembl)
gene.go <- gene.go[which(gene.go$go_id != ""),] # take out genes without GO terms
head(gene.go)

# Create list with element for each gene, containing vectors with all terms for each gene
gene2GO <- tapply(gene.go$go_id, gene.go$ensembl_gene_id, function(x)x)
head(gene2GO)

pcutoff <- 0.05 # cutoff for defining significant genes
DE <- read.delim(infile, stringsAsFactors = F)
# Define gene list as 1's if adjP < cutoff, 0, otherwise
tmp <- ifelse(DE$adj.P.Val < pcutoff, 1, 0)
geneList <- tmp
# geneList needs names that match those for GO terms, so stripping off part after dot again
names(geneList) <- unlist(lapply(strsplit(DE$Gene, split = ".", fixed = T), function(x)x[1]))
head(geneList)
# Create topGOData object
GOdata <- new("topGOdata",
		ontology = "BP",
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
		annot = annFUN.gene2GO, gene2GO = gene2GO)
# Fisher testing
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
tab <- GenTable(GOdata, Fisher = resultFisher, topNodes = length(resultFisher@score),
				numChar = 120)
write.table(tab, file = outfile, quote = F, sep = "\t", row.names = F)









