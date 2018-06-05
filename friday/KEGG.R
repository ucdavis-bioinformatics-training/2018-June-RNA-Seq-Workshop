# install KEGGREST if not already installed
# source("https://bioconductor.org/biocLite.R")
# biocLite("KEGGREST")

library(KEGGREST)

setwd("C:/users/bpdurbin/desktop/biocore/2017 course") # set to where DE file lives

infile <- "I5_v_C_time6.txt"
outfile <- "I5_v_C_time6_KEGG.txt"

DE.table <- read.delim(infile, stringsAsFactors = F)
geneList <- DE.table$P.Value
names(geneList) <- DE.table$Gene
head(geneList)

# Pull all pathways for AT	
pathways.list <- keggList("pathway", "ath")

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list))	
genes.by.pathway <- sapply(pathway.codes,
	function(pwid){
		pw <- keggGet(pwid)
		if (is.null(pw[[1]]$GENE)) return(NA)
		pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
		pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
		return(pw2)
	}
)

# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
	function(pathway) {
		pathway.genes <- genes.by.pathway[[pathway]]
		list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
		list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
		scores.in.pathway <- geneList[list.genes.in.pathway]
		scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
		if (length(scores.in.pathway) > 0){
			p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
		} else{
			p.value <- NA
		}
		return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
	}
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
write.table(outdat, file = outfile, row.names = F, quote = F, sep = "\t")

