#source("http://bioconductor.org/biocLite.R")
#biocLite("RColorBrewer")
#biocLite("cluster")
#biocLite("WGCNA")
#biocLite("lattice")
#biocLite("genefilter")
#biocLite("dplyr"")
#biocLite("metaMA")
install.packages("metaMA")
library(edgeR)
library(ggplot2)
library(lattice)
library(genefilter)
library(RColorBrewer)
library(cluster)
library(WGCNA)

# read in expression counts data and meta data
data <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-June-RNA-Seq-Workshop/master/thursday/all_counts.txt", sep="\t", header=T, stringsAsFactors=F)
pdata <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-June-RNA-Seq-Workshop/master/friday/Figures/pdata.txt", sep="\t", header=T, stringsAsFactors=F)
pdata$Cultivar <- factor(pdata$Cultivar, levels=c("C", "I5", "I8"))
pdata$TimePoint <- factor(pdata$TimePoint, levels=c("6", "9"))

# normalized counts
keep <- rowSums(cpm(data) > 1) >= 1
counts <- data[keep,]
norm.counts <- cpm(counts)

# MDS plot
group <- interaction(pdata$Cultivar, pdata$TimePoint)
labs <- colnames(norm.counts)
plotMDS(norm.counts, col=as.numeric(group), labels=labs, lwd=2, cex=0.7)

# PCA plot
rv <- rowVars(norm.counts)
select <- order(rv, decreasing=TRUE)[seq_len(100)]
pca <- prcomp(t(norm.counts[select,]))
fac <- factor(apply(pdata[,c("Cultivar", "TimePoint")], 1, paste, collapse=":"))
colours <- brewer.pal(nlevels(fac), "Paired")
pcafig <- xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1, aspect="iso",
	col=colours, main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)),
	rep=FALSE)))
print(pcafig)

## Draw colour keys non-default way
pcafig <- xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1, aspect="iso",
	col=colours, main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)),
	rep=FALSE, columns=3)))
print(pcafig)

# Clustering of samples
dist.matrix <- dist(t(norm.counts))
sampleTree <- hclust(dist.matrix)
colours <- data.frame(Cultivar=labels2colors(pdata$Cultivar), TimePoint=labels2colors(pdata$TimePoint))
plotDendroAndColors(sampleTree, colors=colours, groupLabels=c("Cultivar", "TimePoint"), colorHeight=0.1, autoColorHeight=FALSE)

# Visualize differential expression results by volcano plot
#biocLite("dplyr")
library(dplyr)
data <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-June-RNA-Seq-Workshop/master/friday/Figures/I5_v_C_time6.txt", sep="\t", header=T, stringsAsFactors=F)
logFDR <- -log10(data$adj.P.)
## plot -log10(adj.P.Val) ~ logFC
plot(data$logFC, logFDR, xlab="log2(Fold-Change)", ylab="-log10FDR")
## clear plot area
dev.off()
## label significant genes
data <- mutate(data, sig=ifelse(data$adj.P.Val<0.05, "FDR<0.05", "Not Significant"))
## plot volcano with significant genes in red
p <- ggplot(data, aes(logFC, -log10(adj.P.Val))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("red", "black"))
p
## add gene names to selected genes
p + geom_text(data=filter(data, adj.P.Val<0.001), aes(label=Gene))
p + geom_text(data=data[order(data$adj.P.Val, decreasing=FALSE)[1:10],], aes(label=Gene))
dev.off()

# Heatmaps
biocL
library(devtools)
library(gplots)
slt <- order(rv, decreasing=TRUE)[seq_len(20)]
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,7))
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,6), cexRow=0.8, cexCol=0.8)
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,7), cexRow=0.8, cexCol=0.8, ColSideColors=labels2colors(pdata$Cultivar))
rowcols <- rep(brewer.pal(4, 'Set1'), each=5)
names(rowcols) <- rownames(norm.counts[slt,])
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,7), cexRow=0.8, cexCol=0.8, ColSideColors=labels2colors(pdata$Cultivar), RowSideColors=rowcols)

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
rlab <- t(rowcols)
rownames(rlab) <- "GeneType"
clab <- cbind(labels2colors(pdata$Cultivar), labels2colors(pdata$TimePoint))
colnames(clab) <- c("Cultivar", "TimePoint")
heatmap.3(norm.counts[slt,], col=heat.colors, trace="none", cexRow=0.8, cexCol=0.8, ColSideColors=clab, RowSideColors=rlab, ColSideColorsSize=2, RowSideColorsSize=2)

# select genes based on differential expression analysis results
# first, load in your differential expression analysis results
tmp <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-June-RNA-Seq-Workshop/master/friday/Figures/I5_v_C_time6.txt", sep="\t", header=T, stringsAsFactors=F)
sel.genes <- tmp$Gene[1:10]
# then, match the names of your selected genes to the rownames of your counts table
index <- match(sel.genes, rownames(norm.counts))
# then, follow some steps from above to generate necessary colors and labels
rowcols <- rep(brewer.pal(5, 'Set1'), each=2)
names(rowcols) <- rownames(norm.counts[slt,])
rlab <- t(rowcols)
rownames(rlab) <- "GeneType"
clab <- cbind(labels2colors(pdata$Cultivar), labels2colors(pdata$TimePoint))
colnames(clab) <- c("Cultivar", "TimePoint")
heatmap.3(norm.counts[index,], col=heat.colors, trace="none", cexRow=0.8, cexCol=0.8, ColSideColors=clab, RowSideColors=rlab, ColSideColorsSize=2, RowSideColorsSize=2)

# Heatmap using log transformed data.
log.counts <- cpm(counts, log=TRUE)
rv <- rowVars(log.counts)
slt <- order(rv, decreasing=TRUE)[seq_len(20)]
heatmap.2(log.counts[slt,], col=heat.colors, trace="none", margin=c(3,7))


# Visulize pathway enrichment results
#biocLite("pathview")
library(pathview)
DE.paths <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-June-RNA-Seq-Workshop/master/friday/Figures/I5_v_C_time6_KEGG.txt", sep="\t", header=T, stringsAsFactors=F)
head(DE.paths, 1)
pid <- DE.paths$pathway.code[3]
DE.expr <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-June-RNA-Seq-Workshop/master/friday/Figures/I5_v_C_time6.txt", sep="\t", header=T, stringsAsFactors=F)
head(DE.expr, 2)
rownames(DE.expr) <- DE.expr$Gene
gene.data <- subset(DE.expr, select="logFC")
pv.out <- pathview(gene.data=gene.data, pathway.id=pid, species="ath", gene.idtype="KEGG", kegg.native=T)
head(pv.out)

