---
title: "Plots in R"
author: "Stefania Giacomello"
output: 
  html_document:
    keep_md: true
---

This R Markdown report contains instructions on how to make plots in R.
Set your working directory as first thing:

setwd("your/path/to/working/directory")

####1. As first thing, load the necessary libraries and set the working directory.

```r
#Load required library
suppressWarnings(library(metaMA))
suppressMessages(library(lattice))
suppressMessages(suppressWarnings(library(genefilter)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressWarnings(library(cluster))
suppressMessages(suppressWarnings(library(WGCNA)))
suppressMessages(suppressWarnings(library(matrixStats)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressWarnings(library(devtools))
suppressMessages(library(gplots))
suppressMessages(suppressWarnings(library(pathview)))

#Set the working directory where you have downloaded the count table
#setwd("Documents/Conferences/UCDavis")
```

####2. Load the count table and the relative metafile.

```r
d <- read.table(file="all_counts.txt", sep="\t", header=T, stringsAsFactors=F)

m <- read.table(file="metafile.txt", sep="\t", header=T, stringsAsFactors=F)

#Inspect the count table and check its dimensions
head(d)
```

```
##           C61 C62 C63 C64 C91  C92 C93 C94 I561 I562 I563 I564 I591 I592
## AT1G01010 341 371 275 419 400  542 377 372  677  522  455  508  821  466
## AT1G01020 164  94 176 155 200  183 166 115  172  157  122  152  189  171
## AT1G03987   0   0   0   0   0    0   0   0    0    0    0    0    0    0
## AT1G01030  20  34  40  27  28   36  22  40   20    7   57   38   25   10
## AT1G01040 738 487 610 690 945 1033 836 908  857  821  770  751  848  607
## AT1G03993   1   0   0   0   0    0   0   1    3    0    1    1    1    0
##           I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  691  500  157  473  459  228  590  491  565  496
## AT1G01020  163  185   46  162  119   53  172  212  169  157
## AT1G03987    0    0    0    0    0    0    0    0    0    0
## AT1G01030   17   26   49   17   24   48   27   28   47   32
## AT1G01040  871  756  361  618  641  439  783  692  768  625
## AT1G03993    0    0    0    1    0    1    3    2    1    1
```

```r
dim(d)
```

```
## [1] 32833    24
```

How many genes and samples does the count table contain?

####3. Normalize the data using count-per-million.

```r
keep <- rowSums(cpm(d) > 1) > 1  ##What does this command do?
counts <- d[keep,]
norm.counts <- cpm(counts)

#Inspect the normalized count table and check its dimensions
head(norm.counts)
```

```
##                  C61        C62        C63       C64       C91        C92
## AT1G01010  23.778291  29.959781  22.560571  30.08792  30.13728  36.203972
## AT1G01020  11.435894   7.590888  14.438765  11.13038  15.06864  12.223850
## AT1G01030   1.394621   2.745640   3.281538   1.93884   2.10961   2.404692
## AT1G01040  51.461521  39.327260  50.043448  49.54813  71.19933  69.001297
## AT1G01050 105.642554  88.102753  93.195666  92.27442 111.88467 105.606051
## AT1G01060 199.918945 216.098045 252.596353 204.72714 415.36712 443.932839
##                 C93        C94       I561        I562       I563
## AT1G01010  28.56118  31.961613  46.189096  37.5523602  32.580994
## AT1G01020  12.57601   9.880606  11.734896  11.2944838   8.736003
## AT1G01030   1.66670   3.436733   1.364523   0.5035757   4.081575
## AT1G01040  63.33460  78.013830  58.469801  59.0622370  55.137067
## AT1G01050  91.28971  85.402805  90.467860 114.7433228  74.900483
## AT1G01060 575.46610 628.750227 166.403553 139.9940477 175.006482
##                 I564       I591        I592       I593       I594
## AT1G01010  36.510032  58.048466  36.5329508  47.177734  35.590816
## AT1G01020  10.924262  13.363167  13.4058682  11.128756  13.168602
## AT1G01030   2.731065   1.767615   0.7839689   1.160668   1.850722
## AT1G01040  53.974476  59.957490  47.5869124  59.467158  53.813314
## AT1G01050  78.841545 102.733765 111.6371717 103.913908 118.588599
## AT1G01060 178.812912 359.037892 264.0407262 322.119461 285.153619
##                 I861       I862       I863       I864       I891      I892
## AT1G01010  18.371591  34.025482  36.542401  22.393049  40.941044  31.60514
## AT1G01020   5.382759  11.653548   9.473956   5.205402  11.935355  13.64621
## AT1G01030   5.733809   1.222903   1.910714   4.714326   1.873573   1.80233
## AT1G01040  42.242958  44.456127  51.031980  43.116441  54.333623  44.54329
## AT1G01050  51.721295 102.364188  88.927804  65.705920 104.920100 105.95124
## AT1G01060 182.428731 117.110962 124.116783 153.117382 287.697577 324.41935
##                 I893       I894
## AT1G01010  44.181417  37.844782
## AT1G01020  13.215326  11.979094
## AT1G01030   3.675268   2.441599
## AT1G01040  60.055448  47.687478
## AT1G01050 102.672921  95.298656
## AT1G01060 384.886609 292.076264
```

```r
dim(norm.counts)
```

```
## [1] 20646    24
```

What are the main differences between the raw count table and the normilized one?

####4. Investigate the data using Principal Component Analysis (PCA).

```r
rv <- rowVars(norm.counts)
select <- order(rv, decreasing=TRUE)[seq_len(100)]
pca <- prcomp(t(norm.counts[select,]))
fac <- factor(apply(m[,c("Cultivar", "TimePoint")], 1, paste, collapse=":"))
colours <- brewer.pal(nlevels(fac), "Paired")
pcafig <- xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1, aspect="iso",
    col=colours, main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)),
    rep=FALSE)))
print(pcafig)
```

![](script_plots_files/figure-html/PCA-1.png)<!-- -->

```r
#We can play around with visualization optsions:
pcafig <- xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=18, cex=1, aspect="iso",
    col=colours, main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)),
    rep=FALSE)))
print(pcafig)
```

![](script_plots_files/figure-html/PCA-2.png)<!-- -->

```r
pcafig <- xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=6, cex=1, aspect="iso",
    col=colours, main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)),
    rep=FALSE)))
print(pcafig)
```

![](script_plots_files/figure-html/PCA-3.png)<!-- -->

```r
#We can draw colour keys in a non-default way
pcafig <- xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1, aspect="iso",
    col=colours, main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)),
    rep=FALSE, columns=3)))
print(pcafig)
```

![](script_plots_files/figure-html/PCA-4.png)<!-- -->

####5. Cluster the data.

```r
dist.matrix <- dist(t(norm.counts))
sampleTree <- hclust(dist.matrix)
colours <- data.frame(Cultivar=labels2colors(m$Cultivar), TimePoint=labels2colors(m$TimePoint))
plotDendroAndColors(sampleTree, colors=colours, groupLabels=c("Cultivar", "TimePoint"), colorHeight=0.1, autoColorHeight=FALSE)
```

![](script_plots_files/figure-html/clustering-1.png)<!-- -->

####6. Visualize the differential gene expression analysis results by using volcano plots.

```r
#We need to import a new table with logFC values.
fc <- read.table(file="I5_v_C_time6.txt", sep="\t", header=T, stringsAsFactors=F)
logFDR <- -log10(fc$adj.P.Val)
plot(fc$logFC, logFDR, xlab="log2(Fold-Change)", ylab="-log10FDR")
```

![](script_plots_files/figure-html/volcano_plots-1.png)<!-- -->

```r
#Label significant genes
fc <- mutate(fc, sig=ifelse(fc$adj.P.Val<0.05, "FDR<0.05", "Not Significant")) #What does the command do?
```

```
## Warning: package 'bindrcpp' was built under R version 3.4.4
```

```r
#Make plot
p <- ggplot(fc, aes(logFC, -log10(adj.P.Val))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("red", "black"))

p + geom_text(data=filter(fc, adj.P.Val<0.001), aes(label=Gene))
```

![](script_plots_files/figure-html/volcano_plots-2.png)<!-- -->

```r
p + geom_text(data=fc[order(fc$adj.P.Val, decreasing=FALSE)[1:10],], aes(label=Gene))
```

![](script_plots_files/figure-html/volcano_plots-3.png)<!-- -->

```r
p + geom_text(data=fc[order(fc$adj.P.Val, decreasing=FALSE)[1:10],], aes(label=Gene)) + geom_vline(xintercept=c(-2,2), colour="black") + geom_hline(yintercept=1.3, colour="black")
```

![](script_plots_files/figure-html/volcano_plots-4.png)<!-- -->

####7. Look at patterns using heatmaps.

```r
slt <- order(rv, decreasing=TRUE)[seq_len(20)]
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,7))
```

![](script_plots_files/figure-html/heatmaps-1.png)<!-- -->

```r
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,6), cexRow=0.8, cexCol=0.8)
```

![](script_plots_files/figure-html/heatmaps-2.png)<!-- -->

```r
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,7), cexRow=0.8, cexCol=0.8, ColSideColors=labels2colors(m$Cultivar))
```

![](script_plots_files/figure-html/heatmaps-3.png)<!-- -->

```r
rowcols <- rep(brewer.pal(4, 'Set1'), each=5)
names(rowcols) <- rownames(norm.counts[slt,])
heatmap.2(norm.counts[slt,], col=heat.colors, trace="none", margin=c(3,7), cexRow=0.8, cexCol=0.8, ColSideColors=labels2colors(m$Cultivar), RowSideColors=rowcols)
```

![](script_plots_files/figure-html/heatmaps-4.png)<!-- -->

```r
#Heatmaps with multiple side bars using heatmap.3()
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
```

```
## Warning in strptime(x, fmt, tz = "GMT"): unknown timezone 'zone/tz/2018c.
## 1.0/zoneinfo/Europe/Stockholm'
```

```
## SHA-1 hash of file is 015fc0457e61e3e93a903e69a24d96d2dac7b9fb
```

```r
## SHA-1 hash of file is 015fc0457e61e3e93a903e69a24d96d2dac7b9fb
rlab <- t(rowcols)
rownames(rlab) <- "GeneType"
clab <- cbind(labels2colors(m$Cultivar), labels2colors(m$TimePoint))
colnames(clab) <- c("Cultivar", "TimePoint")

# The plot will be saved to a pdf file because of the size of the figure
pdf("test_heatmap3.pdf")
heatmap.3(norm.counts[slt,], col=heat.colors, trace="none", cexRow=0.8, cexCol=0.8, ColSideColors=clab, RowSideColors=rlab, ColSideColorsSize=2, RowSideColorsSize=2, margin=c(5,5))
dev.off()
```

```
## quartz_off_screen
##                 2
```

```r
#select genes from the differential expression analysis results
# first, load in your differential expression analysis results
sel.genes <- fc$Gene[1:10]
# then, match the names of your selected genes to the rownames of your counts table
index <- match(sel.genes, rownames(norm.counts))
# then, follow some steps from above to generate necessary colors and labels
rowcols <- rep(brewer.pal(5, 'Set1'), each=2)
names(rowcols) <- rownames(norm.counts[index,])
rlab <- t(rowcols)
rownames(rlab) <- "GeneType"
clab <- cbind(labels2colors(m$Cultivar), labels2colors(m$TimePoint))
colnames(clab) <- c("Cultivar", "TimePoint")
#Using log transformed data.
log.counts <- cpm(counts, log=TRUE)
rv <- rowVars(log.counts)
slt <- order(rv, decreasing=TRUE)[seq_len(20)]
# use non-default color scheme
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(log.counts[slt,], col=morecols, trace="none", margin=c(3,7))
```

![](script_plots_files/figure-html/heatmaps-5.png)<!-- -->

####8. Visualize pathways.

```r
#Visulize pathway enrichment results using bioconductor package "pathview"
DE.paths <- read.table(file="I5_v_C_time6_KEGG.txt", sep="\t", header=T, stringsAsFactors=F)
head(DE.paths, 1)
```

```
##   pathway.code                                  pathway.name      p.value
## 1     ath03010 Ribosome - Arabidopsis thaliana (thale cress) 5.378589e-35
##   Annotated
## 1       301
```

```r
pid <- DE.paths$pathway.code[3]

head(fc, 2)
```

```
##        Gene      logFC   AveExpr      P.Value    adj.P.Val      sig
## 1 AT4G12520 -10.254556 0.3581132 2.206726e-10 4.651998e-06 FDR<0.05
## 2 AT3G30720   5.817438 3.3950689 9.108689e-10 9.601014e-06 FDR<0.05
```

```r
rownames(fc) <- fc$Gene
colnames(fc)
```

```
## [1] "Gene"      "logFC"     "AveExpr"   "P.Value"   "adj.P.Val" "sig"
```

```r
## [1] "Gene"      "logFC"     "AveExpr"   "P.Value"   "adj.P.Val"
gene.data <- subset(fc, select="logFC")
head(gene.data)
```

```
##                logFC
## AT4G12520 -10.254556
## AT3G30720   5.817438
## AT5G26270   2.421030
## AT3G33528  -4.780814
## AT1G64795  -4.872595
## AT3G05955  -4.158939
```

```r
pv.out <- pathview(gene.data=gene.data, pathway.id=pid, species="ath", gene.idtype="KEGG", kegg.native=T)
```

```
## Info: Working in directory /Users/stefaniagiacomello/Documents/Conferences/UCDavis
```

```
## Info: Writing image file ath04141.pathview.png
```

```r
#By default, running pathview() will create an image file named by the pathway id (for example, in this case there should be a file named "ath04141.pathview.png" in the current directory).
#Another package for ploting is "piano". It generates network style graphs. "Cytoscape" is another possible software to generate graphs for enrichment analysis results.
#Visulize GO enrichment results use revigo.irb.hr web application. In a web browser (Safari, explorer, chrome, firefox), open revigo.irb.hr.
```
