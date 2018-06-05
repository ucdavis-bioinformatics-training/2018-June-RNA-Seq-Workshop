# install edgeR and limma if not already installed
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR") # also installs limma as a dependency 

library(edgeR) # loads both edgeR and limma

setwd("C:/users/bpdurbin/desktop/course")

# read in counts table
counts <- read.delim("all_counts.txt", row.names = 1)
head(counts)

# Create object for use in limma and edgeR analysis
d <- DGEList(counts)

# calculate normalization factors
d <- calcNormFactors(d)
d

# filter genes below X cpms in all samples 
cutoff <- 1
drop <- which(apply(cpm(d), 1, max) < cutoff)
d <- d[-drop,] 
dim(d) # number of genes left

# Group information
cultivar <- rep(c("C","I5","I8"), each = 8)
time <- rep(rep(c("6","9"), each = 4),3)
cultivar
time

group <- interaction(cultivar, time)
group

# Quick MDS plot
plotMDS(d, col = as.numeric(group))


###################### Limma-voom analysis 
mm <- model.matrix(~0 + group) # specify model with no intercept for easier contrasts
y <- voom(d, mm, plot = T)
mtext(side = 3, line = 0.5, text = "1-Factor Model Without Intercept")
fit <- lmFit(y, mm)

head(coef(fit))

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp2, file = "I5_v_C_time6.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 20)
results.mod1 <- tmp2 # for illustration only

# Comparison between times 6 and 9 for cultivar I5	
contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp2, file = "time9_v_time6_I5.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 20)

##########################
# Other ways of fitting the exact same model (different parameterizations)
# Model with intercept
mm <- model.matrix(~group)
y <- voom(d, mm, plot = T)
mtext(side = 3, line = 0.5, text = "1-Factor Model With Intercept")
fit <- lmFit(y, mm)

head(coef(fit))

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(groupI5.6, levels = make.names(colnames(coef(fit)))) # C.6 is now reference group
rownames(contr) <- colnames(coef(fit))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
results.mod2 <- tmp2

# Two-factor model
mm <- model.matrix(~cultivar*time)
y <- voom(d, mm, plot = T)
mtext(side = 3, line = 0.5, text = "2-Factor Model")
fit <- lmFit(y, mm)

head(coef(fit))

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(cultivarI5, levels = make.names(colnames(coef(fit)))) # Cultivar C at time6 is reference group
rownames(contr) <- colnames(coef(fit))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
results.mod3 <- tmp2

head(results.mod1)
head(results.mod2)
head(results.mod3) # fits are identical

#################################################################################
# edgeR analysis
mm <- model.matrix(~0 + group) 
d <- estimateDisp(d, mm) 
d

fit <- glmFit(d, mm) 

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- glmLRT(fit, contrast = contr)
tmp2 <- topTags(tmp, sort.by = "P", n = Inf)$table
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","logCPM","PValue","FDR")]
length(which(tmp2$FDR < 0.05)) # number of DE genes
write.table(tmp2, file = "I5_v_C_time6_EDGER.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 20)

head(results.mod1, 20) # compare to limma results

#################################################################################
# edgeR quasi-likelihood analysis
mm <- model.matrix(~0 + group) 
d <- estimateDisp(d, mm) 

fit <- glmQLFit(d, mm)

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- glmQLFTest(fit, contrast = contr)
tmp2 <- topTags(tmp, sort.by = "P", n = Inf)$table
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","logCPM","PValue","FDR")]
length(which(tmp2$FDR < 0.05)) # number of DE genes
write.table(tmp2, file = "I5_v_C_time6_EDGERQL.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 20)

head(results.mod1, 20) # compare to limma results



