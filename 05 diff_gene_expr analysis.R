##############################
######## Fifth Step ##########
#### DGE analysis and pix ####
##############################

###  set up the data matrix ####

rm(list = ls())
options(stringAsFactors = F)
a=read.table('count.txt', header = T)
meta=a[,1:6]  # the basic info for corresponding gene
exprSet=a[,7:ncol(a)]  # the expression matrix

###  correlation analysis coeffiency ####

library(corrplot)
png('cor.png')
corrplot(cor(log2(exprSet+1)))
dev.off()

library(pheatmap)
png('pheatmap.png')
pheatmap(scale(cor(log2(exprSet+1))))
dev.off()

### hclust to clust the different group #####

## clairify the group_list: treat untreat/ tumor normal
colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
hc=hclust(dist(t(log2(exprSet+1))))
par(mar=c(5,5,5,10))
png('hclust.png',res=120)
plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
dev.off()

######################
## PCA
library(ggfortify)
df=as.data.frame(t(exprSet))
df$group=group_list
png('pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()
######################


#####################################
######## DGE analysis  ##############
#####################################

library(DESeq2)
library(edgeR)
library(limma)
library(airway)
library(pheatmap)

#####################################
#######      DEseq2      ############
#####################################
suppressMessages(library(DESeq2)) 
group_list=c(rep('IMP',3),rep('NC',3),rep('PC',3))
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
dds <- DESeq(dds)
#png("qc_dispersions.png", 1000, 1000, pointsize=20)
#plotDispEsts(dds, main="Dispersion plot")
#dev.off()

res <- results(dds, contrast=c("group_list","IMP","PC"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG=as.data.frame(resOrdered)
DEG=na.omit(DEG)
write.csv(DEG,"../pix/IMP_PC_na_results.csv")


###############################
########     edgeR     ########
###############################

d <- DGEList(counts=exprSet,group=factor(group_list))
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples
dge=d

design <- model.matrix(~0+factor(group_list))
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))

dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)

lrt <- glmLRT(fit,  contrast=c(-1,1,0))
nrDEG=topTags(lrt, n=nrow(exprSet))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
write.csv(nrDEG,"DEG_treat_12_edgeR.csv",quote = F)

lrt <- glmLRT(fit, contrast=c(-1,0,1) )
nrDEG=topTags(lrt, n=nrow(exprSet))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
write.csv(nrDEG,"DEG_treat_2_edgeR.csv",quote = F)
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

################################
########    limma/voom    ######
################################

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)

dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('treat_12-control','treat_2-control'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='treat_12-control', n=Inf)
DEG_treat_12_limma_voom = na.omit(tempOutput)
write.csv(DEG_treat_12_limma_voom,"DEG_treat_12_limma_voom.csv",quote = F)

tempOutput = topTable(fit2, coef='treat_2-control', n=Inf)
DEG_treat_2_limma_voom = na.omit(tempOutput)
write.csv(DEG_treat_2_limma_voom,"DEG_treat_2_limma_voom.csv",quote = F)

png("limma_voom_RAWvsNORM.png",height = 600,width = 600)
exprSet_new=v$E
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(as.matrix(exprSet))
hist(exprSet_new)
dev.off()







