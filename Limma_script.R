# Juliet M Wong
# limma (Ritchie et al. 2015)
# Wong et al. 2019 - Transcriptional profiles of early stage red sea urchins (Mesocentrotus
#franciscanus) reveal differential regulation of gene expression across
#development

# Download limma and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install("edgeR") # analysis performed in limma, but uses a few of the functions in edgeR
library("edgeR")
BiocManager::install("limma")
BiocManager::install("statmod")
library("limma")
library("statmod")
library("ggplot2")

# Read in counts data (.matrix files)
GeneCounts <- read.delim("results_gene.matrix", header=TRUE, row.names=1)

# Only analyzing samples from this study (no including extra sequences used for the assembly)
All<-GeneCounts[ ,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)] 

# Filtering for more that 0.5 cpm across at least 3 samples
keep <- rowSums(cpm(All)>0.5) >=3
y<- All[keep,] 
# Check how many sequences remain after filtering
dim(y)

# Read in the model (contains Sample name, Parent treatment, Offspring treatment, Biological replicate #s)
model<-read.delim("model.txt", header=TRUE)
# Create a treatment factor (Stage)
f <- paste(model$Stage)
f <- factor(f)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f)
colnames(design) <- levels(f)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f)
dge <- calcNormFactors(dge)
dge$samples
# Create MDS plot (not yet voom-transformed)
colors <- rep(c("yellow","orange","green","red","blue","black","purple"))
plotMDS(dge, labels=1:21,col=colors[f])

# voom transformation using voomWithQualityWeights function
# combine observational-level weighting from voom with sample-specific quality weights
v1 <- voomWithQualityWeights(dge, design=design, lib.size=dge$samples$lib.size, plot = TRUE)
labelMDS <- rep(c("64cell","8cell","blastula","egg","gastrula","pluteus","prism"))
plotMDS(v1, labels=labelMDS[f],col=colors[f])

# Run voom again, this time adding a blocking variable and estimating correlation
# Samples are blocked by biological replicate
corfit <- duplicateCorrelation(v1,design,block=model$BioRep) 
corfit$consensus
v2 <- voomWithQualityWeights(dge,design,plot=TRUE,lib.size=dge$samples$lib.size,block=model$BioRep,correlation=corfit$consensus)
boxplot(v2$E, xlab = "", ylab = "Log2 counts per million", las=2, main="Voom transformed logCPM")

# MDS plot of voom-corrected data
par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(v2,pch=16,col=colors[f])
plot(mdsplot, pch=16, col=colors[f],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("egg","8 to 16 cell","64 cell to morula","blastula","gastrula","prism","pluteus"), col=c("red","orange","yellow","green","blue","purple","black"), pch=16, bty="n")
dev.off()

# PCA across voom-corrected data
pcavoom<- prcomp(t(na.omit(v2$E)))
summary(pcavoom)
pca <- as.data.frame(prcomp(t(na.omit(v2$E)))$x)
pca$f<-f
pcaplot <- qplot(x=PC1, y=PC2, data=pca,colour=pca$f, size=I(4))
pcaplot <- pcaplot + scale_colour_manual(values=c("yellow","orange","green","red","blue","black","purple"), name=NULL)
pcaplot <- pcaplot + theme_bw()
pcaplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.text=element_text(size=14),
               panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90), legend.key = element_blank())

# Perform Hierarchical Clustering on Principal Components using HCPC.R script

# Heatmap of top 500 maternal transcripts expressed at EG
maternal <- read.delim("maternal500.txt", header = FALSE)
full <- read.delim("voom-corrected counts.txt", header=TRUE, row.names=1)
dataSub <- full[rownames(full) %in% maternal$V1, ]
library("RColorBrewer")
library("gplots")
# Heatmap of all samples
heatmap.2(as.matrix(dataSub), col=brewer.pal(11,"RdBu"), trace="none",Colv = FALSE, scale="row", dendrogram = "row")
# Heatmap average expression by stage
maternal_ave <- read.delim("maternal500_average.txt", header = TRUE, row.names=1)
heatmap.2(as.matrix(maternal_ave),col=colorRampPalette(c("dodgerblue","white","firebrick3"))(200), trace="none",Colv = FALSE, scale="row", dendrogram = "row", density.info="none")

# Heatmap of transcripts with negative expression values at EG
low_maternal_neg <- read.delim("low_maternal_neg.txt", header = TRUE, row.names=1)
heatmap.2(as.matrix(low_maternal_neg),col=colorRampPalette(c("dodgerblue","white","firebrick3"))(200), trace="none",Colv = FALSE, scale="row", dendrogram = "none", density.info="none",labRow=FALSE)


