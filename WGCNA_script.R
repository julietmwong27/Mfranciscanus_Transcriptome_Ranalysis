# Juliet M Wong
# Weighted Gene Co-Expression Analysis (WGCNA) (Langfelder and Horvath 2008)
# Using limma (Ritchie et al. 2015) - voom-corrected data to determine clustering of similarly expressed genes and correlating them to treatments

# Wong et al. 2019 - Transcriptional profiles of early stage red sea urchins (Mesocentrotus
#franciscanus) reveal differential regulation of gene expression across
#development

# Download WGCNA
# Only run the following commands once to install WGCNA and flashClust on your computer
install.packages("BiocManager")
BiocManager::install("WGCNA")
install.packages("flashClust")

# Load WGCNA and flashClust libraries every time you open R
library(WGCNA)
library(flashClust)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# formatting data for WGCNA
# Take an object that contains the normalized counts file 
head(v2$E)
datExpr <- v2$E
head(datExpr)
dim(datExpr)

# Manipulate file so it matches the format WGCNA needs
datExpr = as.data.frame(t(datExpr))
dim(datExpr)
names(datExpr)


# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg#allOK
# The last statement returned TRUE, so all genes have passed the cuts

# Cluster the samples to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main =2)
# Plot a line to show the cut
abline(h = 15, col = "red")

# No cluster is an outlier (below the line)

#Create an object called "datTraits" that contains your trait data
datTraits = read.csv("Traits.csv")
head(datTraits)
row.names(datExpr) = datTraits$Sample
#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) = datTraits$Sample
datTraits$Sample = NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order
collectGarbage()

# now expression data is datExpr, and the clinical traits is datTraits
# visualize how the clinical traits relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
#Plot the sample dendrogram and the colors underneath
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "SamplesAndTraits.RData")
#load("SamplesAndTraits.RData")



#OLD # Cluster samples by expression to see if there are any obvious outliers
A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity 
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits, signed = FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor, traitColors)
plotDendroAndColors(sampleTree, groupLabels=names(datColors), colors=datColors, main="SampleDendrogram and Trait Heatmap")

# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor, traitColors)
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, main = "Sample Dendrogram and Trait Heatmap")
# No visible outliers

# Contruct a gene co-expression matrix and generate modules
# Do not run in Rstudio enableWGCNAThreads()
# allowWGCNAThreads()
lnames = load(file = "SamplesAndTraits.RData")
# The variable lnames containes the names of loaded variables
lnames

###1-STEP NETWORK CONSTRUCTION
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Choosing 6 as the soft power value because it the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value
help("blockwiseModules")
dynamicMergeCut(21) # calculates the threshold for module merging using the number of samples (in this case, 24)
# 0.2748987
# Constructing the gene network and identifying modules
net = blockwiseModules(datExpr, power = 6, TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.2748987, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMS = TRUE, saveTOMFileBase = "stageTOM", verbose = 3)
# See how many modules were identified and what their sizes are
table(net$colors)
# Open a graphics window
sizeGrWindow(12,9)
#Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
help("plotDendroAndColors")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Save the module assignment and module eigengene information necessary for subsequent analysis
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "Stage_FIESTA_networkConstruction.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs (MEs0)
moduleTraitCor = cor(MEs, datTraits, use ="p")
moduleTraitPValue = corPvalueStudent(moduleTraitCor, nSamples)
# Color code each association by the correlation value
# Display correlations and their p-values
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPValue, 1), ")", sep = "")
dim(textMatrix) = dim (moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
#Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))

# Heirarchical clustering tree of modules
# Caluclate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
MEtree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7,6)
plot(MEtree, main = "Clustering of module eigengenes",xlab = "", sub = "")


### END - 1-step network construction

# Gene relationship to trait and important modules: Gene Significance (GS) and Module Membership (MM)
# Quantify associations of individual genes with our trait of interest (any stage) by defining Gene Significance GS as the absolute value of the correlation between the gene and trait
# For each module, we also define a quantitative measure of module membership MM as the correlationo f the module eigengene and the gene expression profile
# This allows us to quantify the similarity of all genes on the array at every module
# Define variable egg containing the egg column of datTrait
X8cell = as.data.frame(datTraits$X8cell)
names(X8cell) = "8cell"
# names (colors) of the modules
modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, X8cell, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(X8cell), sep = "")
names(GSPvalue) = paste("p.GS.", names(X8cell), sep="")

# Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for egg as well as high module membership in interesting modules
# Look at yellow module, which has a high negative association with egg
# Plot a scatterplot of Gene Significance vs. Module Membership in the yellow module
module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes,1]), xlab = paste("Module Membership in", module, "module"), ylab = "Gene significane for 8 cell", main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Summary output of network analysis results
# merge statistical info with gene annotation and write out file that summarizes the most important results for Excel
names(datExpr)


### START Step-by-step module construction
softPower = 6
adjacency = adjacency(datExpr, power = softPower) #specify network type

# Construct Networks: This may need to be run on a cluster (generated RData files can be transferred back to Rstudio)
# translate the adjacency into topological overlap matrix and caluclate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Generate Modules
# Generate a clustered gene tree; call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method="average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Plot the module assignment under the gene dendrogram
# Convert numeric labels into colors
dynamicColors= labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHand = 0.05, main = "Gene dendrogram and module colors")
# Merging of modules whose expression profiles are very similar
# To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors= dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hClust(as.dist(MEDiss), method= "average")
# Plot the result; how the eigengenes cluster together
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
# set a threshold for merging modules
dynamicMergeCut(21) # calculates the threshold for module merging using the number of samples (in this case, 24)
#MEDissThres = 0.2559876
# Plot the cut line into the dendrogram; adjust threshold if necessary
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight=MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs
# To see what the merging did to our module colors, we plat the gene dendrogram again, with both original and merged modules
sizeGrWindow(12,9)
# pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang=0.05)
dev.off()
# use the merged module colors in mergedColors; save relevant variables for next analysis
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEss
# Save module colors and labels for use in subsequent part
save(MEs, moduleLabels, moduleColors, geneTree, file= "FIESTA-networkConstruction-stepBystep.RData")

### END STEP-BY-STEP

# Relate gene expression modules to traits
# Correlate traits
# Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Print correlation heatmap between modules and traits
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) = dim(moduleTraitCor)
dev.off()
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values with a heatmap plot
#labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= c("NL","UL","NH","UH"),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.8,
               zlim= c(-1,1),
               main= paste("Module-Treatment relationships"))
dev.off()


# Obtain lists of genes for each module eigengene - for GO enrichment analyses 
# Grab names in the dataset and set them as probes
probes = names(datExpr)
#select the modules of interest
intModules = c("brown","red","tan","green","salmon","purple","pink","cyan","magenta","yellow","midnightblue","black","turquoise","blue","greenyellow","grey")
intModules
# Write out the gene lists for each module of interest
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their gene IDs
  modLLIDs = probes[modGenes];
  # Write them into a file
  fileName = paste("GeneList-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Bar graph of number of genes within each module eigengene
counts <- c(2912,2166,169,2613,124,294,908,97,340,2745,55,1216,17610,3562,229,86)
counts
wgcnacolors <- c("brown","red","tan","green","salmon","purple","pink","cyan","magenta","yellow","midnightblue","black","turquoise","blue","greenyellow","grey")
wgcnacolors
barplot(counts, col=wgcnacolors, xlab="Module Eigengenes")

