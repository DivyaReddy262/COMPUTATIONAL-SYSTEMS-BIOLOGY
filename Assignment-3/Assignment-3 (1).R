# Install and load required package
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

setwd("C:/Users/kdivy/OneDrive/Desktop/COMP SYS BIOLOGY/ASSIGNMENT-3")
# Read expression data with duplicated gene names fixed
expr_raw = read.csv("HW03_expression.csv", check.names = FALSE)
rownames(expr_raw) = make.unique(as.character(expr_raw[, 1]))  # Fix duplicate names
exprData = expr_raw[, -1]  # Drop the name column
datExpr0 = as.data.frame(t(exprData))  # Transpose for WGCNA
# Load trait data
datTraits = read.csv("HW03_Traits.csv", row.names = 1)

# Match sample names in expression and trait data
samplesExpr = rownames(datExpr0)
samplesTraits = rownames(datTraits)
commonSamples = intersect(samplesExpr, samplesTraits)

# Subset to common samples
datExpr = datExpr0[commonSamples, ]
datTraits = datTraits[commonSamples, ]

# Check for missing samples/genes
gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

powers = c(1:20)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plotting
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.9, col="blue")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")

softPower = 6  # Choose based on scale-free fit above
adjacency = adjacency(datExpr, power = softPower)


# Turn adjacency into TOM
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
# Hierarchical clustering
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
# Dynamic tree cut to detect modules
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")


MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
# Plot module eigengene dendrogram
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# Merge modules with correlation > 0.75 (height = 0.25)
MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
# Plot updated module colors
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Save final colors
moduleColors = mergedColors


# Recalculate eigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# Correlate eigengenes with traits
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))
# Text for heatmap
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
# Heatmap
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = "Module–trait relationships")
