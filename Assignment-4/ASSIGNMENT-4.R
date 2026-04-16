install.packages("R.utils")
# Set your working directory
setwd("C:/Users/kdivy/OneDrive/Desktop/COMP SYS BIOLOGY/ASSIGNMENT-4")

# Extract the tar file
untar("GSE11352_RAW.tar", exdir = "GSE11352_RAW")

# Go into the folder
setwd("GSE11352_RAW")

# Unzip the .gz files
gz_files <- list.files(pattern = "*.CEL.gz")
sapply(gz_files, R.utils::gunzip)


library(affy)
library(limma)
setwd("C:/Users/kdivy/OneDrive/Desktop/COMP SYS BIOLOGY/ASSIGNMENT-4/GSE11352_RAW")
raw_data <- ReadAffy()
eset <- rma(raw_data)
exprs_data <- exprs(eset)

# Check number of samples (should be 18)
ncol(exprs_data)

# Assign sample groups based on GEO sample metadata
group <- factor(c(
  rep("Control_12h", 3),
  rep("Estrogen_12h", 3),
  rep("Control_24h", 3),
  rep("Estrogen_24h", 3),
  rep("Control_48h", 3),
  rep("Estrogen_48h", 3)
))

# Design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit model and apply contrasts for 12h and 24h comparisons
fit <- lmFit(exprs_data, design)

contrast.matrix <- makeContrasts(
  Estrogen_12h_vs_Control = Estrogen_12h - Control_12h,
  Estrogen_24h_vs_Control = Estrogen_24h - Control_24h,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Top DEGs - all genes
deg_12h <- topTable(fit2, coef = "Estrogen_12h_vs_Control", number = Inf)
deg_24h <- topTable(fit2, coef = "Estrogen_24h_vs_Control", number = Inf)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu133plus2.db")

library(hgu133plus2.db)
deg_12h$Gene.symbol <- mapIds(hgu133plus2.db,
                              keys = rownames(deg_12h),
                              column = "SYMBOL",
                              keytype = "PROBEID",
                              multiVals = "first")

deg_24h$Gene.symbol <- mapIds(hgu133plus2.db,
                              keys = rownames(deg_24h),
                              column = "SYMBOL",
                              keytype = "PROBEID",
                              multiVals = "first")
# Remove rows with missing gene symbols
rnk_12h <- deg_12h[!is.na(deg_12h$Gene.symbol), c("Gene.symbol", "t")]
rnk_24h <- deg_24h[!is.na(deg_24h$Gene.symbol), c("Gene.symbol", "t")]

# Remove duplicate gene symbols (optional: keep max |t| for each gene)
rnk_12h <- rnk_12h[!duplicated(rnk_12h$Gene.symbol), ]
rnk_24h <- rnk_24h[!duplicated(rnk_24h$Gene.symbol), ]

# Write to .rnk files
write.table(rnk_12h, "GSEA_12h.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rnk_24h, "GSEA_24h.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
