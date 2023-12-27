
# Change working directory.
# change this to your own directory
path  = "/Users/aminaw/Desktop/bio assignment in R"
file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"
# extract the files into folders.
untar(file_name)

# change directory to the extracted folders 
setwd(paste(getwd() , "/brca_tcga_pan_can_atlas_2018", sep = ""))


# data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, data_mutations.txt and data_cna.txt
clinical = read.delim("data_clinical_patient.txt")
rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")
# in this assignment we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

# This is more for simplicity.If you keep your analysis would still be correct so no worries.
keep = !duplicated(rnaseq[,1])
rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols
rownames(rnaseq)  = rnaseq[,1]

# Read CNA Data
cna = read.delim('data_cna.txt')

# find ERBB2 in cna
erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.
hist(as.numeric(cna[erbb2_indx,-c(1,2)]), main = "Histogram of ERBB2 CNA Levels")
# match patients in rnaseq to patients in cna.
rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.
rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.
sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

# Pre-allocate memory for ERBB2
meta_erbb2 = matrix(0,length(rna_cna_id),1)

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
}

print(meta_erbb2)

# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.


# simple checks to make sure. 
col_i = colnames(rna_cna_sub)[1]
col_cna = which(colnames(cna)==col_i)

# sanity check
(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.
pos_example = which(meta_erbb2==1)[1]
print(pos_example)
col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# We will add a title to the metadata.
colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers
rna_cna_sub = round(rna_cna_sub)
print(rna_cna_sub)

# Install DESeq2.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2
BiocManager::install("DESeq2")
library(DESeq2)

#normaliza data 
dds <- DESeqDataSetFromMatrix(countData = round(rna_cna_sub),
                              colData = meta_erbb2,
                              design = ~ERBB2Amp)

dds <- DESeq(dds)
View(dds)


#-------------Differential Expression Analysis------------------------

res <- results(dds)
print(res)
# using the Summary funcftion
summary(res)
rownames(res) = rnaseq[keep,1]


# Column Explanation
# log2FoldChange: Logarithmic two-fold difference, indicating the fold change of the expression level of the experimental group relative to the control group. upregularization is represented with positive valuses and negative values represents downregulation.
# lfcSE: standard error of log2FoldChange, a measure of the uncertainty of log2FoldChange
# stat: Statistical value, calculated by devinding log2FoldChange by the standard error, used to test whether log2FoldChange is significantly different from zero
# pvalue: Used to indicate if there is significant difference; we want to set the significant value <0.05 for our hypothesis test
# padj: The adjusted p-value (or false discovery rate, FDR) after multiple testing correction. It takes into account the increased chance of false positives when conducting multiple statistical tests.


# p-value <0.05
signif = which(res$padj<0.05)
print(signif)
deg = res[signif,]
print(deg)

#Sort the results to find the top ten differentially expression genes
logFC_initial <- res$log2FoldChange
abs_logFC <- abs(res$log2FoldChange)
print(abs_logFC)
sorted_index <- order(abs_logFC, decreasing = TRUE)
sorted_res <- res[sorted_index, ]
print(sorted_index)
# select top ten sorted data
# taking the absolute value of log2FoldChange and then the top 10 
top_ten <- sorted_res[1:10, ]
print(top_ten)

# Separate deg
dup = deg[deg[,2]>0.,] # up regulated
ddown = deg[deg[,2]<0.,] # down regulated

# plot a Pie chart to show the percentage of dup and ddown in deg
up_ratio <- nrow(dup)/nrow(deg)
down_ratio <- nrow(ddown)/nrow(deg)
label <- c(paste('dup:',round(up_ratio * 100),"%"),paste('ddown:',round(down_ratio * 100),"%")) 
colors <- c("coral3", "skyblue")
pie(c(up_ratio, down_ratio),labels = label,col = colors, main = "Proportion of dup and ddown in deg")


# ---- HER+ vs HER- -----------

meta_erbb2 <- data.frame(ERBB2Amp = meta_erbb2) # Convert meta_erbb2 to data frame
# filter genes based on the p-value threshold
significant_genes <- res[complete.cases(res) & res$padj < 0.05, ]

# subset significant genes for the samples
her2_pos_genes <- significant_genes[significant_genes$log2FoldChange > 0, ] # HER2 positive samples
her2_neg_genes <- significant_genes[significant_genes$log2FoldChange < 0, ]# HER2 negative samples

# to obtain the counts of differentially expressed genes
num_her2_pos_genes <- nrow(her2_pos_genes)
num_her2_neg_genes <- nrow(her2_neg_genes)

# to finally get the counts of HER2+ and HER2- samples
num_her2_pos_samples <- sum(meta_erbb2$ERBB2Amp == 1, na.rm = TRUE)
num_her2_neg_samples <- sum(meta_erbb2$ERBB2Amp == 0, na.rm = TRUE)

#the results
cat("Number of differentially expressed genes from HER2+ vs HER2- samples:\\n")
cat("HER2+ samples:", num_her2_pos_samples, "\\n")
cat("HER2- samples:", num_her2_neg_samples, "\\n\\n")

cat("Number of differentially expressed genes:\\n")
cat("HER2+ vs HER2-:", num_her2_pos_genes + num_her2_neg_genes, "\\n")
cat("HER2+:", num_her2_pos_genes, "\\n")
cat("HER2-:", num_her2_neg_genes, "\\n")




# Heatmap
# Select the expression values of the significant genes
expression_values <- assay(dds)[rownames(significant_genes), ]
normalized_expression_values <- log2(expression_values + 1)

color_palette <- colorRampPalette(c("lightblue", "azure", "coral"))(100) # select the colours; red for upregulated genes, blue for downregulated genes

heatmap(normalized_expression_values, col = color_palette, scale = "none", main = "Differentially Expressed Genes Heatmap")


#------------------------For pathhway enrichment-------------------------------------------

# first we will need Entrez IDs
entrez_ids = rnaseq[keep,2]
entrez_all = entrez_ids[signif]
print(entrez_all)
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]
print(entrez_up)
print(entrez_down)


# Pathway Enrichment

BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
library(clusterProfiler)
library(pathview)
# Do a KEGG pathway over-representation analysis
all_paths = enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.5)


# plot different gene expression in hippo signaling pathway 
data <- c(logFC_initial)
df <- data.frame(entrezID = data)
pathview(gene.data = df, 
         pathway.id = "hsa04390",
         species = "hsa",
         out.suffix = "4")

print(all_paths)
dotplot(all_paths)
barplot(all_paths)

#the top 10 pathway enrichment results
top_10_paths <- head(all_paths, 10)
print(top_10_paths)



# -------------Principal Components Analysis---------------------------



# first transform the data to visualize
rld <- vst(dds, blind=FALSE)

pc = prcomp(assay(rld))

# extract first principal component 
loadings_pc1 <- pc$rotation[, 2]
# combine the results with initial variation
loadings_data <- data.frame(Variable = colnames(assay(rld)), Loadings_PC1 = loadings_pc1)
# sort by absolute value
loadings_data <- loadings_data[order(abs(loadings_data$Loadings_PC1), decreasing = TRUE), ]
print(head(loadings_data))

#PCA plot
plot(pc$rotation[,1], pc$rotation[,2], col = ifelse(meta_erbb2 == 1, "red", "darkblue"), pch = 19) # red identify amplified genes, black identify not amplified genes

#Clustering
dist_matrix <- dist(t(assay(rld)))
cluster_results <- hclust(dist_matrix)
clusters <- cutree(cluster_results, k = 11)  # K for the 11 columns
plot(pc$rotation[, 1], pc$rotation[, 2], col = clusters, pch = 19, main = "PCA Plot with Clustering")
