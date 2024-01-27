###  
### STEP 1 - DeSEQ2 on LUMC RNA-seq
###

# Set seed
set.seed(83539)

# Load fragments_per_gene data from cohort 1
cohort_1_fpg_path = "data/cohort-1.fragments_per_gene.tsv"
cohort_1_fpg = read.csv(cohort_1_fpg_path, stringsAsFactors = F, check.names = F, sep = "\t")

# Plot PCA of raw data cohort 1
cohort_1_pca = prcomp(t(cohort_1_fpg[,2:ncol(cohort_1_fpg)]))
plot(cohort_1_pca$x, pch=16)

# Load fragments_per_gene data from cohort 2
cohort_2_fpg_path = "data/cohort-2.fragments_per_gene.tsv"
cohort_2_fpg = read.csv(cohort_2_fpg_path, stringsAsFactors = F, check.names = F, sep = "\t")

# Plot PCA of raw data cohort 2
cohort_2_pca = prcomp(t(cohort_2_fpg[,2:ncol(cohort_2_fpg)]))
plot(cohort_2_pca$x, pch=16)

# Check if features are the same for both cohorts
(length(which(cohort_1_fpg[,1] == cohort_2_fpg[,1])) == length(cohort_1_fpg[,1]))

# Combine raw data
combined_fpg = cbind(cohort_1_fpg[,2:ncol(cohort_1_fpg)], cohort_2_fpg[,2:ncol(cohort_2_fpg)])
rownames(combined_fpg) = cohort_1_fpg[,1]
colnames(combined_fpg)

# Select samples to include
selection_fpg = combined_fpg[,c(1:12,14,16,18:41,43:83,85)]
colnames(selection_fpg)

# Plot PCA of combined raw data
selection_pca = prcomp(t(selection_fpg))
plot(selection_pca$x, col=c(rep("#BE1E2D",32), rep("#2EAADE",48)), pch=16)
text(selection_pca$x, labels=rownames(selection_pca$x), pos=2, cex=0.5)

# Load samplesheet
samplesheet = read.csv("data/lumc-samplesheet.tsv", sep="\t", stringsAsFactors = F)

# Replace sample_IDs by sample_names
names = samplesheet$SAMPLE[match(rownames(selection_pca$x), samplesheet$GS_ID)]
colnames(selection_fpg) = names

# Install DESeq2
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# install.packages("colorspace")
# install.packages("tibble")

# Load DESeq2
library(DESeq2)

# Build samples_table
samples_table = as.matrix(factor(c(rep(("cohort_1"),32),rep(("cohort_2"),48))), ,drop=F)
rownames(samples_table) = names
colnames(samples_table) = "cohort"

# Build DESeqData
dds = DESeqDataSetFromMatrix(selection_fpg,
                             colData = samples_table,
                             design = ~ cohort)

# Only keep ENSG genes
keep = which(substr(rownames(dds),1,4) == "ENSG")
dds = dds[keep,]

# Normalise and save
lumc_normalised = assay(vst(dds, blind=F))
write.table(lumc_normalised, "data/normalised-80.tsv", sep="\t")
