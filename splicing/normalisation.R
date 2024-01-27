### Analyze how samples differ in their usage of splice junctions
### Input: SJ.out.tab files derived from STAR alignment (GRCh38)
### Output: table + images *****

# Cohort 1 SJ.out.tab files
dir_cohort1 = "splicing/data/cohort1/"
splicing_files_cohort1 = list.files(path = dir_cohort1, pattern="SJ.out.tab", recursive = T)

# Cohort 2 SJ.out.tab files
dir_cohort2 = "splicing/data/cohort2/"
splicing_files_cohort2 = list.files(path = dir_cohort2, pattern="SJ.out.tab", recursive = T)

# Merge cohorts
merged_files = c(paste0(dir_cohort1, splicing_files_cohort1), paste0(dir_cohort2, splicing_files_cohort2))
merged_names = c(splicing_files_cohort1, splicing_files_cohort2)

# Initialize the list of all splicing_loci
all_splicing_loci = NULL

# Initialize the list of all sample_names
sample_names = NULL

# Iterate through all files
for (f in 1:length(merged_files)) {
  splicing_data = read.csv(merged_files[f], sep="\t", header = F, stringsAsFactors = F)
  splicing_loci = paste0(splicing_data$V1, "-", splicing_data$V2, "-", splicing_data$V3)
  all_splicing_loci = unique(c(all_splicing_loci, splicing_loci))
  sample_names = c(sample_names, strsplit(merged_names[f], "/")[[1]][1])
}

# Create a matrix
splicing_loci_usage = matrix(data = 0, nrow = length(all_splicing_loci), ncol = length(merged_files))
rownames(splicing_loci_usage) = all_splicing_loci

# Load sampleshet
samplesheet = read.csv("data/lumc-samplesheet.tsv", sep="\t", stringsAsFactors = F)

# Add sample names to matrix
colnames(splicing_loci_usage) = samplesheet$SAMPLE[match(sample_names, samplesheet$GS_ID)]

# Iterate through all files
for (f in 1:length(merged_files)) {
  splicing_data = read.csv(merged_files[f], sep="\t", header = F, stringsAsFactors = F)
  splicing_loci = paste0(splicing_data$V1, "-", splicing_data$V2, "-", splicing_data$V3)
  splicing_loci_usage[splicing_loci, f] = splicing_data$V7
}

splicing_loci_usage[1:10,1:20]

###

# PCA on unnormalized data
combined_pca = prcomp(t(splicing_loci_usage))
plot(combined_pca$x, col=c(rep("#BE1E2D",32), rep("#2EAADE",48)), pch=16)
text(combined_pca$x, labels=rownames(combined_pca$x), pos=2, cex=0.5)

# Load DESeq2
library(DESeq2)

# Build samples_table
samples_table = as.matrix(factor(c(rep(("cohort_1"),32),rep(("cohort_2"),48))),,drop=F)
rownames(samples_table) = colnames(splicing_loci_usage)
colnames(samples_table) = "cohort"

# Build DESeqData
DESeqData = DESeqDataSetFromMatrix(splicing_loci_usage,
                                   colData = samples_table,
                                   design = ~ cohort)

# Perform vst, a variance stabilizing transformation 
vsd = vst(DESeqData, blind=F)
combined_normalized = assay(vsd)
combined_normalized_genes = rownames(combined_normalized)

# Illustrate outcome
head(combined_normalized, 3)

# Build PCA of combined normalized data
combined_normalized_pca = prcomp(t(combined_normalized))
plot(combined_normalized_pca$x, col=c(rep("#BE1E2D",32), rep("#2EAADE",48)), pch=16)
text(combined_normalized_pca$x, labels=rownames(combined_normalized_pca$x), pos=4, cex=0.5)


write.table(combined_normalized, "data/splicing-normalised-80.tsv", row.names = T, sep="\t")
