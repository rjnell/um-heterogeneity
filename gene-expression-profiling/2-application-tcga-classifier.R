###  
### STEP 2 - Apply TCGA signature to LUMC data
###

# Set seed
set.seed(83539)

# Load LUMC data
lumc_normalised = as.matrix(read.csv("data/normalised-80.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))

# Load classifier data
classifier = read.csv("../um-rna/gene-expression-profiling/classifier.tsv", stringsAsFactors = F, sep = "\t")

# Plot 
lumc_classifier = lumc_normalised[classifier$ENSG,]
lumc_classifier_ph = pheatmap::pheatmap(lumc_classifier,
                                        scale="row", 
                                        clustering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                                        cutree_cols = 2, 
                                        cutree_rows = 2, 
                                        color = RColorBrewer::brewer.pal(9,"RdBu"), 
                                        show_rownames = T, 
                                        fontsize_col = 5,
                                        fontsize_row = 2,
                                        breaks = seq(from=-7, to=7, length.out = 10))
lumc_classifier_cutree = cutree(lumc_classifier_ph$tree_col, k = 2)

m1 = rowMeans(lumc_classifier[,names(which(lumc_classifier_cutree == 1))])
m2 = rowMeans(lumc_classifier[,names(which(lumc_classifier_cutree == 2))])

# Calculate percentage cutree_1 equals "Class II"
higherin1 = names(which(m1>m2))
lowerin1 = names(which(m1<m2))
c1isII = (length(which(higherin1 %in% classifier$ENSG[which(classifier$class == "Class II higher")])) + length(which(lowerin1 %in% classifier$ENSG[which(classifier$class == "Class I higher")])))/200

# Calculate percentage cutree_2 equals "Class II"
higherin2 = names(which(m2>m1))
lowerin2 = names(which(m2<m1))
c2isII = (length(which(higherin2 %in% classifier$ENSG[which(classifier$class == "Class II higher")])) + length(which(lowerin2 %in% classifier$ENSG[which(classifier$class == "Class I higher")])))/200

# Determine which tumours are Class I and II
if (c1isII < c2isII) {
  lumc_classifier_cutree[which(lumc_classifier_cutree==1)] = "Class I"
  lumc_classifier_cutree[which(lumc_classifier_cutree==2)] = "Class II"
} else {
  lumc_classifier_cutree[which(lumc_classifier_cutree==1)] = "Class II"
  lumc_classifier_cutree[which(lumc_classifier_cutree==2)] = "Class I"
}
lumc_classifier_cutree
table(lumc_classifier_cutree)

# Save classification
write.table(lumc_classifier_cutree, "gene-expression-profiling/lumc_classes.tsv", sep="\t")

# Save as image
png("gene-expression-profiling/tcga-classifier.png",res=600,width=5000,height=5000)
pt = pheatmap::pheatmap(lumc_classifier, 
                        annotation_col = data.frame(lumc_classifier_cutree), 
                        scale="row", 
                        clustering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                        cutree_cols = 2, 
                        cutree_rows = 2, 
                        color = RColorBrewer::brewer.pal(9,"RdBu"), 
                        show_rownames = T, 
                        fontsize_col = 5,
                        fontsize_row = 2,
                        breaks = seq(from=-7, to=7, length.out = 10))
dev.off()
