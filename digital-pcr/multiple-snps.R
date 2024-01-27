###
### COMPLETE CODE TO PROCESS COPY NUMBER DATA
###

# Load libraries
library(readxl)
library(digitalPCRsimulations)

# Load data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)
ids = lumc_data$Tumor_ID
lumc_ids = lumc_data$LUMC_ID
names(lumc_ids) = ids
references = lumc_data$Stable_reference
names(references) = ids

# Function to process raw data from SNP-based experiments
# If multiple SNPs are available, the SNP with smallest CI is chosen
process_snp_data_extra = function(ids, snp_data) {
  
  # Count ids
  n = length(ids)
  
  # Create matrix with final results
  results = matrix(data=NA, nrow=n, ncol=7)
  rownames(results) = ids
  
  # Iterate through ids to collect and process available data
  for (i in 1:n) {
    
    # Load id
    id = rownames(results)[i]
    
    # Select rows from snp_data (id match + OK mark)
    snp_data_selection = snp_data[which(snp_data$`File name` == id & snp_data$Mark == "OK"), , drop=F]
    
    # If no rows available: no data, no results to add
    if (nrow(snp_data_selection) > 0) {
      
      # If 1 row available: only data, to add
      if (nrow(snp_data_selection) == 1) {
        selected_row = 1  
      }
      # If more rows available: add data with smallest CI
      else if (nrow(snp_data_selection) > 1) {
        confidence_interval_lengths = as.numeric(snp_data_selection$`99%-CI_high`) - as.numeric(snp_data_selection$`99%-CI_low`)
        selected_row = which(confidence_interval_lengths == min(confidence_interval_lengths))
        selected_row2 = which(confidence_interval_lengths == max(confidence_interval_lengths))
        
        # Add row data to results
        results[id,1] = lumc_ids[id]
        results[id,2] = as.numeric(snp_data_selection$Result[selected_row])
        results[id,3] = as.numeric(snp_data_selection$`99%-CI_low`[selected_row])
        results[id,4] = as.numeric(snp_data_selection$`99%-CI_high`[selected_row])
        results[id,5] = as.numeric(snp_data_selection$Result[selected_row2])
        results[id,6] = as.numeric(snp_data_selection$`99%-CI_low`[selected_row2])
        results[id,7] = as.numeric(snp_data_selection$`99%-CI_high`[selected_row2])
      }
    }
  }
  
  # Return results
  return(results)
}

# Load and process SNP-based results
chr3p_snp_data = read.csv("digital-pcr/chr3p-snp.txt", sep="\t", stringsAsFactors = F, check.names=F)
chr3p_snp_processed_extra = process_snp_data_extra(ids, chr3p_snp_data)
chr3p_snp_processed_extra = chr3p_snp_processed_extra[which(!is.na(chr3p_snp_processed_extra[,1])),]

chr8q_snp_data = read.csv("digital-pcr/chr8q-snp.txt", sep="\t", stringsAsFactors = F, check.names=F)
chr8q_snp_processed_extra = process_snp_data_extra(ids, chr8q_snp_data)
chr8q_snp_processed_extra = chr8q_snp_processed_extra[which(!is.na(chr8q_snp_processed_extra[,1])),]

# Load data
chr3p_normalised = read.table("digital-pcr/chr3p-normalised.tsv", sep="\t")
chr8q_normalised = read.table("digital-pcr/chr8q-normalised.tsv", sep="\t")

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

x=as.numeric(c(chr3p_snp_processed_extra[,2], chr8q_snp_processed_extra[,2]))
y=as.numeric(c(chr3p_snp_processed_extra[,5], chr8q_snp_processed_extra[,5]))
cor.test(x[which(!is.na(x)&!is.na(y))],y[which(!is.na(x)&!is.na(y))], method = "spearman")
plot(x,y)
abline(0,1)
file = "digital-pcr/figure-multiple-snps.png"
png(file, res=600, 3200, 2300)  
par(mar=c(5,4,5,6))
xlim = c(0,13)
ylim = c(0,7.5)

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")
xlim = c(0,7.5)
mtext(text = "SNP assay 1", cex=1.1, col="#333333", side = 1, line = 2.75, at=3.5)
mtext(text = "SNP assay 2", cex=1.1, col="#333333", side = 2, line = 2.5, at=3.5)

xat = seq(xlim[1], xlim[2], by=1)
yat = seq(ylim[1], ylim[2], by=1)
axis(side = 1, at = xat, labels = xat, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
axis(side = 2, at = yat, labels = yat, las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")

segments(xat, ylim[1], xat, ylim[2], col="#eeeeee", lwd=1.4, xpd=T)
segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
segments(xlim[1], ylim[1], xlim[2], ylim[1], col="#b1b1b1", lwd=1.4, xpd=T)

for (i in 1:nrow(chr3p_snp_processed_extra)) {
  col = "#bfbfbf"
  if (chr3p_normalised[rownames(chr3p_snp_processed_extra)[i],5] != "normal") {
    col = "#93D095"
  }
  points(chr3p_snp_processed_extra[i,2], chr3p_snp_processed_extra[i,5], pch=16, cex=0.9, col=col)
}

for (i in 1:nrow(chr8q_snp_processed_extra)) {
  col = "#bfbfbf"
  if (chr8q_normalised[rownames(chr8q_snp_processed_extra)[i],5] != "normal") {
    col = "#F88E86"
  }
  if (stringr::str_detect(chr8q_normalised[rownames(chr8q_snp_processed_extra)[i],5],"gain")) {
    col = "#F0B27A"
  }
  points(chr8q_snp_processed_extra[i,2], chr8q_snp_processed_extra[i,5], pch=15, cex=0.8, col=col)
}


yy=0.525*2-0.35
text(8,5.6+yy,labels="Chromosome 3p", xpd=T, col="#333333",pos=4, cex=1.1)
points(8.6,4.9+yy, pch=16, cex=0.9, col="#bfbfbf", xpd=T)
text(8.8,4.9+yy,labels="normal", xpd=T, col="#333333",pos=4, cex=1.1)
points(8.6,4.2+yy, pch=16, cex=0.9, col="#93D095", xpd=T)
text(8.8,4.2+yy,labels="loss", xpd=T, col="#333333",pos=4, cex=1.1)

text(8,2.8+yy,labels="Chromosome 8q", xpd=T, col="#333333",pos=4, cex=1.1)
points(8.6,2.1+yy, pch=15, cex=0.8, col="#bfbfbf", xpd=T)
text(8.8,2.1+yy,labels="normal", xpd=T, col="#333333",pos=4, cex=1.1)
points(8.6,1.4+yy, pch=15, cex=.8, col="#F0B27A", xpd=T)
text(8.8,1.4+yy,labels="gain", xpd=T, col="#333333",pos=4, cex=1.1)
points(8.6,0.7+yy, pch=15, cex=.8, col="#F88E86", xpd=T)
text(8.8,0.7+yy,labels="amplification", xpd=T, col="#333333",pos=4, cex=1.1)

rect(0.2,5.2,3.8,6.8,col="white",border=NA)
text(0,6.35,labels="rho = 1.00", xpd=T, col="#333333",pos=4, cex=1.1)
text(0,5.65,labels="p < 0.001", xpd=T, col="#333333",pos=4, cex=1.1)

dev.off()
system(paste("open",file))


