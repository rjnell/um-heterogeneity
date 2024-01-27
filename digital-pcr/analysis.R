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
process_snp_data = function(ids, snp_data) {
  
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
      }
      
      # Add row data to results
      results[id,1] = lumc_ids[id]
      results[id,2] = as.numeric(snp_data_selection$Result[selected_row])
      results[id,3] = as.numeric(snp_data_selection$`99%-CI_low`[selected_row])
      results[id,4] = as.numeric(snp_data_selection$`99%-CI_high`[selected_row])
      results[id,5] = snp_data_selection$`Result interpretation`[selected_row]
      if (results[id,4] < 2) {
        results[id,6] = "loss"
      }
      else if (results[id,3] > 2) {
        results[id,6] = "gain/amplification"
      }
      else {
        results[id,6] = "not significant"
      }
      results[id,7] = paste0("#",snp_data_selection$`File ID`[selected_row])
      
    }
    else {
      results[id,1] = lumc_ids[id]
    }
  }
  
  # Return results
  return(results)
}

# Function to process raw data from Classic CNV experiments
# Stable reference is chosen from multiple/multiplex experiments
process_cnv_data = function(ids, cnv_data, target, references) {
  
  # Count ids
  n = length(ids)
  
  # Create matrix with final results
  results = matrix(data=NA, nrow=n, ncol=7)
  rownames(results) = ids
  
  # Iterate through ids to collect and process available data
  for (i in 1:n) {
    
    # Load id
    id = rownames(results)[i]
    
    # Load reference
    reference = references[id]
    
    # Select rows from cnv_data (id match + OK mark + target + reference)
    cnv_data_selection = cnv_data[which(stringr::str_detect(cnv_data$`File name`,id) & 
                                          stringr::str_detect(tolower(cnv_data$`Result interpretation`), tolower(target)) & 
                                          stringr::str_detect(tolower(cnv_data$`Result interpretation`), tolower(reference)) &
                                          cnv_data$Mark == "OK"), , drop=F]
    
    # If no rows available: no data, no results to add
    if (nrow(cnv_data_selection) > 0) {
      
      # If 1 row available: only data, to add
      if (nrow(cnv_data_selection) == 1) {
        selected_row = 1  
      }
      # If more rows available: add data with smallest CI
      else if (nrow(cnv_data_selection) > 1) {
        confidence_interval_lengths = as.numeric(cnv_data_selection$`99%-CI_high`) - as.numeric(cnv_data_selection$`99%-CI_low`)
        selected_row = which(confidence_interval_lengths == min(confidence_interval_lengths))
      }
      
      # Add row data to results
      results[id,1] = lumc_ids[id]
      results[id,2] = as.numeric(cnv_data_selection$Result[selected_row])
      results[id,3] = as.numeric(cnv_data_selection$`99%-CI_low`[selected_row])
      results[id,4] = as.numeric(cnv_data_selection$`99%-CI_high`[selected_row])
      results[id,5] = cnv_data_selection$`Result interpretation`[selected_row]
      if (results[id,4] < 2) {
        results[id,6] = "loss"
      }
      else if (results[id,3] > 2) {
        results[id,6] = "gain/amplification"
      }
      else {
        results[id,6] = "not significant"
      }
      results[id,7] = paste0("#",cnv_data_selection$`File ID`[selected_row])
      
    }
    else {
      results[id,1] = lumc_ids[id]
    }
  }
  
  # Return results
  return(results)
}

# Load and process SNP-based results
chr3p_snp_data = read.csv("digital-pcr/chr3p-snp.txt", sep="\t", stringsAsFactors = F, check.names=F)
chr3p_snp_processed = process_snp_data(ids, chr3p_snp_data)
chr8q_snp_data = read.csv("digital-pcr/chr8q-snp.txt", sep="\t", stringsAsFactors = F, check.names=F)
chr8q_snp_processed = process_snp_data(ids, chr8q_snp_data)

# Load and process Classic CNV results
cnv_data = read.csv("digital-pcr/chr3p-8q-cnv.txt", sep="\t", stringsAsFactors = F, check.names=F)
chr3p_cnv_processed = process_cnv_data(ids, cnv_data, "PPARG", references)
chr8q_cnv_processed = process_cnv_data(ids, cnv_data, "PTK2", references)

# Calculate overall correlation
cor.test(as.numeric(c(chr3p_cnv_processed[,2], chr8q_cnv_processed[,2])), 
         as.numeric(c(chr3p_snp_processed[,2], chr8q_snp_processed[,2])),
         method="spearman")

# Save processed CNV data
write.table(rbind(chr3p_cnv_processed, chr8q_cnv_processed), "digital-pcr/classic-cnv-processed.tsv", sep="\t")
write.table(rbind(chr3p_snp_processed, chr8q_snp_processed), "digital-pcr/snp-cnv-processed.tsv", sep="\t")

# Create figure for comparison Classic CNV and SNP-based results
{
  x=as.numeric(c(chr3p_cnv_processed[,2], chr8q_cnv_processed[,2]))
  y=as.numeric(c(chr3p_snp_processed[,2], chr8q_snp_processed[,2]))
  cor.test(x[which(!is.na(x)&!is.na(y))],y[which(!is.na(x)&!is.na(y))])
  plot(x,y)
  abline(0,1)
  file = "digital-pcr/figure-comparison-cnv.png"
  png(file, res=600, 2300, 2300)  
  par(mar=c(5,4,5,6))
  xlim = c(0,7)
  ylim = c(0,7)
  
  plot(xlim, 
       ylim, 
       type = "n", 
       axes = F, 
       xlab = "Classic CNV",
       ylab = "SNP-based CNV",
       xaxs = "i", 
       yaxs = "i")
  
  xat = seq(xlim[1], xlim[2], by=1)
  yat = seq(ylim[1], ylim[2], by=1)
  axis(side = 1, at = xat, labels = xat, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
  axis(side = 2, at = yat, labels = yat, las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
  
  segments(xat, ylim[1], xat, ylim[2], col="#eeeeee", lwd=1.4, xpd=T)
  segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
  segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
  segments(xlim[1], ylim[1], xlim[2], ylim[1], col="#b1b1b1", lwd=1.4, xpd=T)
  
  for (i in 1:nrow(chr3p_cnv_processed)) {
    
    if (!is.na(chr3p_cnv_processed[i,2]) & !is.na(chr3p_snp_processed[i,2])) {
      col = "#b1b1b1"
      points(chr3p_cnv_processed[i,2], chr3p_snp_processed[i,2], pch=16, cex=0.8, col=col)
    }
  }
  
  for (i in 1:nrow(chr8q_cnv_processed)) {
    
    if (!is.na(chr8q_cnv_processed[i,2]) & !is.na(chr8q_snp_processed[i,2])) {
      col = "#b1b1b1"
      if (!is.na(lumc_data$chr8q_CNA_multiallelic[which(lumc_data$Tumor_ID == ids[i])])) {
        col="red"
      }
      points(chr8q_cnv_processed[i,2], chr8q_snp_processed[i,2], pch=15, cex=0.8, col=col)
    }
  }

  segments(1,1,6.3,6.3,lwd=1.4,lty=3, col="#b1b1b1")
  
  points(7.5,3.25, pch=15, cex=1, col="#b1b1b1", xpd=T)
  text(7.5,3.25,labels="8q", xpd=T, col.axis="#333333",pos=4)
  points(7.5,3.75, pch=16, cex=1, col="#b1b1b1", xpd=T)
  text(7.5,3.75,labels="3p", xpd=T, col.axis="#333333",pos=4)
  
  dev.off()
  system(paste("open",file))
}

# Select best chromosome 3p measurement (SNP-based > Classic CNV)
chr3p = NULL
for (id in ids) {
  if (!is.na(chr3p_snp_processed[id,2])) {
    SNP = ""
    SNP = strsplit(chr3p_snp_processed[id,5], " ")[[1]][2]
    SNP = substr(SNP,5,nchar(SNP))
    chr3p = rbind(chr3p, c(chr3p_snp_processed[id,], paste0("SNP-based CNV (",SNP,")")))
  }
  else {
    CNV = ""
    CNV = chr3p_cnv_processed[id,5]
    CNV = substr(CNV,13,nchar(CNV))
    CNV = stringr::str_replace(CNV, "\\?", "/")
    chr3p = rbind(chr3p, c(chr3p_cnv_processed[id,], paste0("Classic CNV (",CNV,")")))
  }
}
rownames(chr3p) = ids
write.table(chr3p, "digital-pcr/chr3p-unnormalised.tsv", sep="\t")

# Select best chromosome 8q measurement (SNP-based > Classic CNV, except for multi-allelic amplifications)
chr8q = NULL
for (id in ids) {
  if (!is.na(lumc_data$chr8q_CNA_multiallelic[which(lumc_data$Tumor_ID == id)])) {
    CNV = ""
    CNV = chr8q_cnv_processed[id,5]
    CNV = substr(CNV,13,nchar(CNV))
    CNV = stringr::str_replace(CNV, "\\?", "/")
    chr8q = rbind(chr8q, c(chr8q_cnv_processed[id,], paste0("Classic CNV (",CNV,") [multiallelic]"))) 
  }
  else if (!is.na(chr8q_snp_processed[id,2])) {
    SNP = ""
    SNP = strsplit(chr8q_snp_processed[id,5], " ")[[1]][1]
    SNP = substr(SNP,14,nchar(SNP))
    chr8q = rbind(chr8q, c(chr8q_snp_processed[id,], paste0("SNP-based CNV (",SNP,")")))
  }
  else {
    CNV = ""
    CNV = chr8q_cnv_processed[id,5]
    CNV = substr(CNV,13,nchar(CNV))
    CNV = stringr::str_replace(CNV, "\\?", "/")
    chr8q = rbind(chr8q, c(chr8q_cnv_processed[id,], paste0("Classic CNV (",CNV,")"))) 
  }
}
rownames(chr8q) = ids
write.table(chr8q, "digital-pcr/chr8q-unnormalised.tsv", sep="\t")

# Load purity data
purity = read_xlsx("data/Supplementary Tables.xlsx", sheet=2, skip=0)

# Function to give an interpretation of the corrected CNA
check_cna_corrected = function(cna_corrected) {
  
  cna_interpretation = "clonal loss (-1)"
  
  if (cna_corrected[2] < 0 & cna_corrected[3] > 0) {
    cna_interpretation = "normal"
  }
  
  if (cna_corrected[2] < -1 & cna_corrected[3] > -1) {
    cna_interpretation = "clonal loss (-1)"
  }
  else if (cna_corrected[2] > -1 & cna_corrected[3] < 0) {
    cna_interpretation = "subclonal loss (-1)"
  }
  
  if (cna_corrected[2] > 0 & cna_corrected[3] < 1) {
    cna_interpretation = "subclonal gain (+1)"
  }
  else if (cna_corrected[2] < 1 & cna_corrected[3] > 1) {
    cna_interpretation = "clonal gain (+1)"
  }
  
  for (v in 2:7) {
    if (cna_corrected[2] > v-1 & cna_corrected[3] < v) {
      cna_interpretation = paste0("subclonal amplification (+",v,")")
    }
    else if (cna_corrected[2] < v & cna_corrected[3] > v) {
      cna_interpretation = paste0("clonal amplification (+",v,")")
    }
  }
  
  return(cna_interpretation)
}

# Normalise chromosome 3p measurements for purity
chr3p_normalised = matrix(nrow = 80, ncol = 5)
for (i in 1:80) {
  id = lumc_data$Tumor_ID[i]
  chr3p_tumour = as.numeric(chr3p[id,2:4])-2
  purity_tumour = as.numeric(purity[which(purity$Tumor_ID == id)[1],5:7])/100
  cna_corrected = calc_ratio(chr3p_tumour, purity_tumour)
  cna_interpretation = check_cna_corrected(cna_corrected)
  chr3p_normalised[i,] = c(chr3p[id,8], cna_corrected, cna_interpretation)
}
rownames(chr3p_normalised) = lumc_data$Tumor_ID
write.table(chr3p_normalised, "digital-pcr/chr3p-normalised.tsv", sep="\t")

# Normalise chromosome 8q measurements for purity
chr8q_normalised = matrix(nrow = 80, ncol = 5)
for (i in 1:80) {
  id = lumc_data$Tumor_ID[i]
  chr8q_tumour = as.numeric(chr8q[id,2:4])-2
  purity_tumour = as.numeric(purity[which(purity$Tumor_ID == id)[1],5:7])/100
  cna_corrected = calc_ratio(chr8q_tumour, purity_tumour)
  cna_interpretation = check_cna_corrected(cna_corrected)
  chr8q_normalised[i,] = c(chr8q[id,8], cna_corrected, cna_interpretation)
}
rownames(chr8q_normalised) = lumc_data$Tumor_ID
write.table(chr8q_normalised, "digital-pcr/chr8q-normalised.tsv", sep="\t")
