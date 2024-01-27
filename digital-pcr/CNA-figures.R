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

# Figure template
figure_8q = function(file, ym, ids) {
  
ym= ym

data = NULL#"08-017", "09-028", 
e_ids = ids #c("08-017", "09-028")
#e_ids = c("08-028", "08-012")#, "08-003")
###e_ids = c("10-006", "05-020")#, "08-038", "09-033")
#e_ids = c("09-033", "08-003")
for (id in e_ids) {
  
  cnv = as.numeric(chr8q_cnv_processed[id,2:4])
  var1 = 1/(as.numeric(chr8q_snp_processed[id,2:4]))
  var2 = 1 - 1/(as.numeric(chr8q_snp_processed[id,2:4]))
  
  raw = chr8q_snp_data[which(chr8q_snp_data$`File ID` == as.numeric(substr(chr8q_snp_processed[id,7],2,5))), ]
  if (raw$`Cchannel 1` > raw$`Cchannel 2`) {
    var2 = 1/(as.numeric(chr8q_snp_processed[id,2:4]))
    var1 = 1 - 1/(as.numeric(chr8q_snp_processed[id,2:4]))
  }
  
  data = rbind(data, c(id, calc_ratio(cnv, 1/var1), calc_ratio(cnv, 1/var2)))
}

#ym= 5

file = file
png(file, res=600, 3200, 2300)  
par(mar=c(5,6,5,4))
xlim = c(0,5)
ylim = c(0,7.5)

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

ymax=ym+.5
# Plot y axis
yat = seq(00, ymax, by=1)
segments(0,yat,nrow(data)*1.5,lwd=1.4,col="#EEEEEE", xpd=T)
axis(side = 2, at = yat, labels=yat,las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
#mtext(side = 2, text = "Reference abundance", line = 5.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "Copy number", line = 3.4, cex=1.1, col="#333333",at=ym/2)
mtext(side = 2, text = "chromosome 8q", line = 2.5, cex=1.1, col="#333333",at=ym/2)

segments(0,0,0,ymax,lwd=1.4,col="#b1b1b1",xpd=T)

segments((1:nrow(data))*1.5,0,(1:nrow(data))*1.5,ymax,lwd=1.4,col="#eeeeee", xpd=T)

cols = c("#333333","#333333","#333333","#b1b1b1","#b1b1b1", "#333333")#$ RColorBrewer::brewer.pal(9,"Paired")[c(2,4,6)]

k = 1
for (i in c(0.85,0.75,0.35)) {
  cols[k] = rgb(i,i,i)
  k = k +1
}
#col="#808080"

for (i in 1:nrow(data)) {
  j = i*1.5
  rect(j-1.2, as.numeric(data[i,5]), j-0.8, 0, pch=16, cex=0.9, col=cols[1], border=NA, lwd=1.4, xpd=T)
  arrows(j-1, as.numeric(data[i,6]), j-1, as.numeric(data[i,7]), length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  rect(j-0.7, as.numeric(data[i,2]), j-0.3, 0, pch=16, cex=0.9, col=cols[2], border=NA, lwd=1.4, xpd=T)
  arrows(j-0.5, as.numeric(data[i,3]), j-0.5, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  
  #text(j-0.5, as.numeric(data[i,4])+ymax/10, labels=paste0(format(as.numeric(data[i,2]), nsmall = 0), ""), cex=0.8, col="#333333", xpd=T)
  
}

#segments(1-0.7, as.numeric(data[1,2]), 2.5, col=col, lwd=1.4, lty=3)
#segments(2-0.7, as.numeric(data[2,2]), 2.5, col=col, lwd=1.4, lty=3)
#arrows(2.5, as.numeric(data[2,2])+ymax/25, 2.5, as.numeric(data[1,2])-ymax/25, length=0.05, code=1, col="#333333", lwd=1.4)


axis(side = 1, at = (0:nrow(data))*1.5, labels=rep("",nrow(data)+1), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, xpd=T, tick = 0)
for (i in 1:nrow(data)) {
  #axis(side = 3, at = i*1.5-0.75, font=2,labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID==data[i,1])], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0, xpd=T)
  text(i*1.5-0.75, ym+1.5,font=2,labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID==data[i,1])], col="#333333", cex=1.1, xpd=T)
  #axis(side = 1, at = 1-0.5, labels=expression('var'[1]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
  #axis(side = 1, at = 2-0.5, labels=expression('var'[2]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
}

segments(0,0,nrow(data)*1.5,lwd=1.4,col="#b1b1b1", xpd=T)

dev.off()
system(paste("open", file))
}

figure_8q("digital-pcr/CNAs/fig-2.png", 2, c("08-017", "09-028"))
figure_8q("digital-pcr/CNAs/fig-3.png", 5, c("08-028", "08-012"))


05-020
08-003
08-038
09-033
10-006
10-043

figure_8q("digital-pcr/CNAs/fig-4.png", 5, c("09-033", "08-003"))

figure_8q("digital-pcr/CNAs/fig-4-supp.png", 3, c("08-038", "10-043", "10-006", "05-020"))


ym= ym

data = NULL#"08-017", "09-028", 
e_ids = c("08-017", "09-004")

for (id in e_ids) {
  
  cnv = as.numeric(chr3p_cnv_processed[id,2:4])
  var1 = 1/(as.numeric(chr3p_snp_processed[id,2:4]))
  var2 = 1 - 1/(as.numeric(chr3p_snp_processed[id,2:4]))
  
  raw = chr3p_snp_data[which(chr3p_snp_data$`File ID` == as.numeric(substr(chr3p_snp_processed[id,7],2,5))), ]
  if (raw$`Cchannel 1` > raw$`Cchannel 2`) {
    var2 = 1/(as.numeric(chr3p_snp_processed[id,2:4]))
    var1 = 1 - 1/(as.numeric(chr3p_snp_processed[id,2:4]))
  }
  
  data = rbind(data, c(id, calc_ratio(cnv, 1/var1), calc_ratio(cnv, 1/var2)))
}

ym= 2

file = "digital-pcr/CNAs/fig-1.png"
png(file, res=600, 3200, 2300)  
par(mar=c(5,6,5,4))
xlim = c(0,5)
ylim = c(0,7.5)

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

ymax=ym+.5
# Plot y axis
yat = seq(00, ymax, by=1)
segments(0,yat,nrow(data)*1.5,lwd=1.4,col="#EEEEEE", xpd=T)
axis(side = 2, at = yat, labels=yat,las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
#mtext(side = 2, text = "Reference abundance", line = 5.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "Copy number", line = 3.4, cex=1.1, col="#333333",at=ym/2)
mtext(side = 2, text = "chromosome 3p", line = 2.5, cex=1.1, col="#333333",at=ym/2)

segments(0,0,0,ymax,lwd=1.4,col="#b1b1b1",xpd=T)

segments((1:nrow(data))*1.5,0,(1:nrow(data))*1.5,ymax,lwd=1.4,col="#eeeeee", xpd=T)

cols = c("#333333","#333333","#333333","#b1b1b1","#b1b1b1", "#333333")#$ RColorBrewer::brewer.pal(9,"Paired")[c(2,4,6)]

k = 1
for (i in c(0.85,0.75,0.35)) {
  cols[k] = rgb(i,i,i)
  k = k +1
}
#col="#808080"

for (i in 1:nrow(data)) {
  j = i*1.5
  rect(j-1.2, as.numeric(data[i,5]), j-0.8, 0, pch=16, cex=0.9, col=cols[1], border=NA, lwd=1.4, xpd=T)
  arrows(j-1, as.numeric(data[i,6]), j-1, as.numeric(data[i,7]), length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  rect(j-0.7, as.numeric(data[i,2]), j-0.3, 0, pch=16, cex=0.9, col=cols[2], border=NA, lwd=1.4, xpd=T)
  arrows(j-0.5, as.numeric(data[i,3]), j-0.5, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  
  #text(j-0.5, as.numeric(data[i,4])+ymax/10, labels=paste0(format(as.numeric(data[i,2]), nsmall = 0), ""), cex=0.8, col="#333333", xpd=T)
  
}

#segments(1-0.7, as.numeric(data[1,2]), 2.5, col=col, lwd=1.4, lty=3)
#segments(2-0.7, as.numeric(data[2,2]), 2.5, col=col, lwd=1.4, lty=3)
#arrows(2.5, as.numeric(data[2,2])+ymax/25, 2.5, as.numeric(data[1,2])-ymax/25, length=0.05, code=1, col="#333333", lwd=1.4)


axis(side = 1, at = (0:nrow(data))*1.5, labels=rep("",nrow(data)+1), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, xpd=T, tick = 0)
for (i in 1:nrow(data)) {
  #axis(side = 3, at = i*1.5-0.75, font=2,labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID==data[i,1])], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0, xpd=T)
  text(i*1.5-0.75, ym+1.5,font=2,labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID==data[i,1])], col="#333333", cex=1.1, xpd=T)
  #axis(side = 1, at = 1-0.5, labels=expression('var'[1]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
  #axis(side = 1, at = 2-0.5, labels=expression('var'[2]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
}

segments(0,0,nrow(data)*1.5,lwd=1.4,col="#b1b1b1", xpd=T)

dev.off()
system(paste("open", file))
