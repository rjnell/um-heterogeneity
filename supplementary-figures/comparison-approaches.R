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

comparison = cbind(as.numeric(c(chr3p_cnv_processed[,2], chr8q_cnv_processed[,2])), 
         as.numeric(c(chr3p_snp_processed[,2], chr8q_snp_processed[,2])),
         c(rownames(chr3p_cnv_processed),rownames(chr8q_cnv_processed)))

data = rbind(
  #c("rs1062633", 1.16,1.08,1.24,19, "#F88E86", "#2196F3"),
  
  c("Classic approach", as.numeric(chr3p_cnv_processed[1,2:4]), 16, "#797979", "#2196F3"),
  c("SNP-based approach", as.numeric(chr3p_snp_processed[1,2:4]), 15, "#797979", "#2196F3"))

#)

title = "LUMC-01"
ymax = 2.1
ymin = 1.5
ystep = 0.1



file = paste0("supplementary-figures/comparison-approaches-1.png")
png(file,res=600,width=2700,height=2200)
par(mar=c(5,10,5,4))
xlim = c(0,2)
plot(xlim,c(ymin,ymax),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i")

xmax = 2
text(xmax/2,(ymax-ymin)/6+ymax, labels=title, col="#333333", xpd=T, cex=1.1, font=2)

# Plot y axis
yat = seq(ymin, ymax, by=ystep)
segments(0,yat,xmax,lwd=1.4,col="#EEEEEE", xpd=T)
axis(side = 2, at = yat, labels=format(round(yat*100)/100,nsmall=2),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
#mtext(side = 2, text = "Reference abundance", line = 5.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "Copy number value", line = 4.9, at = mean(yat), cex=1.1)
mtext(side = 2, text = "chromosome 3p", line = 4, at = mean(yat), cex=1.1)

segments(0,ymin,0,ymax,lwd=1.4,col="#b1b1b1",xpd=T)

#segments((1:(xmax/2))*2,ymin,(1:(xmax/2))*2,ymax,lwd=1.4,col="#eeeeee")
cols=NULL
k = 1
for (i in c(0.5,0.55,0.35)) {
  cols[k] = rgb(i,i,i)
  k = k +1
}
col="#808080"

for (i in 1:nrow(data)) {
  xm = 0.25
  xx= 0.1
  if(i == 3) {
    arrows(1.95, as.numeric(data[i,3]), 1.95, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col=data[i,6], lwd=1.4)
    #rect(i-0.7, as.numeric(data[i,2]), i-0.3, 0, pch=16, cex=0.9, col=col, border=col, lwd=1.4)
    points(1.95, as.numeric(data[i,2]), pch=as.numeric(data[i,5]), cex=0.8, col="#333333")
    text(1.95+.35, as.numeric(data[i,2]), labels=format(round(as.numeric(data[i,2])*100)/100,nsmall = 2), cex=0.95, col="#333333", xpd=T)
  }
  else if(i %%2 == 1) {
    arrows(i-0.5-xm+xx, as.numeric(data[i,3]), i-0.5-xm+xx, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col=data[i,6], lwd=1.4)
    #rect(i-0.7, as.numeric(data[i,2]), i-0.3, 0, pch=16, cex=0.9, col=col, border=col, lwd=1.4)
    points(i-0.5-xm+xx, as.numeric(data[i,2]), pch=as.numeric(data[i,5]), cex=0.9, col="#333333")
    text(i-0.5+xm-0.15+xx, as.numeric(data[i,2]), labels=format(round(as.numeric(data[i,2])*100)/100,nsmall = 2), cex=0.95, col="#333333", xpd=T)
  }
  else {
    arrows(i-0.5-xm-xx, as.numeric(data[i,3]), i-0.5-xm-xx, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col=data[i,6], lwd=1.4)
    #rect(i-0.7, as.numeric(data[i,2]), i-0.3, 0, pch=16, cex=0.9, col=col, border=col, lwd=1.4)
    points(i-0.5-xm-xx, as.numeric(data[i,2]), pch=as.numeric(data[i,5]), cex=0.8, col="#333333")
    text(i-0.5+xm-0.15-xx, as.numeric(data[i,2]), labels=format(round(as.numeric(data[i,2])*100)/100,nsmall = 2), cex=0.95, col="#333333", xpd=T)
    
  }
}


segments(xmax, ymin, xmax, ymax, lwd = 1.4, col= "#eeeeee", xpd=T)
segments(0, ymin, xmax, ymin, lwd = 1.4, col= "#b1b1b1", xpd=T)

i = 1-0.1
y2 = ymin- (ymax-ymin)/20 - (ymax-ymin)/7.5
y3 = y2 - (ymin - (ymax-ymin)/6 - (ymax-ymin)/7.5)


segments(i-0.55, ymin, i-0.55, y2, lwd = 1.4, col= "#b1b1b1", xpd=T)
segments(i-0.55, y2, i-0.55-0.2, y2, lwd = 1.4, col= "#b1b1b1", xpd=T) 
text(labels=data[i+0.1,1], x = i-0.55-0.2, y= y2, cex=1.1, col='#333333', pos=2, xpd=T)

i = 2-0.3
y2 = ymin - (ymax-ymin)/6 - (ymax-ymin)/7.5
segments(i-0.55, ymin, i-0.55, y2, lwd = 1.4, col= "#b1b1b1", xpd=T)
segments(i-0.55, y2, i-0.55-0.2, y2, lwd = 1.4, col= "#b1b1b1", xpd=T)
text(labels=data[i+0.3,1], x = i-0.5-0.2, y= y2, cex=1.1, col='#333333', pos=2, xpd=T)

x1 = 00
x2 = 2
ylim=c(ymin,ymax)
y2 = ylim[2]+ (ylim[2]-ylim[1])/5*1.25
y1 = ylim[2]+ (ylim[2]-ylim[1])/5*2.

segments(0, 2, nrow(data), lwd = 1.4, lty = 3, col = "#b1b1b1") 

dev.off()
system(paste("open", file))






data = rbind(
  #c("rs1062633", 1.16,1.08,1.24,19, "#F88E86", "#2196F3"),
  
  c("Classic approach", as.numeric(chr8q_cnv_processed[1,2:4]), 16, "#797979", "#2196F3"),
  c("SNP-based approach", as.numeric(chr8q_snp_processed[1,2:4]), 15, "#797979", "#2196F3"))

#)

title = "LUMC-01"
ymax = 2.5
ymin = 1.9
ystep = 0.1



file = paste0("supplementary-figures/comparison-approaches-2.png")
png(file,res=600,width=2700,height=2200)
par(mar=c(5,7,5,7))
xlim = c(0,2)
plot(xlim,c(ymin,ymax),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i")

text(nrow(data)/2,(ymax-ymin)/6+ymax, labels=title, col="#333333", xpd=T, cex=1.1, font=2)

# Plot y axis
yat = seq(ymin, ymax, by=ystep)
segments(0,yat,nrow(data),lwd=1.4,col="#EEEEEE", xpd=T)
axis(side = 2, at = yat, labels=format(round(yat*100)/100,nsmall=2),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
#mtext(side = 2, text = "Reference abundance", line = 5.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "Copy number value", line = 4.9, at = mean(yat), cex=1.1)
mtext(side = 2, text = "chromosome 8q", line = 4, at = mean(yat), cex=1.1)

segments(0,ymin,0,ymax,lwd=1.4,col="#b1b1b1",xpd=T)

segments((1:(nrow(data)/2))*2,ymin,(1:(nrow(data)/2))*2,ymax,lwd=1.4,col="#eeeeee", xpd=T)
segments(0, ymin, nrow(data), ymin, lwd = 1.4, col= "#b1b1b1", xpd=T)
k = 1
for (i in c(0.5,0.55,0.35)) {
  cols[k] = rgb(i,i,i)
  k = k +1
}
col="#808080"

for (i in 1:nrow(data)) {
  xm = 0.25
  xx= 0.1
  if(i %%2 == 1) {
    arrows(i-0.5-xm+xx, as.numeric(data[i,3]), i-0.5-xm+xx, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col=data[i,6], lwd=1.4)
    #rect(i-0.7, as.numeric(data[i,2]), i-0.3, 0, pch=16, cex=0.9, col=col, border=col, lwd=1.4)
    points(i-0.5-xm+xx, as.numeric(data[i,2]), pch=as.numeric(data[i,5]), cex=0.9, col="#333333")
    text(i-0.5+xm-0.15+xx, as.numeric(data[i,2]), labels=format(round(as.numeric(data[i,2])*100)/100,nsmall = 2), cex=0.95, col="#333333", xpd=T)
  }
  else {
    arrows(i-0.5-xm-xx, as.numeric(data[i,3]), i-0.5-xm-xx, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col=data[i,6], lwd=1.4)
    #rect(i-0.7, as.numeric(data[i,2]), i-0.3, 0, pch=16, cex=0.9, col=col, border=col, lwd=1.4)
    points(i-0.5-xm-xx, as.numeric(data[i,2]), pch=as.numeric(data[i,5]), cex=0.8, col="#333333")
    text(i-0.5+xm-0.15-xx, as.numeric(data[i,2]), labels=format(round(as.numeric(data[i,2])*100)/100,nsmall = 2), cex=0.95, col="#333333", xpd=T)
    
  }
}



segments(0, 2, nrow(data), lwd = 1.4, lty = 3, col = "#b1b1b1") 

i = 1-0.1
y2 = ymin- (ymax-ymin)/20 - (ymax-ymin)/7.5
segments(i-0.55, ymin, i-0.55, y2, lwd = 1.4, col= "#b1b1b1", xpd=T)
segments(i-0.55, y2, i-0.55-0.2, y2, lwd = 1.4, col= "#b1b1b1", xpd=T) 
text(labels=data[i+0.1,1], x = i-0.55-0.2, y= y2, cex=1.1, col='#333333', pos=2, xpd=T)

i = 2-0.3
y2 = ymin - (ymax-ymin)/6 - (ymax-ymin)/7.5
segments(i-0.55, ymin, i-0.55, y2, lwd = 1.4, col= "#b1b1b1", xpd=T)
segments(i-0.55, y2, i-0.55-0.2, y2, lwd = 1.4, col= "#b1b1b1", xpd=T)
text(labels=data[i+0.3,1], x = i-0.5-0.2, y= y2, cex=1.1, col='#333333', pos=2, xpd=T)

#segments(1-0.7, as.numeric(data[1,2]), 2.5, col=col, lwd=1.4, lty=3)
#segments(2-0.7, as.numeric(data[2,2]), 2.5, col=col, lwd=1.4, lty=3)
#arrows(2.5, as.numeric(data[2,2])+ymax/25, 2.5, as.numeric(data[1,2])-ymax/25, length=0.05, code=1, col="#333333", lwd=1.4)
segments()


#axis(side = 1, at = (0:(nrow(data)/2))*2, labels=rep("",nrow(data)/2+1), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
for (i in 1:nrow(data)) {
  #axis(side = 1, at = i-0.5, labels=data[i,1], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 1, tick = 0)
  #text(x = i-0.5, labels=data[i,1], srt=0, y = ymin-(ymax-ymin)/7, xpd=T, cex=1.1, col="red", srt=45, adj=1)
  #axis(side = 1, at = 1-0.5, labels=expression('var'[1]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
  #axis(side = 1, at = 2-0.5, labels=expression('var'[2]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
}

x1 = 00
x2 = 2
ylim=c(ymin,ymax)
y2 = ylim[2]+ (ylim[2]-ylim[1])/5*1.25
y1 = ylim[2]+ (ylim[2]-ylim[1])/5*2.

#segments(x1, y2, x1, mean(c(y1,y2)), col='#b1b1b1', lwd=1.4, xpd=T)
#segments(x2, y2, x2, mean(c(y1,y2)), col='#b1b1b1', lwd=1.4, xpd=T)
#segments(x1, mean(c(y1,y2)), x2, col='#b1b1b1', lwd=1.4, xpd=T)
#segments(mean(c(x1,x2)), mean(c(y1,y2)), mean(c(x1,x2)), y1, col='#b1b1b1', lwd=1.4, xpd=T)

dev.off()
system(paste("open", file))




