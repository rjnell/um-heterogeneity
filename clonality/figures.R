# Scripts to visualise individual comparisons

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)
lumc_data = lumc_data[order(lumc_data$LUMC_ID),]
gaq = read_xlsx("data/Supplementary Tables.xlsx", sheet=2)
bse = read_xlsx("data/Supplementary Tables.xlsx", sheet=3)
chr3p_unnormalised = read.table("digital-pcr/chr3p-unnormalised.tsv", sep="\t", stringsAsFactors = F, row.names = 1)
chr8q_unnormalised = read.table("digital-pcr/chr8q-unnormalised.tsv", sep="\t", stringsAsFactors = F, row.names = 1)
chr8q_normalised = read.table("digital-pcr/chr8q-normalised.tsv", sep="\t", stringsAsFactors = F, row.names = 1)

col_subclonal = "#F88E86"
col_clonal = "#D4190C"

# Iterate through tumours
for (id in "###") {
  
  file = paste0("clonality/figures/",id,".png")
  png(file, res=600, width=3500, height=2200)
  par(mar=c(7,7,3,7))
  xlim = c(0,8)
  max = 100
  plot(xlim,c(-50,max),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i",main="")

  # Load data
  id_data = NULL
  
  # Gaq
  id_gaq_w = which(gaq$Tumor_ID == id)
  if (length(id_gaq_w) == 1) {
    id_data = rbind(id_data, c(as.character(gaq[id_gaq_w, 3]), as.numeric(gaq[id_gaq_w, 5:7]), col_clonal))
  }
  if (id == "###") {
    id_gaq_w = which(gaq$Tumor_ID == "### (2nd)")
    id_data = rbind(id_data, c(as.character(gaq[id_gaq_w, 3]), as.numeric(gaq[id_gaq_w, 5:7]), col_clonal))
  }
  
  # BSE = EIF1AX
  id_bse_w = which(bse$Tumor_ID == id & bse$Gene == "EIF1AX")
  if (length(id_bse_w) == 1) {
    r = digitalPCRsimulations::calc_ratio(as.numeric(bse[id_bse_w, 5:7]), as.numeric(id_data[1,2:4]))
    if (r[3] < 1) {
      col = col_subclonal
    }
    else {
      col = col_clonal
    }
    id_data = rbind(id_data, c(paste(as.character(bse[id_bse_w, 2]),as.character(bse[id_bse_w, 3])), as.numeric(bse[id_bse_w, 5:7]), col))
  }
  
  # BSE = SF3B1
  id_bse_w = which(bse$Tumor_ID == id & bse$Gene == "SF3B1")
  if (length(id_bse_w) == 1) {
    r = digitalPCRsimulations::calc_ratio(as.numeric(bse[id_bse_w, 5:7]), as.numeric(id_data[1,2:4]))
    if (r[3] < 1) {
      col = col_subclonal
    }
    else {
      col = col_clonal
    }
    id_data = rbind(id_data, c(paste(as.character(bse[id_bse_w, 2]),as.character(bse[id_bse_w, 3])), as.numeric(bse[id_bse_w, 5:7]), col))
  }
  
  # 3p
  id_3p_w = which(rownames(chr3p_unnormalised) == id)
  if (length(id_3p_w) == 1) {
    if (chr3p_unnormalised[id_3p_w,4] < 2) {
      r = digitalPCRsimulations::calc_ratio(as.numeric(chr3p_unnormalised[id_3p_w,2:4]-2)*-100, as.numeric(id_data[1,2:4]))
      if (r[3] < 1) {
        col = col_subclonal
      }
      else {
        col = col_clonal
      }
      id_data = rbind(id_data, c("Chr.3p loss", as.numeric(chr3p_unnormalised[id_3p_w,2:4]-2)*-100, col))
    }
  }
  
  # BSE = BAP1
  id_bse_w = which(bse$Tumor_ID == id & bse$Gene == "BAP1")
  if (length(id_bse_w) == 1) {
    b = digitalPCRsimulations::calc_ratio(as.numeric(chr3p_unnormalised[id_3p_w,2:4]), 1/as.numeric(bse[id_bse_w, 5:7]/100)[c(1,3,2)])
    r = digitalPCRsimulations::calc_ratio(b*100, as.numeric(id_data[1,2:4]))
    if (r[3] < 1) {
      col = col_subclonal
    }
    else {
      col = col_clonal
    }
    id_data = rbind(id_data, c(paste(as.character(bse[id_bse_w, 2]),as.character(bse[id_bse_w, 3])), b*100, col))
  }
  else {
    if (!is.na(lumc_data$BAP1_mutation_RNA[which(lumc_data$Tumor_ID == id)])) {
      vv = lumc_data$BAP1_mutation_RNA[which(lumc_data$Tumor_ID == id)]
      id_data = rbind(id_data, c("BAP1 mutation", 0, 0, 0, col))
    }
    
  }
  
  # 8q
  id_8q_w = which(rownames(chr8q_unnormalised) == id)
  if (length(id_8q_w) == 1) {
    if (chr8q_normalised[id_8q_w,3] < 0) {
    }
    else if (chr8q_normalised[id_8q_w,3] < 1) {
      r = digitalPCRsimulations::calc_ratio((as.numeric(chr8q_unnormalised[id_8q_w,2:4])-2)*100, as.numeric(id_data[1,2:4]))
      if (r[3] < 1) {
        col = col_subclonal
      }
      else {
        col = col_clonal
      }
      id_data = rbind(id_data, c("Chr.8q gain", (as.numeric(chr8q_unnormalised[id_8q_w,2:4])-2)*100, col))
    }
    else if (chr8q_normalised[id_8q_w,3] > 1) {
      id_data = rbind(id_data, c("Chr.8q amplification", 0, 0, 0, col))
    }
  }
  
  # Gaq imbalance
  id_bse_w = which(bse$Tumor_ID == id & bse$Gene %in% c("GNA11","GNAQ"))
  if (length(id_bse_w) == 1) {
    v = as.numeric(bse[id_bse_w, 5:7])
    if (v[1] > 2) {
      v = (v-2)*100
    } else  {
      v = (2-v)*100
      v = v[c(1,3,2)]
    }
    r = digitalPCRsimulations::calc_ratio(v, as.numeric(id_data[1,2:4]))
    if (r[3] < 1) {
      col = col_subclonal
    }
    else {
      col = col_clonal
    }
    id_data = rbind(id_data, c(paste(as.character(bse[id_bse_w, 2]),as.character(bse[id_bse_w, 3])), v, col))
  }
  
  
  if (id_data[1,1] == "n/a") {
    id_data = id_data[2:nrow(id_data),]
  } 
  
  
  xmax = nrow(id_data) + 0.5

  # Plot y axis
  yat = seq(0, max, by=25)
  segments(0,yat,xmax,lwd=1.4,col="#EEEEEE", xpd=T)
  axis(side = 2, at = yat, labels=paste0(yat,"%"),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
  mtext(side = 2, text = "Fraction cells", line = 4, at = mean(yat), cex=1.1, col.axis="#333333", col="#333333")
  axis(side = 1, at = (1:xmax)-0.25, labels=rep("",xmax), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, pos = 0)

  for (i in 1:nrow(id_data)) {
      rect(i-0.35-0.25, 0, i+0.35-0.25, as.numeric(id_data[i,2]), border=NA, col= id_data[i,5]) 
      arrows(i-0.25, as.numeric(id_data[i,3]), i-0.25, as.numeric(id_data[i,4]), length = 0.05, angle = 90, code = 3, col ="#b1b1b1", lwd = 1.4, xpd = T)
      mut = strsplit(id_data[i,1], " ")[[1]]
      font = 3
      if (mut[1] == "Chr.8q") {
        mut[1] = "Chr. 8q"
        font = 1
      }
      if (mut[1] == "Chr.3p") {
        mut[1] = "Chr. 3p"
        font = 1
      }
  
      if (mut[2] == "amplification") {
        text(i-0.25, 10, xpd=T, labels="n/c",col="#333333", cex=0.9)
      }
      
      if (mut[1] == "BAP1" & as.numeric(id_data[i,4]) == 0) {
        text(i-0.25, 10, xpd=T, labels="n/a",col="#333333", cex=0.9)
      }
      
      text(i-0.7+0.25, -17, xpd=T, labels=mut[1],srt=90,adj=1, col="#333333", cex=1.1, font=font)
      text(i-0.3+0.25, -17, xpd=T, labels=mut[2],srt=90,adj=1, col="#333333", cex=0.9)
  }
  
  
  text(xmax/2, 120, xpd=T, labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID == id)], col="#333333", cex=1.1, font=2)
  
  segments(0,0,xmax,lwd=1.4,col="#b1b1b1", xpd=T)
  dev.off()
  #system(paste("open", file))
  
}
