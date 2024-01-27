###
##  Figure 8 - Metastases- Landscape of transcriptional, molecular and clinical characteristics
###

# Close current images
system("taskkill /F /IM Microsoft.Photos.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=4, skip=0)

m = 82
bar = "darkblue"
s = 87
z = 10

file = "metastases/figure.png"
png(file, res=600, width=5000, height=7500)
par(mar=c(5,5,5,5))
plot(c(-15,95), 
     c(-10,40),
     type="n", 
     axes=F, 
     xaxs="i", 
     yaxs="i", 
     xlab="", 
     ylab="") 


white = "#EEEEEE"
border = NA
xi = 0.4


y = 26
x = 1
for (i in 1:10) {
  col = white
  
  if (i %% 2 == 1) {
    col = "#B5930B"  
  }
  if (i %% 2 == 0) {
    col = "#333333"  
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

y = 25

y = y-1
x = 1
for (i in 1:10) {
  col = white
  val = lumc_data$Gaq_mutation[i]
  if (!is.na(val)) {
    if (str_detect(val, "GNAQ p.Q209")) { 
      val = lumc_data$Gaq_mutation_clonality[i]
      if (!is.na(val)) {
        if (val == "clonal") { col = RColorBrewer::brewer.pal(9,"Paired")[6] }
        if (val == "subclonal") { col = RColorBrewer::brewer.pal(9,"Paired")[5] }
        if (val == "no assay") { col = "#999999" }
      }  
    }
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

y = y-1
x = 1
for (i in 1:10) {
  col = white
  val = lumc_data$Gaq_mutation[i]
  if (!is.na(val)) {
    if (str_detect(val, "GNA11 p.Q209")) { 
      val = lumc_data$Gaq_mutation_clonality[i]
      if (!is.na(val)) {
        if (val == "clonal") { col = RColorBrewer::brewer.pal(9,"Paired")[6] }
        if (val == "subclonal") { col = RColorBrewer::brewer.pal(9,"Paired")[5] }
        if (val == "no assay") { col = "#999999" }
      } 
    }
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

#y = y-1
x = 1
for (i in 1:10) {
  col = white
  #rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

#y = y-1
x = 1
for (i in 1:10) {
  col = white
  #rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}


y = y-2
x = 1
for (i in 1:10) {
  col = white
  #rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

#y = y-1
x = 1
for (i in 1:10) {
  col = white
  val = lumc_data$SF3B1_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { col = RColorBrewer::brewer.pal(9,"Paired")[6] }
    if (val == "subclonal") { col = RColorBrewer::brewer.pal(9,"Paired")[5] }
    if (val == "no assay") { col = "#999999" }
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

y = y-1
x = 1
for (i in 1:10) {
  col = white
  val = lumc_data$BAP1_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { col = RColorBrewer::brewer.pal(9,"Paired")[6] }
    if (val == "subclonal") { col = RColorBrewer::brewer.pal(9,"Paired")[5] }
    if (val == "no assay") { col = "#999999" }
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

y = y-2
x = 1
for (i in 1:10) {
  col = white
  val = lumc_data$chr3p_CNA_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { col = RColorBrewer::brewer.pal(9,"Paired")[6] }
    if (val == "subclonal") { col = RColorBrewer::brewer.pal(9,"Paired")[5] }
    if (val == "no assay") { col = "#999999" }
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  x = x + 2
  if (i %% 2 == 0) { x = x + 2 }
}

y = y-1
x = 1
for (i in 1:10) {
  col = white
  val = lumc_data$chr8q_CNA_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { col = RColorBrewer::brewer.pal(9,"Paired")[6] }
    if (val == "subclonal") { col = RColorBrewer::brewer.pal(9,"Paired")[5] }
    if (val == "no assay") { col = "#999999" }
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  
  if (i == 3) {
    rect(x-xi, y-0.5, x+xi+1, y-0.75, border=border, col="#2196F3", lwd=1.4)
  }
  if (i == 4) {
    rect(x-xi, y-0.5, x+xi+1, y-0.75, border=border, col="#4CAF50", lwd=1.4)
  }
  if (i == 7) {
    rect(x-xi, y-0.5, x+xi+1, y-0.75, border=border, col="#2196F3", lwd=1.4)
  }
  if (i == 8) {
    rect(x-xi, y-0.5, x+xi+1, y-0.75, border=border, col="#4CAF50", lwd=1.4)
  }
  
  x = x + 2
  if (i %% 2 == 0) { 
    text(x-2.5,y-1.5,cex=0.9,srt=90,adj=1,labels=lumc_data$LUMC_ID[i], col="#333333")
    x = x + 2 
  }
}




pos = 1-0.7-2




axis(side = 2, at = 26, labels = c("Lesion type"), las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)


axis(side = 2, at = 24:23, labels = c("GNAQ", "GNA11"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)



axis(side = 2, at = 21:20, labels = c("SF3B1","BAP1"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
axis(side = 2, at = 18:17, labels = c("Chr. 3p loss","Chr. 8q gain/amp."), las=2, cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)

#pos = c(pos, s)
segments(pos, 26+0.4, pos, 26-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 23+1.4, pos, 23-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 20+1.4, pos, 20-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 17+1.4, pos, 17-0.4, lwd=1.4, col="#b1b1b1")



y = 12.4+1
#axis(side = 1, at = c(s,s+2.5,s+5,s+7.5, s+10), labels= c("0%", NA, NA, NA, "100%") , xpd=T, pos = y, cex=1.1, col="#333333",  col.ticks = "#b1b1b1", lwd.ticks = 1.4)
#axis(side = 1, at = 30/2, labels= c("GEP class I") , xpd=T, pos = y, cex.axis=1.1, col="#333333",  col.ticks = "#b1b1b1", lwd.ticks = 1.4)
#axis(side = 1, at = 114/2, labels= c("GEP class II") , xpd=T, pos = y, cex.axis=1.1, col="#333333",  col.ticks = "#b1b1b1", lwd.ticks = 1.4)

#segments(1-0.4,y,4+0.4,y, col="#b1b1b1", lwd=1.4)
#segments(2.5,y,2.5,y-0.5, col="#b1b1b1", lwd=1.4)

xn = 40

yn = 26
y = yn
x =xn
text(x-2.5, y, labels = "Sample type", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#B5930B", xpd=T)
text(x+1, y, labels = "primary tumour", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#333333", xpd=T)
text(x+1, y, labels = "metastatic lesion", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


y =y -2
x =xn
text(x-2.5, y, labels = "Clonality", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#D4190C" , xpd=T)
text(x+1, y, labels = "clonal", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#F88E86" , xpd=T)
text(x+1, y, labels = "subclonal", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#999999" , xpd=T)
text(x+1, y, labels = "not analysed", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=white , xpd=T)
text(x+1, y, labels = "no alteration", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y =y -2

text(x-2.5, y, labels = "Chr. 8q", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#2196F3" , xpd=T)
text(x+1, y, labels = "allele I", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#4CAF50", xpd=T)
text(x+1, y, labels = "allele II", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

dev.off()
system(paste("open", file))





