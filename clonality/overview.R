###
##  Figure 4 - Clonality analysis
###

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)
lumc_data = lumc_data[order(lumc_data$order),]

purity = read_xlsx("data/Supplementary Tables.xlsx", sheet=2)

colors = list()
colors[["female"]] = RColorBrewer::brewer.pal(12, "Paired")[5]
colors[["Class I"]] = RColorBrewer::brewer.pal(10, "Paired")[9]#777777"
colors[["Class II"]] = RColorBrewer::brewer.pal(10, "Paired")[10]#"#b1b1b1"
colors[["PLCB4"]] = RColorBrewer::brewer.pal(8,"YlGnBu")[4]
colors[["CYSLTR2"]] = RColorBrewer::brewer.pal(8,"YlGnBu")[6]
colors[["EIF1AX-exon1"]] = RColorBrewer::brewer.pal(10, "Paired")[10]
colors[["EIF1AX-exon2"]] = RColorBrewer::brewer.pal(9, "Paired")[9]
colors[["SF3B1-hotspot"]] = RColorBrewer::brewer.pal(10, "Paired")[8]
colors[["SF3B1-other"]] = RColorBrewer::brewer.pal(10, "Paired")[7]
colors[["BAP1"]] = RColorBrewer::brewer.pal(12, "Paired")[12]

colors[["loss"]] = RColorBrewer::brewer.pal(12, "Paired")[6]
colors[["gain"]] = RColorBrewer::brewer.pal(11, "RdYlGn")[9]
colors[["amp"]] = RColorBrewer::brewer.pal(11, "RdYlGn")[10]

greens = colorRampPalette(RColorBrewer::brewer.pal(11, "PRGn")[6:11])(12)
colors[["subclonal gain"]] = RColorBrewer::brewer.pal(11, "Paired")[7]
colors[["clonal gain"]] = RColorBrewer::brewer.pal(11, "Paired")[8]
colors[["subclonal amplification"]] = "#F88E86"
colors[["clonal amplification"]] = "#D4190C"

colors[["clonal loss"]] = RColorBrewer::brewer.pal(12, "Paired")[4]
colors[["subclonal loss"]] = RColorBrewer::brewer.pal(12, "Paired")[3]


colors[["mut-subclonal"]] = "#F88E86"#0B72C4"
colors[["mut-clonal"]] = "#D4190C"#7AC0F8"
colors[["mut-unknown"]] = "#B1b1b1"


m = 82
bar = "darkblue"
s = 87
z = 10

file = "clonality/overview.png"
png(file, res=600, width=7250, height=7000)
par(mar=c(5,7,5,7))
plot(c(-40,90), 
     c(-10,42),
     type="n", 
     axes=F, 
     xaxs="i", 
     yaxs="i", 
     xlab="", 
     ylab="") 

#pos=c(s,s+2.5,s+5,s+7.5, s+10)
#segments(pos, 29+0.4, pos, 29-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 14+1.4, pos, 14-0.4, lwd=1.4, col="#eeeeee", xpd=T)

y = 40.4

white = "#EEEEEE"
border = NA
xi = 0.4


y = 27

#y = y-
#for (i in 1:82) {
#  col = white
#  x = i 
#  if (i > 30) { x = x + 2 }
#  val = lumc_data$survival_information[i]
#  if (val == "death by melanoma") { n1=n1+1;col = "#777777"}
#  if (val == "death by unknown") { n2=n2+1;col = "#B1b1b1"}
#  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
#}



get_x = function (i) {
  x = i 
  if (i > 8) { x = x + 2 }
  if (i > 22) { x = x + 2 }
  if (i > 69) { x = x + 2 }
  return(x)
}


# Gaq boxes
y = 29

y = y-1
for (i in 1:80) {
  col = white
  x = get_x(i)
  val = lumc_data$survival_information[i]
  if (!is.na(val)) {
    if (val == "death by melanoma") { col = colors[["Class I"]]}
    else if (val == "2") { col = colors[["Class II"]]}
  } 
  #rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col, lwd=1.4)
}

y = y-1



  y = y-4+5
  for (i in 1:80) {
    col = white
    x = get_x(i)
    val = lumc_data$Gaq_mutation_clonality[i]
    if (!is.na(val)) {
      if (val == "clonal") { col = colors[["mut-clonal"]]}
      else if (val == "subclonal") { col = colors[["mut-subclonal"]]}
      else { col = colors[["mut-unknown"]]}
    } 
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col, lwd=1.4)
  }



# Gaq axis (24.6-28.4)
axis(side = 2, at = 28, labels = c(expression(italic("Gaq")~"mutation")), las=2, cex=1.1, cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, 28.4, pos, 27.6, lwd=1.4, col="#b1b1b1")

#axis(side = 2, at = 23.5, labels = c("Gaq CNA"), las=2, cex=1.1, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
#segments(pos, 22.9, pos, 24.1, lwd=1.4, col="#b1b1b1")

# Gaq CNA (23.4-19.6)
#y = 21.5
#segments(pos, y+c(-1.5,0,1.5), pos+89, lwd=1.4, col="#eeeeee")
#axis(side = 2, at = y+c(-1.5,0,1.5), labels = c("-1",0,"+1"), las=2, cex=1.1, cex.axis=1.1,col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
#segments(pos, y-1.9, pos, y+1.9, lwd=1.4, col="#b1b1b1")

y = 17+5
axis(side = 2, at = (y):y, labels = c(
                                        expression(italic("BAP1")~"mutation")), las=2, cex=1.1,cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, y+.4, pos, y-1+.6, lwd=1.4, col="#b1b1b1")


y = 16+5
axis(side = 2, at = (y-2):y, labels = c(expression(italic("Gaq")~"imbalance"),"Chromosome 8q amp.",
                                        "Chromosome 8q gain"), las=2, cex=1.1,cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, y+1.4, pos, y-3+.6, lwd=1.4, col="#b1b1b1")


# BSE (18.4)
y = 21+5
axis(side = 2, at = (y-2):y, labels = c("Chromosome 3p loss",expression(italic("SF3B1")~"mutation"),expression(italic("EIF1AX")~"mutation")), las=2, cex=1.1,cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, y+.4, pos, y-3+.6, lwd=1.4, col="#b1b1b1")
for (i in 1:80) {
  col = white
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$EIF1AX_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { 
      col = colors[["mut-clonal"]]
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
    else if (val == "subclonal") { 
      col = colors[["mut-subclonal"]]
      rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
    }
    else { 
      col = "#b1b1b1"
      points(x,y,pch=8,cex=0.5,col=col)
    }
  } 
}

y = y-1
for (i in 1:80) {
  col = white
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$SF3B1_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { 
      col = colors[["mut-clonal"]]
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
    else if (val == "subclonal") { 
      col = colors[["mut-subclonal"]]
      rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
    }
    else { 
      #col = "#333333"
      #points(x,y,pch=8,cex=0.5,col=col)
      col = "#b1b1b1"
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
  } 
}

y = y-1

for (i in 1:80) {
  col = white
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$chr3p_CNA_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { col = colors[["mut-clonal"]] 
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)}
    
    else if (val == "subclonal") { col = colors[["mut-subclonal"]]
    rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)}
    else { col = colors[["mut-unknown"]]}
  } 
  #rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}

















y=y-2

for (i in 1:80) {
  col = white
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$BAP1_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { 
      col = colors[["mut-clonal"]]
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
    else if (val == "subclonal") { 
      col = colors[["mut-subclonal"]]
      rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
    }
    else { 
      col = "#b1b1b1"
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
  } 
}

y = y-1

for (i in 1:80) {
  col = white
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$chr8q_CNA_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { col = colors[["mut-clonal"]] 
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)}
    
    else if (val == "subclonal") { col = colors[["mut-subclonal"]]
    rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)}
    else { col = colors[["mut-unknown"]]}
  } 
  if (stringr::str_detect(lumc_data$chr8q_CNA[i],"amplification")) {
    col = colors[["mut-clonal"]] 
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  }
}



#y = y-2
y = y-1
for (i in 1:80) {
  col = white
  
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$chr8q_CNA_clonality[i]
  val2 = substr(lumc_data$chr8q_CNA[i],1,13)
  if (!is.na(val)) {
    if (val2 == "gain") {
      if (val == "clonal") { col = colors[["clonal gain"]]
      #rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
      }
      else if (val == "subclonal") { col = colors[["subclonal gain"]]
      #rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
      }
    }
    if (val2 == "amplification") {
      if (val == "clonal") { col = colors[["mut-clonal"]]
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)}
      else if (val == "subclonal") { col = colors[["mut-subclonal"]]
      rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)}
      #text(x,y, col="white",labels="3")
    }
  } 
  #
}

y=y-1

for (i in 1:80) {
  col = white
  x = get_x(i)
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  val = lumc_data$Gaq_imbalance_clonality[i]
  if (!is.na(val)) {
    if (val == "clonal") { 
      col = colors[["mut-clonal"]]
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
    else if (val == "subclonal") { 
      col = colors[["mut-subclonal"]]
      rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
    }
    else { 
      col = "#b1b1b1"
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    }
  } 
}


###




pos = 1-0.7-2

#axis(side = 2, at = 35, labels = "GEP classifier", las=2, cex=1.0, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
#segments(pos, 40+0.4, pos, 29-0.4, lwd=1.4, col="#b1b1b1")

y_purity = 30
segments(pos, y_purity+c(1,2,3,4), pos+89, lwd=1.4, col="#eeeeee")
axis(side = 2, at = y_purity:(y_purity+4), labels = c("0%","25%","50%","75%","100%"), las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, y_purity+4.4, pos, y_purity, lwd=1.4, col="#b1b1b1")

text(-15, y_purity+2.5, labels= "Purity", col="#333333",cex=1.1, xpd=T, srt=90)


for (i in 1:80) {
  col = white
  x = get_x(i)
  purity_tumour = min(c(as.numeric(purity$Result[which(purity$Tumor_ID==lumc_data$Tumor_ID[i])])/100,1))
  rect(x-xi, y_purity, x+xi, y_purity+purity_tumour*4, border=NA, col="#777777")
  #segments(x, y_purity+as.numeric(lumc_data$ccf_min[i])/100*4,x, y_purity+as.numeric(lumc_data$ccf_max[i])/100*4, lwd=1.4, col="#777777")
}
segments(pos, y_purity, pos+89, lwd=1.4, col="#b1b1b1")

y_purity = 12
segments(pos, seq(from = y_purity, to = (y_purity+5), length.out=7), pos+89, lwd=1.4, col="#eeeeee")
axis(side = 2, at = seq(from = y_purity, to = (y_purity+5), length.out=7), labels = c("0","+1","+2","+3","+4","+5", "+6"), las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
axis(side = 2, at = seq(from = y_purity, to = (y_purity+5), length.out=7)[c(2,4,6)], labels = c("+1","+3","+5"), las=2, cex.axis=1.1, col="#333333", pos = pos, tick = 0)
segments(pos, y_purity+5.4, pos, y_purity, lwd=1.4, col="#b1b1b1")

text(-16+3, y_purity+2.5, labels= "Chr. 8q", col="#333333",cex=1.1, xpd=T, srt=90)
text(-14+3, y_purity+2.5, labels= "copy number", col="#333333",cex=1.1, xpd=T, srt=90)

chr8q_normalised = read.table("digital-pcr/chr8q-normalised.tsv", sep="\t")

for (i in 1:80) {
  col = white
  x = get_x(i)
  tumour_8q = as.numeric(chr8q_normalised[which(rownames(chr8q_normalised)==lumc_data$Tumor_ID[i]),2:4])/1
  if (chr8q_normalised$V5[which(rownames(chr8q_normalised)==lumc_data$Tumor_ID[i])] == "normal") {
  tumour_8q = c(0,0,0)
  }
  
  val = lumc_data$chr8q_CNA_clonality[i]
  val2 = substr(lumc_data$chr8q_CNA[i],1,13)
  
  col = "green"
  if (!is.na(val)) {
    if (val2 == "gain") {
      if (val == "clonal") { col = colors[["mut-clonal"]]
      #rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
      }
      else if (val == "subclonal") { col = colors[["mut-subclonal"]]
      #rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
      }
    }
    if (val2 == "amplification") {
      if (val == "clonal") { col = colors[["mut-clonal"]] }
      #rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)}
      else if (val == "subclonal") { col = colors[["mut-subclonal"]] }
      #rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)}
      #text(x,y, col="white",labels="3")
    }
  } 
  
  
  rect(x-xi, y_purity, x+xi, y_purity+tumour_8q[1]/6*5, border=NA, col=col)
  segments(x, y_purity+tumour_8q[3]/6*5, x, y_purity+tumour_8q[2]/6*5, col="#b1b1b1", lwd=1.4)
  #segments(x, y_purity+as.numeric(lumc_data$ccf_min[i])/100*4,x, y_purity+as.numeric(lumc_data$ccf_max[i])/100*4, lwd=1.4, col="#777777")
}
segments(pos, y_purity, pos+89, lwd=1.4, col="#b1b1b1")


y = 36.5
x1 = 1+0
x2 = 1+7
segments(x1-0.4,y,x2+0.4,y, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x1-0.4,y,x1-0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x2+0.4,y,x2+0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(mean(c(x1,x2)),y,mean(c(x1,x2)),y+0.7, col="#b1b1b1", lwd=1.4, xpd=T)
text(mean(c(x1,x2)), y+4.95, labels= substitute(italic("EIF1AX")), pos=1,col="#333333",cex=1.1,xpd=T)
text(mean(c(x1,x2)), y+3.85, labels= substitute("mutation"), pos=1,col="#333333",cex=1.1)
#text(mean(c(x1,x2)), y+3.85, labels= "mutant", pos=1,col="#333333",cex=1.1)
text(mean(c(x1,x2)), y+2.75, labels= "(n=8)", pos=1,col="#333333",cex=1.1)

x1 = 11
x2 = 10+14
segments(x1-0.4,y,x2+0.4,y, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x1-0.4,y,x1-0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x2+0.4,y,x2+0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(mean(c(x1,x2)),y,mean(c(x1,x2)),y+0.7, col="#b1b1b1", lwd=1.4, xpd=T)
text(mean(c(x1,x2)), y+4.95, labels= substitute(italic("SF3B1")), pos=1,col="#333333",cex=1.1,xpd=T)
text(mean(c(x1,x2)), y+3.85, labels= substitute("mutation"), pos=1,col="#333333",cex=1.1)
#text(mean(c(x1,x2)), y+3.85, labels= "mutant", pos=1,col="#333333",cex=1.1)
text(mean(c(x1,x2)), y+2.75, labels= "(n=14)", pos=1,col="#333333",cex=1.1)

x1 = 27
x2 = 73
segments(x1-0.4,y,x2+0.4,y, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x1-0.4,y,x1-0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x2+0.4,y,x2+0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(mean(c(x1,x2)),y,mean(c(x1,x2)),y+0.7, col="#b1b1b1", lwd=1.4, xpd=T)
text(mean(c(x1,x2)), y+4.95, labels= substitute("Chr. 3p"), pos=1,col="#333333",cex=1.1,xpd=T)
text(mean(c(x1,x2)), y+3.85, labels= substitute("loss"), pos=1,col="#333333",cex=1.1)
#text(mean(c(x1,x2)), y+3.85, labels= "driven", pos=1,col="#333333",cex=1.1)
text(mean(c(x1,x2)), y+2.75, labels= "(n=47)", pos=1,col="#333333",cex=1.1)


x1 = 76
x2 = 86
segments(x1-0.4,y,x2+0.4,y, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x1-0.4,y,x1-0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(x2+0.4,y,x2+0.4,y-0.7, col="#b1b1b1", lwd=1.4, xpd=T)
segments(mean(c(x1,x2)),y,mean(c(x1,x2)),y+0.7, col="#b1b1b1", lwd=1.4, xpd=T)
text(mean(c(x1,x2)), y+4.95, labels= substitute("Alternative"), pos=1,col="#333333",cex=1.1,xpd=T)
text(mean(c(x1,x2)), y+3.85, labels= substitute("evolution"), pos=1,col="#333333",cex=1.1)
#text(mean(c(x1,x2)), y+3.85, labels= "driven", pos=1,col="#333333",cex=1.1)
text(mean(c(x1,x2)), y+2.75, labels= "(n=11)", pos=1,col="#333333",cex=1.1)
#segments(33-0.4,y,84+0.4,y, col="#b1b1b1", lwd=1.4)
#segments(87,y,87+10,y, col="#b1b1b1", lwd=1.4, xpd=T)

#text(117/2, y+3.25, labels= "GEP class II", pos=1,col="#333333",cex=1.0)


dev.off()
system(paste("open", file))

xn = -22.5
y = 3
x =xn
x=xn
text(x-2.5, y, labels = "Clonality of CNAs", cex=1.0, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["subclonal loss"]], xpd=T)
text(x, y, labels = "subclonal loss", cex=1.0, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["clonal loss"]], xpd=T)
text(x, y, labels = "clonal loss", cex=1.0, col="#333333", font=1, pos=4, xpd=T)
y = y-1.5
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["subclonal gain"]], xpd=T)
text(x, y, labels = "subclonal gain", cex=1.0, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["clonal gain"]], xpd=T)
text(x, y, labels = "clonal gain", cex=1.0, col="#333333", font=1, pos=4, xpd=T)
y = y-1.5
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["subclonal amplification"]], xpd=T)
text(x, y, labels = "subclonal amplification", cex=1.0, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["clonal amplification"]], xpd=T)
text(x, y, labels = "clonal amplification", cex=1.0, col="#333333", font=1, pos=4, xpd=T)




