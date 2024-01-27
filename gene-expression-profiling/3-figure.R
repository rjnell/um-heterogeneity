###
##  Figure 1 - Landscape of transcriptional, molecular and clinical characteristics
###

# Close current images
system("taskkill /F /IM Microsoft.Photos.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)

colors = list()
colors[["female"]] = RColorBrewer::brewer.pal(12, "Paired")[5]
colors[["Class I"]] = "#777777"
colors[["Class II"]] = "#b1b1b1"
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

m = 82
bar = "darkblue"
s = 87
z = 10

file = "gene-expression-profiling/figure-1.png"
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

# Set seed
set.seed(83539)

# Load LUMC data
lumc_normalised = as.matrix(read.csv("data/normalised-80.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))

# Load classifier data
classifier = read.csv("../um-rna/gene-expression-profiling/classifier.tsv", stringsAsFactors = F, sep = "\t")

# Plot 
lumc_classifier = lumc_normalised[classifier$ENSG,]
lumc_classifier_ph = pheatmap::pheatmap(lumc_classifier, scale="row", custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                                        cutree_cols = 2, cutree_rows = 2, color = RColorBrewer::brewer.pal(9,"RdBu"), show_rownames=F,breaks = seq(from=-7,to=7,length.out = 10),
                                        silent = T)
#clustering_distance_cols = "correlation", clustering_method = "ward.D2", 


lumc_classifier = lumc_classifier[lumc_classifier_ph$tree_row$labels[lumc_classifier_ph$tree_row$order],lumc_data$Tumor_ID]
rows = which(matrixStats::rowSds(lumc_classifier) > quantile(matrixStats::rowSds(lumc_classifier), 0.05))


# Plot 
lumc_classifier = lumc_normalised[classifier$ENSG,]
lumc_classifier = lumc_classifier[which(matrixStats::rowSds(lumc_classifier)>0.5,),]
lumc_classifier_ph = pheatmap::pheatmap(lumc_classifier, 
                                        scale="row", 
                                        custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                                        cutree_cols = 2, 
                                        cutree_rows = 2,
                                        silent = T)
ord = (lumc_classifier_ph$tree_col$labels[lumc_classifier_ph$tree_col$order])
lumc_data = lumc_data[match(ord, lumc_data$Tumor_ID),]
lumc_classifier = lumc_classifier[lumc_classifier_ph$tree_row$order,match(ord, colnames(lumc_classifier))]


rows = 1:nrow(lumc_classifier)

pos=c(s,s+2.5,s+5,s+7.5, s+10)
#segments(pos, 29+0.4, pos, 29-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 14+1.4, pos, 14-0.4, lwd=1.4, col="#eeeeee", xpd=T)

y = 40.4
xi=0.5
yi = 11.8/length(rows)
palette = colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(101)
for (row in rows) {
  for (i in 1:80) {
    x = i 
    if (i > 29) { x = x + 2 }
    col = (lumc_classifier[row, i] - min(lumc_classifier[row, ])) / (max(lumc_classifier[row, ]) - min(lumc_classifier[row, ])) * 100 + 1
    rect(x-xi, y-yi/2, x+xi, y+yi/2, border=NA, col=palette[col],xpd=T)
  }
  y = y - yi
}

white = "#EEEEEE"
border = NA
xi = 0.4


y = 28

y = y-1
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$gender[i]
  if (val == "F") { col = colors[["female"]] }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
n = length(which(lumc_data$gender == "F"))
rect(s, y-0.4, s+n/m*z, y+0.4, border=border, col=colors[["female"]], xpd=T)
text(s+n/m*z, y, labels=paste0(round(n/m*100),"%"), pos=4, xpd=T)



y = y-1
n1 = 0
n2 = 0
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$melanoma_death[i]
  if (val == "yes") { n1=n1+1;col = "#777777"}
  #if (val == "death by unknown") { n2=n2+1;col = "#B1b1b1"}
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
rect(s, y-0.4, s+n1/m*z, y+0.4, border=NA, col="#777777", xpd=T)
rect(s+n1/m*z, y-0.4, s+(n1+n2)/m*z, y+0.4, border=NA, col="#B1b1b1", xpd=T)
rect(s, y-0.4, s+(n1+n2)/m*z, y+0.4, border=border , xpd=T)
text(s+(n1+n2)/m*z, y, labels=paste0(round((n1+n2)/m*100),"%"), pos=4, xpd=T)



y = 25

y = y-1
n1 = 0
n2 = 0
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$Gaq_mutation[i]
  if (!is.na(val)) {
    if (str_detect(val, "GNAQ p.Q209")) { n1 = n1 +1; col = RColorBrewer::brewer.pal(9,"Paired")[2]}
    else if (str_detect(val, "GNAQ p.G48")) { n2 = n2+1; col = RColorBrewer::brewer.pal(9,"Paired")[9]}
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
rect(s, y-0.4, s+n1/m*z, y+0.4, border=NA, col=RColorBrewer::brewer.pal(9,"Paired")[2], xpd=T)
rect(s+n1/m*z, y-0.4, s+(n1+n2)/m*z, y+0.4, border=NA, col=RColorBrewer::brewer.pal(9,"Paired")[9], xpd=T)
rect(s, y-0.4, s+(n1+n2)/m*z, y+0.4, border=border , xpd=T)
text(s+(n1+n2)/m*z, y, labels=paste0(round((n1+n2)/m*100),"%"), pos=4, xpd=T)

y = y-1
n1 = 0
n2 = 0
for (i in 1:82) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$Gaq_mutation[i]
  if (!is.na(val)) {
    if (str_detect(val, "GNA11 p.Q209")) { n1 = n1+1; col = RColorBrewer::brewer.pal(9,"Paired")[2]}
    else if (str_detect(val, "GNA11 p.R183")) { n2 = n2+1; col = RColorBrewer::brewer.pal(9,"Paired")[1]}
  } 
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
rect(s, y-0.4, s+n1/m*z, y+0.4, border=NA, col=RColorBrewer::brewer.pal(9,"Paired")[2], xpd=T)
rect(s+n1/m*z, y-0.4, s+(n1+n2)/m*z, y+0.4, border=NA, col=RColorBrewer::brewer.pal(9,"Paired")[1], xpd=T)
rect(s, y-0.4, s+(n1+n2)/m*z, y+0.4, border=border , xpd=T)
text(s+(n1+n2)/m*z, y, labels=paste0(round((n1+n2)/m*100),"%"), pos=4, xpd=T)


y = y-1
for (i in 1:82) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$Gaq_mutation[i]
  if (!is.na(val)) {
    if (str_detect(val, "PLCB4")) { col = colors[["PLCB4"]]}
  } 
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
n = length(which(str_detect(lumc_data$Gaq_mutation,"PLCB4")))
rect(s, y-0.4, s+n/m*z, y+0.4, border=border, col=colors[["PLCB4"]], xpd=T)
text(s+n/m*z, y, labels=paste0(round(n/m*100),"%"), pos=4, xpd=T)

y = y-1
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$Gaq_mutation[i]
  if (!is.na(val)) {
    if (str_detect(val, "CYSLTR2")) { col = colors[["CYSLTR2"]]}
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
n = length(which(str_detect(lumc_data$Gaq_mutation,"CYSLTR2")))
rect(s, y-0.4, s+n/m*z, y+0.4, border=border, col=colors[["CYSLTR2"]], xpd=T)
text(s+n/m*z, y, labels=paste0(round(n/m*100),"%"), pos=4, xpd=T)

#rect(1-0.6, y+3.5, 82+0.6, y-0.5)

y = y-2
n1 = 0
n2 = 0
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$EIF1AX_mutation[i]
  if (is.na(val)) { val = "n/a" }
  if (str_detect(val, "p.P2") | str_detect(val, "p.K3") | str_detect(val, "p.N4")) { n1 = n1+1; col = colors[["EIF1AX-exon1"]]}
  else if (str_detect(val, "splicing") | str_detect(val, "p.G8") | str_detect(val, "p.G9")) { n2 = n2+1; col = colors[["EIF1AX-exon2"]]}
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
rect(s, y-0.4, s+n1/m*z, y+0.4, border=NA, col=colors[["EIF1AX-exon1"]], xpd=T)
rect(s+n1/m*z, y-0.4, s+(n1+n2)/m*z, y+0.4, border=NA, col=colors[["EIF1AX-exon2"]], xpd=T)
rect(s, y-0.4, s+(n1+n2)/m*z, y+0.4, border=border , xpd=T)
text(s+(n1+n2)/m*z, y, labels=paste0(round((n1+n2)/m*100),"%"), pos=4, xpd=T)


y = y-1 
n1 = 0
n2 = 0
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$SF3B1_mutation[i]
  if (is.na(val)) { val = "n/a" }
  if (str_detect(val, "p.R625")) { n1 = n1+1; col = colors[["SF3B1-hotspot"]] }
  else if (val != "n/a") { n2 = n2+1; col = colors[["SF3B1-other"]]}
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
rect(s, y-0.4, s+n1/m*z, y+0.4, border=NA, col=colors[["SF3B1-hotspot"]], xpd=T)
rect(s+n1/m*z, y-0.4, s+(n1+n2)/m*z, y+0.4, border=NA, col=colors[["SF3B1-other"]], xpd=T)
rect(s, y-0.4, s+(n1+n2)/m*z, y+0.4, border=border , xpd=T)
text(s+(n1+n2)/m*z, y, labels=paste0(round((n1+n2)/m*100),"%"), pos=4, xpd=T)



y = y-1
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$BAP1_mutation[i]
  if (is.na(val)) { val = "n/a" }
  if (val == "mutant") { col = colors[["BAP1"]]}
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
n = length(which(!is.na(lumc_data$BAP1_mutation)))
rect(s, y-0.4, s+n/m*z, y+0.4, border=border, col=colors[["BAP1"]], xpd=T)
text(s+n/m*z, y, labels=paste0(round(n/m*100),"%"), pos=4, xpd=T)



y = y-2
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$chr3p_CNA[i]
  if (stringr::str_detect(val,"loss")) { col = colors[["loss"]]}
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
n = length(which(stringr::str_detect(lumc_data$chr3p_CNA,"loss")))
rect(s, y-0.4, s+n/m*z, y+0.4, border=border, col=colors[["loss"]], xpd=T)
text(s+n/m*z, y, labels=paste0(round(n/m*100),"%"), pos=4, xpd=T)



y = y-1
n1=0
n2=0
for (i in 1:80) {
  col = white
  x = i 
  if (i > 29) { x = x + 2 }
  val = lumc_data$chr8q_CNA[i]
  if (stringr::str_detect(val,"gain")) { n1=n1+1;col = colors[["gain"]]}
  else if (stringr::str_detect(val,"amplification")) { n2=n2+1;col = colors[["amp"]]}
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
rect(s, y-0.4, s+n1/m*z, y+0.4, border=NA, col=colors[["gain"]], xpd=T)
rect(s+n1/m*z, y-0.4, s+(n1+n2)/m*z, y+0.4, border=NA, col=colors[["amp"]], xpd=T)
rect(s, y-0.4, s+(n1+n2)/m*z, y+0.4, border=border , xpd=T)
text(s+(n1+n2)/m*z, y, labels=paste0(round((n1+n2)/m*100),"%"), pos=4, xpd=T)

pos = 1-0.7-2

axis(side = 2, at = 35, labels = "GEP classifier", las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, 40+0.4, pos, 29-0.4, lwd=1.4, col="#b1b1b1",xpd=T)


axis(side = 2, at = 27:26, labels = c("Gender", "UM death"), las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)


axis(side = 2, at = 24:21, labels = c("GNAQ", "GNA11", "PLCB4", "CYSLTR2"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)



axis(side = 2, at = 19:17, labels = c("EIF1AX","SF3B1","BAP1"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
axis(side = 2, at = 15:14, labels = c("Chromosome 3p","Chromosome 8q"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)

#pos = c(pos, s)
segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 14+1.4, pos, 14-0.4, lwd=1.4, col="#b1b1b1")



y = 12.4
axis(side = 1, at = c(s,s+2.5,s+5,s+7.5, s+10), labels= c("0%", NA, NA, NA, "100%") , xpd=T, pos = y, cex=1.1, col="#333333",  col.ticks = "#b1b1b1", lwd.ticks = 1.4)
axis(side = 1, at = 30/2, labels= c("GEP class I") , xpd=T, pos = y, cex.axis=1.1, col="#333333",  col.ticks = "#b1b1b1", lwd.ticks = 1.4)
axis(side = 1, at = 114/2, labels= c("GEP class II") , xpd=T, pos = y, cex.axis=1.1, col="#333333",  col.ticks = "#b1b1b1", lwd.ticks = 1.4)

segments(1-0.4,y,29+0.4,y, col="#b1b1b1", lwd=1.4)
segments(32-0.4,y,82+0.4,y, col="#b1b1b1", lwd=1.4)
segments(87,y,87+10,y, col="#b1b1b1", lwd=1.4, xpd=T)
text(30/2, y-1.75, labels= "(n=29)", pos=1,col="#333333",cex=1.1)
text(114/2, y-1.75, labels= "(n=51)", pos=1,col="#333333",cex=1.1)
text(92, y-1.75, labels= "frequency", pos=1,col="#333333",cex=1.1,xpd=T)

xn = -22.5
y = 7
x =xn
text(x-2.5, y, labels = "GEP classifier", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1 
text(x-2.5, y, labels = "low", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
text(x+25.5, y, labels = "high", cex=1.1, col="#333333", font=1, pos=2, xpd=T)

palette = colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(101)
for (k in 1:101) {
  rect(x+5+12/101*(k-1), y-0.4, x+5+12/101*(k), y+0.4, border=NA, col=palette[k],xpd=T)
}
rect(x+5, y-0.4, x+17, y+0.4, border=border, xpd=T)
text(x+11, y-1, labels = "expression", cex=1.1, col="#333333", font=1, xpd=T)




y =y -3
x =xn
text(x-2.5, y, labels = "Gender", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["female"]], xpd=T)
text(x, y, labels = "female", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x = x + 15
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#eeeeee", xpd=T)
text(x, y, labels = "male", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


y =y -2
x =xn
text(x-2.5, y, labels = "Death (metastatic UM)", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#y = y-1
#text(x-2.5, y, labels = "", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["Class I"]], xpd=T)
text(x, y, labels = "yes", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x = x+ 15
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["Class II"]], xpd=T)
text(x, y, labels = "no", cex=1.1, col="#333333", font=1, pos=4, xpd=T)



yn = 7
y = yn
xn = 19
x = xn
text(x-2.5, y, labels = "GNAQ/GNA11", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=RColorBrewer::brewer.pal(9,"Paired")[2])
text(x, y, labels = "p.Q209", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x = x + 16
#y=y-1
#x = xn
#y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=RColorBrewer::brewer.pal(9,"Paired")[1])
text(x, y, labels = "p.R183", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=RColorBrewer::brewer.pal(9,"Paired")[9])
text(x, y, labels = "p.G48", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

#y = yn-3
#text(x-2.5, y, labels = "GNA11", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
#y = y-1
#rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=RColorBrewer::brewer.pal(9,"Paired")[2])
#text(x, y, labels = "p.Q209", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#x = x + 16
#y=y-1
#rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=RColorBrewer::brewer.pal(9,"Paired")[1])
#text(x, y, labels = "p.R183", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#x = x + 16
#y=y+1
#rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=RColorBrewer::brewer.pal(9,"Paired")[9])
#text(x, y, labels = "p.G48", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


y = yn-4
x = xn
text(x-2.5, y, labels = "PLCB4", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["PLCB4"]], xpd=T)
text(x, y, labels = "p.D630", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y = yn-7
x = xn
text(x-2.5, y, labels = "CYSLTR2", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["CYSLTR2"]], xpd=T)
text(x, y, labels = "p.L129", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y = yn
xn = 59
x = xn
text(x-2.5, y, labels = "EIF1AX", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["EIF1AX-exon1"]], xpd=T)
text(x, y, labels = "exon 1", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["EIF1AX-exon2"]], xpd=T)
text(x, y, labels = "exon 2", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y = yn - 4
x = xn
text(x-2.5, y, labels = "SF3B1", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["SF3B1-hotspot"]], xpd=T)
text(x, y, labels = "p.R625", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["SF3B1-other"]], xpd=T)
text(x, y, labels = "other", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y = yn - 8
x= xn
text(x-2.5, y, labels = "BAP1", cex=1.1, col="#333333", font=3, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["BAP1"]], xpd=T)
text(x, y, labels = "mutant", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y = yn
xn= 82
x=xn
text(x-2.5, y, labels = "Chromosome 3p", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["loss"]], xpd=T)
text(x, y, labels = "loss", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y=y-2
text(x-2.5, y, labels = "Chromosome 8q", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["gain"]], xpd=T)
text(x, y, labels = "gain", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["amp"]], xpd=T)
text(x, y, labels = "amplification", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

dev.off()
system(paste("open", file))






# Load dependencies (pheatmap)
source("../um-rna/bin/pheatmap.r")
library(grid)
library(RColorBrewer)
library(scales)
library(gtable)
library(stats)
library(grDevices)
library(graphics)

# Plot and save dendrogram
file = "gene-expression-profiling/figure-1-dendrogram.png"
png(file, res=600, width=5000, height=1500)
dendrogram = pheatmap(lumc_classifier, 
                      scale="row", 
                      custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                      cutree_cols = 2, 
                      cutree_rows = 2,
                      treeheight_col=35)
dev.off()
system(paste("open", file))

