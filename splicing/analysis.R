# Load dependencies (pheatmap)
source("../um-rna/bin/pheatmap.R")
library(grid)
library(RColorBrewer)
library(scales)
library(gtable)
library(stats)
library(grDevices)
library(graphics)

#lumc_splicing = data.table::fread("data/splicing-normalised-80.tsv", sep="\t")
signature = read.table("../um-rna/data/splice-junctions-paris-GRCh38.tsv", sep="\t")

w_low = which(signature$Log2FC<0)
selection_low = which(lumc_splicing$V1 %in% 
                        paste0("chr",signature$chr_38,"-",signature$pos1_38,"-",signature$pos2_38)[w_low]
                      |
                        lumc_splicing$V1 %in% 
                        paste0("chr",signature$chr_38,"-",signature$pos1_38-1,"-",signature$pos2_38-1)[w_low])


w_high = which(signature$Log2FC>0)
selection_high = which(lumc_splicing$V1 %in% 
                         paste0("chr",signature$chr_38,"-",signature$pos1_38,"-",signature$pos2_38)[w_high]
                       |
                         lumc_splicing$V1 %in% 
                         paste0("chr",signature$chr_38,"-",signature$pos1_38-1,"-",signature$pos2_38-1)[w_high])

data_selection = as.matrix(lumc_splicing[c(selection_low,selection_high),2:81])
rownames(data_selection) = lumc_splicing$V1[c(selection_low,selection_high)]
annotation = c(rep("low", length(selection_low)), rep("high", length(selection_high)))

data_sds = matrixStats::rowSds(data_selection)
data_selection = data_selection[which(data_sds>quantile(data_sds,0.1)),]

# Plot and save dendrogram
file = "splicing/figure-dendrogram.png"
png(file, res=600, width=5000, height=1500)
dendrogram = pheatmap(data_selection, 
                      scale='row',
                      treeheight_col = 35,
                      custering_callback = function(hc, ...) { 
                        dendsort::dendsort(hc, type = "min") 
                      })
dev.off()
system(paste("open", file))

###

# Plot figure

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# Load LUMC data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)
lumc_data = lumc_data[match(dendrogram$tree_col$labels[dendrogram$tree_col$order], lumc_data$Tumor_ID),]

colors = list()
colors[["DNA + RNA"]] = "#999999"
colors[["RNA"]] = "#2196F3"
colors[["DNA"]] = "#F44336"

m = 82
bar = "darkblue"
s = 87
z = 10

i_stop = 63

file = "splicing/figure.png"
png(file, res=600, width=5000, height=7500)
par(mar=c(5,3,5,7))
plot(c(-15,95), 
     c(-7,43),
     type="n", 
     axes=F, 
     xaxs="i", 
     yaxs="i", 
     xlab="", 
     ylab="") 

# Set seed
set.seed(83539)

# Plot 
lumc_classifier = data_selection[,lumc_data$Tumor_ID]
rows = dendrogram$tree_row$labels[dendrogram$tree_row$order]


y = 40.4+3
xi=0.5
yi = 9.8/length(rows)
palette = (colorRampPalette(RColorBrewer::brewer.pal(9,"BuPu")[1:7])(101))
for (row in rows) {
  for (i in 1:80) {
    x = i 
    if (i > i_stop) { x = x + 2 }
    col = (lumc_classifier[row, i] - min(lumc_classifier[row, ])) / (max(lumc_classifier[row, ]) - min(lumc_classifier[row, ])) * 100 + 1
    rect(x-xi, y-yi/2, x+xi, y+yi/2, border=NA, col=palette[col],xpd=T)
  }
  #if (row == 206) {
  #y=y-0.1
  #}
  y = y - yi
}

white = "#EEEEEE"
border = NA
xi = 0.4

y = 33
for (i in 1:80) {
  col = white
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = lumc_data$SF3B1_mutation[i]
  if (!is.na(val)) {
    col = colors[["DNA + RNA"]]
    if (val == "p.A638S (c.G1912T)") {
      col = "red"
    }
    else if (val == "p.R625C (c.C1873T)") {
      col = "#E67E22"
    }
    else if (val == "p.G664delinsEC (c.1990_1991insAAT)") {
      col = "brown"
    }
    else if (val == "p.R625H (c.G1874A)") {
      col = "#E67E22"
    }
    else if (val == "p.E622D (c.G1866C)") {
      col = "#F0B27A"
    }
    else if (val == "p.K666T (c.A1997C)") {
      col = "gold"
    }
  }
  else { 
    col = white
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}

pos = 1-0.7-2

limits=c(40+3,34)
axis(side = 2, at = mean(limits), labels = "Splicing", las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, limits[1]+0.4, pos, limits[2]-0.4, lwd=1.4, col="#b1b1b1", xpd=T)

limits=c(33,33)
axis(side = 2, at = mean(limits), labels = "SF3B1", las=2, font=3, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, limits[1]+0.4, pos, limits[2]-0.4, lwd=1.4, col="#b1b1b1", xpd=T)

xn = 90
y = 40+3
x =xn
text(x-2.5, y, labels = "Splicing", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1 
text(x-2.5, y, labels = "low", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
text(x+25.5, y, labels = "high", cex=1.1, col="#333333", font=1, pos=2, xpd=T)

palette = colorRampPalette(RColorBrewer::brewer.pal(9,"BuPu"))(101)
for (k in 1:101) {
  rect(x+5+12/101*(k-1), y-0.4, x+5+12/101*(k), y+0.4, border=NA, col=palette[k],xpd=T)
}
rect(x+5, y-0.4, x+17, y+0.4, border=border, xpd=T)
text(x+11, y-1, labels = "expression", cex=1.1, col="#333333", font=1, xpd=T)


yn = 36+3
y = yn
xn = 90
x = xn
text(x-2.5, y, labels = substitute(italic("SF3B1")~"mutation"), cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#E67E22")
text(x, y-0.05, labels = "p.R625C/H", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#F0B27A")
text(x, y-0.05, labels = "p.E622D", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="brown")
text(x, y-0.05, labels = "p.G664delinsEC", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="gold")
text(x, y-0.05, labels = "p.K666T", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="red")
text(x, y-0.05, labels = "p.A638S", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white)
text(x, y-0.05, labels = "not detected", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


dev.off()
system(paste("open", file))



