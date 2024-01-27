yshift = 1

plot_point = function(pos,n,x,label,col=RColorBrewer::brewer.pal(9,"Paired")[6]) {
  segments(pos, 0, pos, n/ymax*5+yshift, lwd=1.4,col="#B1b1b1",xpd=T)
  points(pos, n/ymax*5+yshift, cex=0.9, col=RColorBrewer::brewer.pal(9,"Paired")[6], pch=16,xpd=T)
  text(pos+x, n/ymax*5+yshift+.75, labels=label, col="#333333",xpd=T, srt=0, cex=0.8)
}

plot_point_inv = function(pos,n,x,label,col=RColorBrewer::brewer.pal(9,"Paired")[6]) {
  segments(pos, 0, pos, n/ymax*-5-yshift, lwd=1.4,col="#B1b1b1",xpd=T)
  points(pos, n/ymax*-5-yshift, cex=0.9, col=RColorBrewer::brewer.pal(9,"Paired")[6], pch=16,xpd=T)
  text(pos+x, n/ymax*-5-yshift-.75, labels=label, col="#333333",xpd=T, srt=0, cex=0.75)
}

plot_point = function(pos,n,x,label,col=RColorBrewer::brewer.pal(9,"Paired")[2],srt=90) {
  segments(pos, 0, pos, 0.75, lwd=1.4,col="#B1b1b1",xpd=T)
  segments(pos, 0.75, pos+x, 1, lwd=1.4,col="#B1b1b1",xpd=T)
  segments(pos+x, 1, pos+x, n/ymax*5+yshift, lwd=1.4,col="#B1b1b1",xpd=T)
  points(pos+x, n/ymax*5+yshift, cex=1, col=col, pch=16,xpd=T)
  if (srt == 0) {
    text(pos+x, n/ymax*5+yshift+.75, labels=label, col="#333333",xpd=T, srt=0, cex=1.1)
  }
  else {
    text(pos+x, n/ymax*5+yshift+0.25, labels=label, col="#333333",xpd=T, srt=90, adj=0, cex=1)
  }
}

### Settings
system("taskkill /F /IM Microsoft.Photos.exe /T")
graphics.off()
gene = "SF3B1"
length = 1304
ymax = 14
yby = 7

# Init plot

file = paste0("mutations/",gene,".png")
png(file,res=600,width=3000,height=2750)
par(mar=c(5,7,5,7))
xlim = c(0,length)
plot(xlim,c(-3,6.5),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i")

# Plot y axis
yat = seq(0, 5, length.out = yby)+yshift
axis(side = 2, at = yat, labels=c(0,rep("",yby-2),ymax),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 1)
mtext(side = 2, text = "mutations in", line = 4.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "LUMC cohort", line = 3.5, at = mean(yat), cex=1.1)

### Points
xs = xlim[2]/25
#plot_point(pos = 48, n = 1, x = 0*xs, label = "G48L")
#plot_point(pos = 209, n = 33, x = 0*xs, label = "Q209P/L")

plot_point(pos = 622, n = 1, x = -1*xs, label = "p.E622D", col="#B1B1b1")
plot_point(pos = 625, n = 14, x = 0*xs, label = "p.R625C/H", srt=0, col="#2196F3")
plot_point(pos = 638, n = 1, x = +1*xs, label = "p.A638S", col="#B1B1b1")
plot_point(pos = 666, n = 1, x = +3.75*xs, label = "p.K666T", col=RColorBrewer::brewer.pal(9,"Paired")[8])
plot_point(pos = 664, n = 1, x = +2.25*xs, label = "p.G664delinsEC", col=RColorBrewer::brewer.pal(9,"Paired")[8])


#ymax = 37
#yby = 7
# Plot y axis
#yat = seq(0, 5, length.out = yby)+yshift
#axis(side = 2, at = -yat, labels=c(0,rep("",yby-2),ymax),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex=0.9, line = 1)
#mtext(side = 2, text = "mutations", line = 3.5, at = mean(-yat))

#plot_point_inv(pos = 209, n = 37, x = 0*xs, label = "Q209P/L")
#plot_point_inv(pos = 183, n = 2, x = -1*xs, label = "R183Q")
#plot_point_inv(pos = 48, n = 1, x = 0*xs, label = "G48L")


# Basic bar
rect(xlim[1],-0.25,xlim[2],0.25,xpd=T,border=NA,col="#B1b1b1")

# Item
pos = c(400,1300)
rect (pos[1],-0.5,pos[2],0.5,xpd=T,border=NA,col=RColorBrewer::brewer.pal(9,"Paired")[4])
text(mean(pos),0,labels="HEAT domains", col = "white", cex=1.)

# Item
pos = c(225,310)
rect (pos[1],-0.5,pos[2],0.5,xpd=T,border=NA,col=RColorBrewer::brewer.pal(9,"Paired")[6])
text(mean(pos),0,labels="S", col = "white", cex=1.1)

# Axis
pos = -1.25
xat = c(seq(xlim[1], xlim[2], by=500), xlim[2])
segments(xlim[1],pos,xlim[2], col = "#b1b1b1", lwd = 1.4,xpd=T)
axis(side = 1, pos = pos,at = xat[1:(length(xat)-1)], labels = xat[1:(length(xat)-1)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
axis(side = 1, pos = pos,at = xat[length(xat)], labels = paste0("  ",xat[length(xat)],"aa"), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
#text(length/2, -6.5, labels=gene, font=4, xpd=T)

# Axis
#text(xlim[2], 0, labels=paste0(xlim[2],"aa"), pos=4, xpd=T, cex=0.8, col="#333333")
#text(xlim[2]/2, 9, labels=gene, font=4, xpd=T, col="#333333")

#text(xlim[2]*1.055, 4.7, pos=4, labels="Digital PCR", cex=1.1, xpd=T, font=2)
#points(xlim[2]*1.123, 3.9, cex=1, col=RColorBrewer::brewer.pal(9,"Paired")[2], pch=16,xpd=T)
#text(xlim[2]*1.1, 3.9, pos=4, labels="targeted assay", cex=1.1, xpd=T)
#points(xlim[2]*1.123, 3.1, cex=1, col=RColorBrewer::brewer.pal(9,"Paired")[8], pch=16,xpd=T)
#text(xlim[2]*1.1, 3, pos=4, labels="drop-off assay", cex=1.1, xpd=T)
#points(xlim[2]*1.123, 2.3, cex=1, col="#b1b1b1", pch=16,xpd=T)
#text(xlim[2]*1.1, 2.1, pos=4, labels="not analysed", cex=1.1, xpd=T)

dev.off()
system(paste("open", file))












### Settings
system("taskkill /F /IM Microsoft.Photos.exe /T")
graphics.off()
gene = "EIF1AX"
length = 144
ymax = 5
yby = 6

# Init plot

file = paste0("mutations/",gene,".png")
png(file,res=600,width=3000,height=2750)
par(mar=c(5,7,5,7))
xlim = c(0,length)
plot(xlim,c(-3,6.5),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i")

# Plot y axis
yat = seq(0, 5, length.out = yby)+yshift
axis(side = 2, at = yat, labels=c(0,rep("",yby-2),ymax),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 1)
mtext(side = 2, text = "mutations in", line = 4.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "LUMC cohort", line = 3.5, at = mean(yat), cex=1.1)

### Points
xs = xlim[2]/25*1.4
#plot_point(pos = 48, n = 1, x = 0*xs, label = "G48L")
#plot_point(pos = 209, n = 33, x = 0*xs, label = "Q209P/L")

#plot_point(pos = 622, n = 1, x = -1*xs, label = "E622D", col="#B1B1b1")
#plot_point(pos = 625, n = 14, x = 0*xs, label = "R625C/H", srt=0)
#plot_point(pos = 638, n = 1, x = +1*xs, label = "A638S", col="#B1B1b1")
#plot_point(pos = 666, n = 1, x = +3.5*xs, label = "K666T", col=RColorBrewer::brewer.pal(9,"Paired")[8])
plot_point(pos = 2, n = 3, x = -0.5*xs, label = "p.P2S/R/L", col=RColorBrewer::brewer.pal(9,"Paired")[8])
plot_point(pos = 3, n = 1, x = 0.5*xs, label = "p.K3N", col=RColorBrewer::brewer.pal(9,"Paired")[8])
plot_point(pos = 4, n = 1, x = +1.5*xs, label = "p.N4S", col=RColorBrewer::brewer.pal(9,"Paired")[8])

plot_point(pos = 6, n = 1, x = +3*xs, label = "p.G6splice", col=RColorBrewer::brewer.pal(9,"Paired")[7])
plot_point(pos = 8, n = 2, x = +4*xs, label = "p.G8R", col=RColorBrewer::brewer.pal(9,"Paired")[7])
plot_point(pos = 9, n = 3, x = +5*xs, label = "p.G9R/D", col=RColorBrewer::brewer.pal(9,"Paired")[7])

#ymax = 37
#yby = 7
# Plot y axis
#yat = seq(0, 5, length.out = yby)+yshift
#axis(side = 2, at = -yat, labels=c(0,rep("",yby-2),ymax),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex=0.9, line = 1)
#mtext(side = 2, text = "mutations", line = 3.5, at = mean(-yat))

#plot_point_inv(pos = 209, n = 37, x = 0*xs, label = "Q209P/L")
#plot_point_inv(pos = 183, n = 2, x = -1*xs, label = "R183Q")
#plot_point_inv(pos = 48, n = 1, x = 0*xs, label = "G48L")


# Basic bar
rect(xlim[1],-0.25,xlim[2],0.25,xpd=T,border=NA,col="#B1b1b1")

# Item
pos = c(32,94)
rect (pos[1],-0.5,pos[2],0.5,xpd=T,border=NA,col=RColorBrewer::brewer.pal(9,"Paired")[4])
text(mean(pos),0,labels="RNA binding", col = "white", cex=1.)


# Axis
pos = -1.25
xat = c(seq(xlim[1], xlim[2], by=50), xlim[2])
segments(xlim[1],pos,xlim[2], col = "#b1b1b1", lwd = 1.4,xpd=T)
axis(side = 1, pos = pos,at = xat[1:(length(xat)-1)], labels = xat[1:(length(xat)-1)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
axis(side = 1, pos = pos,at = xat[length(xat)], labels = paste0(" ",xat[length(xat)],"aa"), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
#text(length/2, -6.5, labels=gene, font=4, xpd=T)

# Axis
#text(xlim[2], 0, labels=paste0(xlim[2],"aa"), pos=4, xpd=T, cex=0.8, col="#333333")
#text(xlim[2]/2, 9, labels=gene, font=4, xpd=T, col="#333333")
n = 2+1
#points(xlim[2]*1.123, 6-(5-(n*0.8))/2-0.4-0.8*(  1  ), cex=1, col=RColorBrewer::brewer.pal(9,"Paired")[8], pch=16,xpd=T)
#points(xlim[2]*1.123, 6-(5-(n*0.8))/2-0.4-0.8*(  2  ), cex=1, col=RColorBrewer::brewer.pal(9,"Paired")[7], pch=16,xpd=T)
#points(xlim[2]*1.123, 6-(5-(n*0.8))/2-0.4-0.8*(  3  ), cex=1, col="#b1b1b1", pch=16,xpd=T)

dev.off()
system(paste("open", file))































### Settings
system("taskkill /F /IM Microsoft.Photos.exe /T")
graphics.off()
gene = "BAP1"
length = 729
ymax = 5
yby = 5

# Init plot

file = paste0("mutations/",gene,".png")
png(file,res=600,width=3000,height=2750)
par(mar=c(5,7,5,7))
xlim = c(0,length)
plot(xlim,c(-3,6.5),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i")

# Plot y axis
yat = seq(0, 5, length.out = yby)+yshift
axis(side = 2, at = yat, labels=c(0,rep("",yby-2),ymax),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 1)
mtext(side = 2, text = "mutations in", line = 4.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "LUMC cohort", line = 3.5, at = mean(yat), cex=1.1)

### Points
xs = xlim[2]/25*1.4
#plot_point(pos = 48, n = 1, x = 0*xs, label = "G48L")
#plot_point(pos = 209, n = 33, x = 0*xs, label = "Q209P/L")

#plot_point(pos = 622, n = 1, x = -1*xs, label = "E622D", col="#B1B1b1")
#plot_point(pos = 625, n = 14, x = 0*xs, label = "R625C/H", srt=0)
#plot_point(pos = 638, n = 1, x = +1*xs, label = "A638S", col="#B1B1b1")
#plot_point(pos = 666, n = 1, x = +3.5*xs, label = "K666T", col=RColorBrewer::brewer.pal(9,"Paired")[8])

nr = as.numeric(substr(lumc_data$`BAP1-position`[which(!is.na(lumc_data$`BAP1-position`) & lumc_data$`BAP1-position`!="exon05")],3,6))
tt = table(nr)
for (i in 1:length(tt)) {
  col = "grey"
  if (as.numeric(names(tt[i])) %in% c(88,90,97)) {
    col = "#E67E22"
  }
  else if (as.numeric(names(tt[i])) %in% c(68,126,161,184,185,194,223,226,344,346,631)) {
    col = "#2196F3"
  }
  xxx = 0*xs
  if (as.numeric(names(tt[i])) == 23) {
    #xxx = -0.1*xs
  }
  if (as.numeric(names(tt[i])) == 24) {
    #xxx = +0.1*xs
  }
  
  if (as.numeric(names(tt[i])) == 42) {
    #xxx = -0.1*xs
  }
  if (as.numeric(names(tt[i])) == 43) {
    #xxx = +0.1*xs
  }
  
  if (as.numeric(names(tt[i])) == 86) {
    #xxx = -0.1*xs
  }
  if (as.numeric(names(tt[i])) == 90) {
    #xxx = +0.1*xs
  }
  plot_point(pos = as.numeric(names(tt[i])), n = tt[i], x = xxx, label="", col=col)
}


as.numeric(names(tt))

#text(344, 1/ymax*5+yshift+0.25, labels="p.N344Mfs*17", col="#333333",xpd=T, srt=90, adj=0, cex=1)
#text(631, 1/ymax*5+yshift+0.25, labels="p.E631X", col="#333333",xpd=T, srt=90, adj=0, cex=1)

#ymax = 37
#yby = 7
# Plot y axis
#yat = seq(0, 5, length.out = yby)+yshift
#axis(side = 2, at = -yat, labels=c(0,rep("",yby-2),ymax),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex=0.9, line = 1)
#mtext(side = 2, text = "mutations", line = 3.5, at = mean(-yat))

#plot_point_inv(pos = 209, n = 37, x = 0*xs, label = "Q209P/L")
#plot_point_inv(pos = 183, n = 2, x = -1*xs, label = "R183Q")
#plot_point_inv(pos = 48, n = 1, x = 0*xs, label = "G48L")


# Basic bar
rect(xlim[1],-0.25,xlim[2],0.25,xpd=T,border=NA,col="#B1b1b1")

# Item
pos = c(1,240)
rect (pos[1],-0.5,pos[2],0.5,xpd=T,border=NA,col=RColorBrewer::brewer.pal(9,"Paired")[4])
text(mean(pos),0,labels="UCH", col = "white", cex=1.)
pos = c(594,657)
rect (pos[1],-0.5,pos[2],0.5,xpd=T,border=NA,col=RColorBrewer::brewer.pal(9,"Paired")[6])
text(mean(pos),0,labels="B", col = "white", cex=1.)

# Axis
pos = -1.25
xat = c(seq(xlim[1], xlim[2], by=250), xlim[2])
segments(xlim[1],pos,xlim[2], col = "#b1b1b1", lwd = 1.4,xpd=T)
axis(side = 1, pos = pos,at = xat[1:(length(xat)-1)], labels = xat[1:(length(xat)-1)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
axis(side = 1, pos = pos,at = xat[length(xat)], labels = paste0("  ",xat[length(xat)],"aa"), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
#text(length/2, -6.5, labels=gene, font=4, xpd=T)

# Axis
#text(xlim[2], 0, labels=paste0(xlim[2],"aa"), pos=4, xpd=T, cex=0.8, col="#333333")
#text(xlim[2]/2, 9, labels=gene, font=4, xpd=T, col="#333333")
n = 2+1
#points(xlim[2]*1.123, 6-(5-(n*0.8))/2-0.4-0.8*(  1  ), cex=1, col=RColorBrewer::brewer.pal(9,"Paired")[8], pch=16,xpd=T)
#points(xlim[2]*1.123, 6-(5-(n*0.8))/2-0.4-0.8*(  2  ), cex=1, col=RColorBrewer::brewer.pal(9,"Paired")[7], pch=16,xpd=T)
#points(xlim[2]*1.123, 6-(5-(n*0.8))/2-0.4-0.8*(  3  ), cex=1, col="#b1b1b1", pch=16,xpd=T)

dev.off()
system(paste("open", file))