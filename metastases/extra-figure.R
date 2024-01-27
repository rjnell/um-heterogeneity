
col_subclonal = "#F88E86"
col_clonal = "#D4190C"


# Load data
id_data = NULL

# LUMC-79-mem
id_data = rbind(id_data,
                c("GNA11 ...", 54.2, 31.3, 77.0, "#D4190C"),
                c("Chr.3p loss", 0, 0, 0, "#D4190C"),
                c("Chr.8q amplification", 0,0,0, "#D4190C")  #5.11, 3.55, 9.10
)
id_data = NULL

id = "LUMC-22"
id2 = "metastatic lesion"

# LUMC-2-mem
id_data = rbind(id_data,
                c("GNAQ p.Q209L", 89.4, 86.1, 92.8, "#D4190C"),
                c("SF3B1 p.R625H", 91.5, 87.7, 95.4, "#D4190C"),
                c("Chr.8q gain", 16,9,23, "#93D095")  #5.11, 3.55, 9.10
)

id = "LUMC-22"
id2 = "primary tumour"
id_data = NULL
id_data = rbind(id_data,
                c("GNAQ p.Q209L (c.A626T)","96.915181054773","92.546724362451","101.28477674077","#D4190C"),
                  c("SF3B1 p.R625H (c.G1874A)","98.426429788867","94.318382904742","102.53503142087","#D4190C"),
c("Chr.8q gain","45.4583169","35.3754112","56.4233097","#7AC0F8"))
  



{
  file = paste0("clonality/extra/",id,"-",id2,".png")
  png(file, res=600, width=3500, height=2200)
  par(mar=c(6,7,4,7))
  xlim = c(0,8)
  max = 100
  plot(xlim,c(-50,max),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i",main="")

  
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
    else if (as.numeric(id_data[i,4]) == 0) {
      text(i-0.25, 10, xpd=T, labels="n/a",col="#333333", cex=0.9)
    }
    
    text(i-0.7+0.25, -17, xpd=T, labels=mut[1],srt=90,adj=1, col="#333333", cex=1.1, font=font)
    text(i-0.3+0.25, -17, xpd=T, labels=mut[2],srt=90,adj=1, col="#333333", cex=0.9)
  }
  
  text(xmax/2, 150, xpd=T, labels=id, col="#333333", cex=1.1, font=2)
  text(xmax/2, 130, xpd=T, labels=id2, col="#333333", cex=1.1, font=1)
  
  segments(0,0,xmax,lwd=1.4,col="#b1b1b1", xpd=T)
  dev.off()
  system(paste("open", file))
} 

