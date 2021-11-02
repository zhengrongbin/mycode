rect_plot = function(rose, adjRP, x1Label, x2Label){
  x1 = t(rose)
  x2 = t(adjRP)
  if (is.data.frame(x1)) {
    x1 = as.matrix(x1)
    x2 = as.matrix(x2)
  }
  xlim = c(0,nrow(x1))
  ylim = c(-1,ncol(x1)+1)
  attach(mtcars)  
  opar= par(no.readonly=TRUE) 
  par(mai=c(0.4,1,0.2,1),cex=0.6, xpd = T)
  plot(1:nrow(x1), 1:nrow(x1), main = "", xlim=xlim, ylim=ylim, type="n", xaxs = "r", yaxs="r", axes=F, xlab="", ylab="")
  axis(1, at=1:nrow(x1), labels=row.names(x1), las=2, tick=F, line=NA,outer=F, hadj=0)
  axis(2, at=1:ncol(x1), labels=colnames(x1), las=1, tick=F, line=NA,outer=F, hadj=0.5)
  Dat_x1 <- c(t(x1))
  rank_x1 <- rank(Dat_x1)
  # colors = sapply(Dat, function(x) if (x == 1) {return(NA)}else if (x==0){return("red")
  #   }else if (x == -1){return("blue")})
  colfunc <- colorRampPalette(c("blue", "red"))
  # colors_x1 = colorRamp2(breaks = c(min(Dat_x1),max(Dat_x1)), colors = c("blue","red"))
  colors_x1 = colfunc(length(Dat_x1))
  for (j in 1:nrow(x1)) {
    for (i in 1:ncol(x1)) {
      rect(j-0.5,i-0.5,j+0.5,i+0.5)
      polygon(c(j,j,j+1)-0.5,c(i-1,i,i)+0.5,col=colors_x1[rank_x1[i+ncol(x1)*(j-1)]], border="black")
    }
  }

  Dat_x2 <- c(t(x2))
  rank_x2 <- rank(Dat_x2)
  # colors_x2 = colorRamp2(breaks = c(min(Dat_x2),max(Dat_x2)), colors = c("blue","red"))
  colors_x2 = colfunc(length(Dat_x2))
  for (j in 1:nrow(x2)) {
    for (i in 1:ncol(x2)) {
      polygon(c(j,j+1,j+1)-0.5,c(i-1,i-1,i)+0.5,col=colors_x2[rank_x2[i+ncol(x2)*(j-1)]], border="white")
    }
  }
  # legend(0,73, legend = c(x1Label, x2Label),
  #        col = c('red', 'blue'), pch = 15, cex = 1.3, horiz = T)
  par(opar)  
  detach(mtcars)
}
