
library(RColorBrewer)

color.pal <- brewer.pal(9, "Set1")

p <- function(...) {
  paste(..., sep="")
}

draw.guidelines <- function(positions, draw.at.dyad=TRUE,
                            draw.at.10=TRUE, draw.at.5=FALSE,
                            draw.at.core=TRUE, draw.at.linker=FALSE,
                            y.min=-20, y.max=20) {

  if(draw.at.10) {
    x <- c(seq(from=-10, to=min(positions), by=-10),
           seq(from=10, to=max(positions), by=10))
    segments(x0=x, y0=rep(y.min, length(x)),
             x1=x, y1=rep(y.max, length(x)),
             col="grey60", lty=2)
  }

  if(draw.at.5) {
    x <- c(seq(from=5, to=min(positions), by=-10),
           seq(from=-5, to=max(positions), by=10))
    segments(x0=x, y0=rep(y.min, length(x)),
             x1=x, y1=rep(y.max, length(x)),
             col="grey60", lty=1)
  }

  if(draw.at.dyad) {
    dyad <- c(0)
    segments(x0=dyad, y0=rep(y.min, length(dyad)),
             x1=dyad, y1=rep(y.max, length(dyad)),
             col="black", lty=1, lwd=1)
  }

  if(draw.at.core) {
    # draw lines at edge of "nucleosome core"
    core.edge <- c(-73,73)
    segments(x0=core.edge, y0=rep(y.min, length(dyad)),
             x1=core.edge, y1=rep(y.max, length(dyad)),
             col="black", lty=1, lwd=1)
  }

  if(draw.at.linker) {
    # draw lines at edge of "chromatosome" where linker histone protects DNA
    core.edge <- c(-88,88)
    segments(x0=core.edge, y0=rep(y.min, length(dyad)),
             x1=core.edge, y1=rep(y.max, length(dyad)),
             col="black", lty=1, lwd=1)
  }
}


  
draw.ww.ss.profile <- function(tab, xlim=c(-100, 100),
                               ylim=NULL, ww.col="#E41A1C",
                               ss.col="#377EB8", lwd=1, main="") {
    pos <- tab$POS
    
    n.sites <- tab$A + tab$C + tab$G + tab$T

    ww.dinuc <- c("AA", "AT", "TA", "TT")
    ss.dinuc <- c("GG", "GC", "CG", "CC")
    
    pseudo.count <- 1
    
    ylab <- ""
      
    # use proportion
    ww.vals <- apply(tab[,ww.dinuc], 1, sum) / n.sites
    ss.vals <- apply(tab[,ss.dinuc], 1, sum) / n.sites
    
    ylab <- "dinucleotide proportion"
    
    ylim.cur <- ylim
    if(is.null(ylim.cur)) {
        f <- (tab$POS > -65) & (tab$POS < 65)
        ylim.cur <- range(ww.vals[f], ss.vals[f])
    }
     
    plot(c(0), c(0), ylim=ylim.cur, type="n",
         xlab="distance from dyad (bp)", xlim=xlim, lwd=lwd, las=1,
         ylab=ylab, main=main)
      
    lines(tab$POS + 0.5, ww.vals, col=ww.col, lwd=lwd)
    lines(tab$POS + 0.5, ss.vals, col=ss.col, lwd=lwd)
      
    if((xlim[2] - xlim[1]) < 500) {
          draw.guidelines(pos, y.min=ylim.cur[1], y.max=ylim.cur[2])
    }
}    



xlim <- c(-80, 80)

cell.types <- c("antenna", "eye", "haltere", "leg")

for(cell.type in cell.types) {
  filename <-  p("~/data/Dmel/", cell.type, "_profile.txt")
  tab <- read.table(filename, header=T)
  tab <- tab[tab$POS >= xlim[1] & tab$POS <= xlim[2],]
  
  filename <- p("~/data/Dmel/", cell.type, "_dinuc_profile.pdf")
  pdf(filename, width=6, height=6)
  draw.ww.ss.profile(tab, xlim=xlim, main=cell.type, ylim=c(0.2, 0.32))
  dev.off()
}
