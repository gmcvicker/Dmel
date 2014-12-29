
library(RColorBrewer)

source("/mnt/lustre/home/gmcvicker/proj/script/R/lib/util.R")


draw.heatmap.legend <- function(neg.colors, pos.colors, scale) {

  n.neg.colors <- length(neg.colors)
  n.pos.colors <- length(pos.colors)
  
  colors <- c(rev(neg.colors), pos.colors)
  n.colors <- length(colors)

  # initalize plot, getting rid of axis labels, outer frame and most margins
  #old.mar <- par('mar')
  #old.oma <- par('oma')
  #par(mar=c(4, 2, 1, 4))
  #par(oma=c(0, 0, 0, 0))  
  plot(x=c(0), y=c(0), ylim=c(0, n.colors), xlim=c(0, 1),
       type="n", xaxt="n", yaxt="n", xlab="", ylab="", xpd=NA,
       frame.plot=FALSE)
  #par(mar=old.mar)
  #par(oma=old.oma)

  # draw filled rectangles of the heatmap colors
  ybottom <- seq(0, n.colors-1)
  ytop <- ybottom + 1
  xleft <- 0
  xright <- 1
  rect(xleft, ybottom, xright, ytop, col=colors,
       border=NA)

  # draw the axis on the right side
  axis.vals <- c(-rev(seq(1, n.neg.colors-1) * (scale / n.neg.colors)),
                 0,
                 seq(1, n.pos.colors-1) * (scale / n.pos.colors))
  
  axis(side=2, at=seq(1, n.colors-1),
       labels=sprintf("%.2f", axis.vals),
       outer=FALSE,
       las=1)
}


draw.heatmap <- function(mat,  xlab="", ylab="",
                         labels=NA, scale=NA) {
  library(RColorBrewer)

  m <- mat
  
  # choose matching number of positive / negative colors
  neg.colors <- brewer.pal(9, "Purples")[1:9]
  pos.colors <- brewer.pal(9, "Oranges")[1:9]
  # colors <- c(neg.colors, pos.colors)
  n.colors <- length(pos.colors)

  if(is.na(scale)) {
    scale <- max(abs(quantile(m, c(0.05, 0.95))))
  }
  min.score <- -scale
  max.score <- scale
  
  # determine color each element of matrix should be
  color.unit <- scale / n.colors
  color.idx <- round(abs(m)/color.unit) + 1
  color.idx[color.idx > n.colors] <- n.colors
  
  color.idx[color.idx < 1] <- 1
  color.idx[color.idx > n.colors] <- n.colors

  layout(matrix(c(1, 2), ncol=2),  widths=c(1, 0.2), heights=c(1))
  # layout.show(2)
  
  plot(x=c(0, ncol(m)), y=c(0, nrow(m)),
       type="n", xaxt="n", yaxt="n", ylab=ylab, xlab=xlab,
       frame.plot=FALSE)
  
  # determine vertical location of rows
  # separating clusters with some whitespace
  ybottom <- seq(0, nrow(m)-1)
  ytop <- ybottom + 1
  
  # draw main heatmap, one column at a time
  for(i in 1:ncol(m)) {
    xleft <- rep(i-1, nrow(m))
    xright <- xleft + 1
    
    vals <- m[,i]
    colors <- pos.colors[color.idx[,i]]
    colors[vals < 0] <- neg.colors[color.idx[vals < 0,i]]
    
    rect(xleft, ybottom, xright, ytop, col=colors,
         border=NA)
  }
  
  # draw labels for each phenotype
  mids <- 1:ncol(m) - 0.5
  if(is.na(labels)) {
    labels <- colnames(m)
  }
  # text(x=mids, y=1, labels=labels, pos=1)
  axis(side=1, at=mids, labels=labels, tick=FALSE, las=2)
  
  draw.heatmap.legend(neg.colors, pos.colors, scale)
}

