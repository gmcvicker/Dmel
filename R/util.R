
#
# Adds a scale-bar to an image plot. Function by
# Jonathan Rougier (J.C.Rougier@durham.ac.uk), taken from:
# http://tolstoy.newcastle.edu.au/R/help/99b/0483.html
#
image.scale <- function (z, col, x, y = NULL, size = NULL, digits = 2,
                         labels = c("breaks", "ranges")) {
    # sort out the location
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05*dx/2 # default x to right of image
    else if (is.list(x)) {
        if (length(x$x) == 2)
          size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    } else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
          size <- 0.618*dy/n # default size, golden ratio
          y <- my + 0.618*dy/2 # default y to give centred scale
        } else size <- (y-my)*2/n
    if (length(size)==1)
        size <- rep(size, 2) # default square boxes
    if (is.null(y))
        y <- my + n*size[2]/2
    # draw the image scale
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
        col = rev(col), xpd = TRUE)
    # sort out the labels
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format="f", digits=digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj =
        ifelse(size[1]>0, 0, 1), xpd = TRUE)
}



#
# Places the Z data into bins that are divided linearly across the
# the range of X and Y
#
grid.data <- function(x, y, z, n.y.bin=10, n.x.bin=10,
                      bin.type="MEAN") {
  y.r <- range(y, na.rm=T)
  x.r <- range(x, na.rm=T)

  x.bin.sz <- (x.r[2] - x.r[1]) / n.x.bin
  
  y.bin.sz <- (y.r[2] - y.r[1]) / n.y.bin

  x.bin.start <- seq(0, n.x.bin-1, by=1) * x.bin.sz + x.r[1]
  x.bin.end   <- seq(1, n.x.bin, by=1) * x.bin.sz + x.r[1]
  y.bin.start <- seq(0, n.y.bin-1, by=1) * y.bin.sz + y.r[1]
  y.bin.end   <- seq(1, n.y.bin, by=1) * y.bin.sz + y.r[1]


  x.lab <- paste(signif(x.bin.start,2), signif(x.bin.end,2), sep="-")
  y.lab <- paste(signif(y.bin.start,2), signif(y.bin.end,2), sep="-")
  
  
  z.vals <- matrix(ncol=n.x.bin, nrow=n.y.bin, dimname=list(x.lab,y.lab))
  
  for(i in 1:n.x.bin) {
    x.filter <- x >= x.bin.start[i] & x < x.bin.end[i]
    for(j in 1:n.y.bin) {

      y.filter <- y >= y.bin.start[j] & y < y.bin.end[j]

      if(bin.type == "MEAN") {
        z.vals[i,j] <- mean(z[x.filter & y.filter])
      }
      else if(bin.type == "SUM") {
        z.vals[i,j] <- sum(z[x.filter & y.filter])
      }
      else {
        warning("unknown bin type, assuming MEAN")
        z.vals[i, j] <- mean(z[x.filter & y.filter])
      }        
    }
  }

  return(z.vals)
}


#
# Places data into bins with same number of datapoints in each
#
bin.data <- function(x, y=NA, n.bin=100) {
  ord <- order(x)
  x.ord <- x[ord]

  y.data <- length(y) > 1
  
  if(y.data) {
    y.ord <- y[ord]
  } else {
    y.ord <- NA
  }

  n.points <- length(x.ord)

  bin.sz <- n.points / n.bin

  bin.start <- 1

  bins <- list(x.start=rep(NA, n.bin),
               x.end=rep(NA, n.bin),
               x.mean=rep(NA, n.bin),
               y.mean=rep(NA, n.bin),
               y.sem=rep(NA, n.bin),
               y.sum=rep(NA, n.bin),
               n=rep(NA, n.bin))
                      
  for(i in 1:n.bin) {
    bin.end <- bin.start + bin.sz - 1
    if(bin.end > n.points) {
      bin.end <- n.points
    }
    
    bin.x <- x.ord[bin.start:bin.end]
    
    n <- length(bin.x)
    bins$n[i] <- n
    bins$x.mean[i] <- mean(bin.x)
    bins$x.start[i] <- bin.x[1]
    bins$x.end[i] <- bin.x[length(bin.x)]

    if(y.data) {    
      bin.y <- y.ord[bin.start:bin.end]
      bins$y.mean[i] <- mean(bin.y)
      bins$y.sem[i] <- sd(bin.y) / sqrt(n)
      bins$y.sum[i] <- sum(bin.y)
    }

    bin.start <- bin.end + 1
  }
  return(bins)
}



#
# Barplot function that is somewhat more flexible than Rs
#
brplt <- function(data.matrix, col="grey", labels=c(),
                  ylab="", xlab="", legend.text=c(),
                  main="", density=NULL, angle=NULL,
                  xhatch=FALSE, bty="n", las=2,
                  ylim=c(0, max(data.matrix)), mcex=1.0) {
  has.legend <- length(legend.text) > 0
  
  n.bars <- nrow(data.matrix)
  n.grps <- ncol(data.matrix)

  x.positions.mtrx <- matrix(nrow=n.bars, ncol=n.grps)

  
  bar.width <- 1

  clrs <- rep(col, n.bars/length(col))
  
  # calculate width of bar-group based on single bar width
  grp.width <- (n.bars+1) * bar.width

  # calculate start points of each group
  grp.starts <- seq(0,(n.grps-1)) * grp.width + bar.width
  grp.mids <- (grp.starts + grp.width/2 - bar.width/2)
  
  bh <- bar.width/2 # half bar width

  plt.width <- max(grp.starts) + grp.width
  
  plot(x=c(0), y=c(0), xlim=c(0, plt.width), ylim=ylim, las=las,
       bty=bty, main=main, type="n", xaxt="n", xlab=xlab, ylab=ylab)

  if(has.legend) {
    legend("topright", legend.text, fill=clrs, cex=mcex,
           inset=0.02, bty="n")
  }
  
  for(i in 1:n.grps) {
    for(j in 1:n.bars) {
      bar.mid <- grp.starts[i] + (j-1)*bar.width + bh
      x.positions.mtrx[j,i] <- bar.mid
      
      bar.x <- c(bar.mid - bh, # left
                 bar.mid - bh,
                 bar.mid + bh, # right
                 bar.mid + bh)
      bar.y <- c(0, data.matrix[j,i], data.matrix[j,i], 0)
      
      polygon(bar.x, bar.y, col=clrs[j])

      if(!is.null(density)) {
        polygon(bar.x, bar.y, col="black", density=density[j],
                angle=angle[j])

        if(xhatch[j]) {
          polygon(bar.x, bar.y, col="black", density=density[j],
                  angle=angle[j]+90)
        }
      }
      
    }
  }

  mtext(labels, side=1, at=grp.mids, cex=mcex) 

  return(x.positions.mtrx)  
}


#
# less-verbose wrapper of paste function
#
p <- function(..., sep="") {
  paste(..., sep=sep)
}



#
# less-verbose wrapper of length function
#
len <- function(...) {
  length(...)
}


count.lines <- function(path) {
  is.gzipped <- length(grep("gz$", path)) > 0

  if(is.gzipped) {
    cmd <- paste("gunzip -c ",  path, "| wc -l | awk '{print $1}'")
  } else {
    cmd <- paste("wc -l", path, "| awk '{print $1}'")
  }

  res <- system(cmd, intern=TRUE)

  return(as.numeric(res))
}

count.cols <- function(path) {
  # return number of columns on first line
  tab <- read.table(path, nrows=1)
  return(ncol(tab))
}

#
# convenience wrapper which reads a matrix in from a file
#
read.matrix <- function(path, max.val=NULL, nrow=NULL, ncol=NULL) {
  if(is.null(ncol)) {
    ncol <- count.cols(path)
  }
  
  colClasses <- rep("numeric", ncol)

  if(is.null(nrow)) {
    nrow <- count.lines(path)
  }
  
  m <- as.matrix(read.table(path, header=F, colClasses=colClasses,
                            nrows=nrow))
  
  if(!is.null(max.val)) {
    m[!is.na(m) & (m > max.val)] <- max.val
  }
  return(m)
}


read.big.table <- function(path, header=F, comment.char="",
                           ...) {
  # read the first 50 rows to guess datatypes
  tab50rows <- read.table(path, header=header, nrows=50)
  classes <- sapply(tab50rows, class)

  cat("  table has ", length(classes), " columns\n")
  
  # count number of lines in table using wc
  nrow <- count.lines(path)

  if(header) {
    nrow <- nrow - 1
  }
  
  cat("  table has ", nrow, " rows\n")
  
  # now read entire table
  tab <- read.table(path, header=header, colClasses=classes, nrows=nrow,
                    comment.char=comment.char, ...)
  
  return(tab)
}


#
# convenience wrapper which reads a numeric vector in from a file
#
read.vector <- function(path, nrow=NULL) {  
  if(is.null(nrow)) {
    nrow <- count.lines(path)
  }

  tab <- read.table(path, header=F, colClasses="numeric",
                    nrow=nrow)
  
  return(as.numeric(tab$V1))
}



#
# smooths provided values using a sliding window of specified size
# performs no smoothing if win.sz is <= 1
#
sliding.window <- function(vals, win.sz) {
  if(win.sz <= 1) {
    # just return a copy, do no smoothing
    return(c(vals))
  }

  mid <- round(win.sz / 2)

  start <- 1
  end <- length(vals) - mid + 1

  smoothed.vals <- rep(NA, length(vals))
  
  for(i in start:end) {
    smoothed.vals[i + mid - 1] <- mean(vals[i:(i+win.sz-1)])
  }

  return(smoothed.vals)
}




#
# Returns the central %-range of the
# provided values
#
central.range <- function(vec, pct) {
  ord <- order(vec)
  
  n <- length(vec)
  from <- ((n - (n * pct)) / 2) + 1
  to   <- (n - from) + 1

  return(range(vec[ord][from:to]))
}



#
# Performs N chisq tests on 1x2 tables
# where N is the length of the provided vectors.
# x[1] vs v[1] and x[2] vs v[2] is tested and so on.
# The p-values from the resultant tests are returned
#
pairwise.chisq.pvalues <- function(x, y) {
  pvals <- rep(1, length(x))
  
  for(i in seq(1:length(x))) {
    if(x[i] == 0 && y[i] == 0) {
      pvals[i] <- NaN;
    } else {
      cs <- chisq.test(c(x[i], y[i]));
      pvals[i] <- cs$p.value;
    }
  }

  return(pvals);
}

