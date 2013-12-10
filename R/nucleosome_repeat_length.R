
###################################################################
# Write the regions stratified into HMGD/H1 ratio
#

tab <- read.table("~/data/Dmel/DNase/hmgd_h1_by_dhs_dist.txt",
                  header=T)

out.dir <- "~/data/Dmel/nuc_repeat_len_by_region/"


get.n.mapped <- function(type) {
  mapped.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt",
                           header=F)
  return(as.numeric(mapped.tab$V2[mapped.tab$V1==type]))
}


h1.count <- tab$H1.COUNT
hmgd.count <- tab$HMGD.COUNT
mnase.count <- tab$MNASE.COUNT

h1.mapped <- get.n.mapped("h1")
hmgd.mapped <- get.n.mapped("hmgd")
mnase.mapped <- get.n.mapped("S2")

region.size <- tab$END - tab$START + 1

h1.rpkm <- ((h1.count+1) * 1e9) / (region.size * h1.mapped)
hmgd.rpkm <- ((hmgd.count+1) * 1e9) / (region.size * hmgd.mapped)
mnase.rpkm <- ((mnase.mapped+1) * 1e9) / (region.size * mnase.mapped)

  
f <- (h1.count > 50) & (hmgd.count > 50) & (mnase.count > 50)

n.quantile <- 5

hmgd.h1.ratio <- log2((hmgd.rpkm) / (h1.rpkm))

hist(hmgd.h1.ratio[f], col="grey50", breaks=100)


quants <- (0:n.quantile)/n.quantile
q <- quantile(hmgd.h1.ratio, quants)
# q <- c(-2, -1, 0, 1, 2, 4)

for(i in 1:n.quantile) {
  q.start <- q[i]
  q.end <- q[i+1]

  # q.name <- paste(quants[i], "-", quants[i+1], sep="")
  q.name <- i

  cat(q.name, "\n")
  filt <- f & (hmgd.h1.ratio >= q.start) & (hmgd.h1.ratio < q.end)
  out.file <- paste(out.dir, "regions_hmgd_h1_ratio_Q", q.name, ".txt",
                    sep="")

  write.table(tab[filt,], file=out.file, quote=F, col.names=T,
              row.names=FALSE)
}


# stratify by H1 count
q <- quantile(h1.count, quants)

for(i in 1:n.quantile) {
  q.start <- q[i]
  q.end <- q[i+1]

  # q.name <- paste(quants[i], "-", quants[i+1], sep="")
  q.name <- i
  cat(q.name, "\n")  
  filt <- f & (h1.count >= q.start) & (h1.count < q.end)

  out.file <- paste(out.dir, "regions_h1_Q", q.name, ".txt",
                    sep="")

  write.table(tab[filt,], file=out.file, quote=F, col.names=T,
              row.names=FALSE)
}


# stratify by HMGD count
q <- quantile(hmgd.count, quants)

for(i in 1:n.quantile) {
  q.start <- q[i]
  q.end <- q[i+1]

  # q.name <- paste(quants[i], "-", quants[i+1], sep="")
  q.name <- i
  cat(q.name, "\n")  
  filt <- f & (hmgd.count >= q.start) & (hmgd.count < q.end)

  out.file <- paste(out.dir, "regions_hmgd_Q", q.name, ".txt",
                    sep="")

  write.table(tab[filt,], file=out.file, quote=F, col.names=T,
              row.names=FALSE)
}






#####################################################################


library(RColorBrewer)

# color.pal <- brewer.pal(9, "Set1")
color.pal <- c(rgb(0,0,0),  # black
               rgb(58.0/255, 51.0/255, 255.0/255), # blue
               rgb(0, 128.0/255, 0), # green
               rgb(143.0/255, 31.0/255, 255.0/255), # purple
               rgb(1.0, 0, 0)) # red


n.quantile <- 5

max.dist <- 1500

counts.matrix <- matrix(nrow=max.dist+1, ncol=n.quantile)
rpkm.matrix <- matrix(nrow=max.dist+1, ncol=n.quantile)
smooth.matrix <- matrix(nrow=max.dist+1, ncol=n.quantile)

repeat.len.tab <- data.frame(hmgd=rep(NA, n.quantile),
                             h1=rep(NA, n.quantile),
                             hmgd_h1_ratio=rep(NA, n.quantile))

types <- c("hmgd", "h1", "hmgd_h1_ratio")
labels <- c("HMGD", "H1", "HMGD / H1 ratio")

for(j in 1:3) {
  type <- types[j]
  label <- labels[j]

  cat(type, "\n")
  
  smooth.win <- 20
  
  for(i in 1:n.quantile) {
    q.name <- i
    
    cat(q.name, "\n")
    
    file <- paste(out.dir, "nrl_", type, "_Q", q.name, ".txt",
                  sep="")
    
    tab <- read.table(file, header=T)
    
    counts.matrix[,i] <- tab$COUNT[1:(max.dist+1)]

    rpkm.matrix[,i] <- (counts.matrix[,i] * 1e9) / (tab$N.SITES[1:(max.dist+1)] * mnase.mapped)
    
    smooth.matrix[,i] <- filter(rpkm.matrix[,i],
                                rep(TRUE, smooth.win))/smooth.win    
  }



  filename <- paste("nrl_by_", type, ".pdf", sep="")
  pdf(filename, width=6, height=6)
  
  dist <- 0:max.dist
  plot(c(0), c(0), type="n",
    #   ylim=c(0.0005, 0.001), xlim=range(dist),
       ylim=c(8, 20), xlim=range(dist),
       ylab="nucleosome midpoints (RPKM)", las=1,
       xlab="distance from anchor midpoint")
  
  for(i in 1:5) {
    f <- !is.na(smooth.matrix[,i])
    points(dist[f], smooth.matrix[f,i], #/
 #          sum(as.numeric(smooth.matrix[f,i])),
           col=color.pal[i], type="l", lwd=2)
  }
  
  legend("topright", legend=paste(label, " quintile ", 1:5, sep=""),
         col=color.pal, lty=1, lwd=2)
  
  dev.off()
  
  
                             
  filename <- paste("nrl_", type, "_analysis.pdf", sep="")
  pdf(filename, width=8, height=8)
  
  par(mfrow=c(3,1))

  repeat.len <- c()
  
  for(i in 1:n.quantile) {
    # take counts starting at position 100
    counts <- counts.matrix[100:nrow(counts.matrix),i]

    # fit quadratic to counts
    x <- 1:length(counts)
    x2 <- x*x
    mdl <- lm(counts ~ x + x2)
    y <- mdl$coef[1] +  x*mdl$coef[2] + x*x*mdl$coef[3]
    
    # regress out quadratic
    plot(x, counts, xlab="position", ylab="midpoint count",
         main=paste("Quadratic fit - quintile ", i))
    lines(x, y, col="red")
    
    plot(x, mdl$resid, xlab="position",
         ylab="residual midpoint count",
         main=paste("Residuals - quintile ", i))
    
    # compute autocorrelation
    res <- acf(mdl$resid, lag.max=length(counts),
               main=paste("Autocorrelation - quintile ", i))
    
    f <- res$lag > 100 & res$lag < 250
    repeat.len[i] <- res$lag[f][res$acf[f] == max(res$acf[f])]
    
    segments(x0=repeat.len[i], x1=repeat.len[i], y0=-1, y1=1, col="red")  
  }

  dev.off()

  # record repeat lengths
  repeat.len.tab[type] <- repeat.len
}

pdf("nrl_summary.pdf", width=5, height=5)

plot(c(0), c(0), type="n",
     xlab="quintile", las=1, ylim=range(repeat.len.tab),
     xlim=c(1, n.quantile),
     ylab="nucleosome repeat length (bp)")

for(i in 1:ncol(repeat.len.tab)) {
  label <- labels[i]

  points(x=1:n.quantile, y=repeat.len.tab[,i],
         col=color.pal[i], lwd=2, type="o", pch=21)
}

legend("topright", legend=labels, 
       col=color.pal[1:ncol(repeat.len.tab)],
       pch=21, lwd=2)

dev.off()










## pdf("nrl_by_hmgd_h1_ratio.pdf", width=6, height=6)

## dist <- 0:max.dist
## plot(c(0), c(0), type="n",
##      xlim=c(0, n.peaks), ylim=range(peak.matrix),
##      ylab="distance from anchor midpoint",
##      xlab="nucleosome number")

## repeat.len <- c()

## for(i in 1:n.quantile) {
##   peak.num <- 1:n.peaks
##   peak.dist <- peak.matrix[,i]
  
##   points(1:n.peaks, peak.matrix[,i], pch=21, bg=color.pal[i], 
##          col=color.pal[i], type="p")

##   mdl <- lm(peak.dist ~ peak.num)

##   repeat.len[i] <- mdl$coef[2]

##   abline(mdl$coef, col=color.pal[i])
  
##   legend("bottomright", legend=paste("HMGD/H1 ratio quintile ", 1:5),
##          col=color.pal, lty=1)
## }

## dev.off()
