

library(RColorBrewer)
source("~/proj/script/R/lib/util.R")





color.pal <-  brewer.pal(9, "Set1")


data.dir <- "~/data/Dmel/MNase/tss_mnase_midpoints_by_gene"

read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")

gene.tab <- read.table(p(data.dir, "/gene_summary.fb_gene_names.txt.gz"), header=T)

h1.matrix <- read.matrix(p(data.dir, "/h1.txt.gz"))
h1.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]

hmgd.matrix <- read.matrix(p(data.dir, "/hmgd.txt.gz"))
hmgd.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "hmgd"]

mnase.matrix <- read.matrix(p(data.dir, "/S2.txt.gz"))
mnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]

dnase.matrix <- read.matrix(p(data.dir, "/dnase.txt.gz"))
dnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "dnase"]

max.count <- 20
mnase.matrix[mnase.matrix > max.count] <- max.count
hmgd.matrix[hmgd.matrix > max.count] <- max.count
h1.matrix[h1.matrix > max.count] <- max.count



flanking <- 1000


ds.repeat.len <- 167
us.repeat.len <- 177

nfr.start <- -117
nfr.end <- 59

us.ends <- rev(seq(nfr.start-1,
                   nfr.start-1-us.repeat.len*4, by=-us.repeat.len))
us.starts <- us.ends - us.repeat.len + 1

ds.starts <- seq(60, 60 + ds.repeat.len*4, by=ds.repeat.len)
ds.ends <- ds.starts + ds.repeat.len - 1

upstream.start <- min(us.starts)
upstream.end <- max(us.ends)

downstream.start <- min(ds.starts)
downstream.end <- max(ds.ends)

region.tab <- data.frame(NAME=c(paste("minus", rev(1:5), sep=""),
                           "NFR",
                           paste("plus", 1:5, sep=""),
                           "upstream", "downstream"),
                         START=c(us.starts, nfr.start,
                           ds.starts, upstream.start, downstream.start),
                         END=c(us.ends, nfr.end, ds.ends,
                           upstream.end, downstream.end))


pos <- seq(-flanking, flanking)


expr.pct <- gene.tab[["EXPR.PCT.S2_R."]]

expr.filters <-
  data.frame(high.expr = !is.na(expr.pct) & (expr.pct >= 0.75),
             mid.expr = !is.na(expr.pct) & (expr.pct >= 0.25) &
                         (expr.pct <= 0.75),
             low.expr = !is.na(expr.pct) & (expr.pct < 0.25))


# aggregate plot for all genes showing region break down
pdf("mnase_h1_hmgd_tss_regions.pdf", width=6, height=12)

par(mfrow=c(3,1))

ylim <- c(0, 50.0)

for(expr.level in names(expr.filters)) {
  f <- expr.filters[[expr.level]]

  mean.h1.mids <- apply(h1.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(h1.ttl))

  mean.hmgd.mids <- apply(hmgd.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(hmgd.ttl))

  mean.mnase.mids <- apply(mnase.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(mnase.ttl))

  mean.dnase.mids <- apply(dnase.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(dnase.ttl))
  
  plot(pos, mean.mnase.mids,
       type="l", las=1, main=expr.level,
       ylab="RPKM",
       ylim=ylim)

  points(pos, mean.hmgd.mids, type="l", col=color.pal[1])

  points(pos, mean.h1.mids, type="l", col=color.pal[2])

  points(pos, mean.dnase.mids, type="l", col=color.pal[4])
  
  segments(x0=region.tab$START, x1=region.tab$START,
           y0=ylim[1], y1=ylim[2], col="grey70")

  segments(x0=region.tab$END, x1=region.tab$END,
           y0=ylim[1], y1=ylim[2], col="grey70")

  legend("topleft", legend=c("MNase", "H1", "HMGD", "DNase"),
         lty=c(1,1,1,1),
         col=c("black", color.pal[2], color.pal[1], color.pal[4]),
         bg="white")

  mids <- region.tab$START + (region.tab$END - region.tab$START) / 2
  n.reg <- nrow(region.tab)
  
  text(x=mids[1:(n.reg-2)], y=c(2),
       labels=region.tab$NAME[1:(n.reg-2)], pos=1)  
}

dev.off()




# aggregate plot for all genes showing H1 and HMGD normalized by total nucleosomes
pdf("mnase_h1_hmgd_ratios_tss_regions.pdf", width=6, height=12)

par(mfrow=c(3,1))

ylim <- c(-2, 2)

for(expr.level in names(expr.filters)) {
  f <- expr.filters[[expr.level]]

  mean.h1.mids <- apply(h1.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(h1.ttl))

  mean.hmgd.mids <- apply(hmgd.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(hmgd.ttl))

  mean.mnase.mids <- apply(mnase.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(mnase.ttl))

  h1.ratio <- log2(mean.h1.mids / mean.mnase.mids)
  hmgd.ratio <- log2(mean.hmgd.mids / mean.mnase.mids)

  
  plot(pos, h1.ratio,
       type="l", las=1, main=expr.level,
       col=color.pal[2],
       ylab="H1 / total nucleosomes (log2 RPKM ratio)",
       ylim=ylim)

  points(pos, hmgd.ratio, type="l", col=color.pal[1])
  
  segments(x0=region.tab$START, x1=region.tab$START,
           y0=ylim[1], y1=ylim[2], col="grey70")

  segments(x0=region.tab$END, x1=region.tab$END,
           y0=ylim[1], y1=ylim[2], col="grey70")

  segments(x0=-1010, x1=1010, y0=0, y1=0, col="grey70", lty=2)

  legend("topleft", legend=c("H1", "HMGD"),
         lty=c(1,1),
         col=c(color.pal[2], color.pal[1]), bg="white")

  mids <- region.tab$START + (region.tab$END - region.tab$START) / 2
  n.reg <- nrow(region.tab)
  
  text(x=mids[1:(n.reg-2)], y=c(-1.8),
       labels=region.tab$NAME[1:(n.reg-2)], pos=1)  
}

dev.off()















#####
# plot by datatype (HMGD1, H1, HMGD1/MNase, H1/MNase) with different lines
# for each expression level
pdf("mnase_h1_hmgd_tss_by_expr.pdf", width=6, height=12)

par(mfrow=c(2,2))

clrs <- c(color.pal[1], color.pal[4], color.pal[2])
  
ylim <- c(0, 0.06)
data.matrix <- hmgd.matrix
data.ttl <- hmgd.ttl

plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     ylim=ylim, las=1,
     main="HMGD1",
     ylab="midpoints per base per million mapped reads",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  mean.mids <- apply(data.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(data.ttl))  
  points(pos, mean.mids, type="l", col=clrs[i])
}

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)

legend("topright", col=clrs, lty=1,
       legend=names(expr.filters))



ylim <- c(0, 0.02)
data.matrix <- h1.matrix
data.ttl <- h1.ttl

plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     las=1,
     ylim=ylim,
     main="H1",
     ylab="midpoints per base per million mapped reads",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  mean.mids <- apply(data.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(data.ttl))  
  points(pos, mean.mids, type="l", col=clrs[i])
}

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)



ylim <- c(-2, 2)

plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     las=1,
     ylim=ylim,
     main="HMGD1 normalized by total nucleosomes",
     ylab="log2(HMGD1 rate / total nucleosome rate)",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  hmgd.mids <- apply(hmgd.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(hmgd.ttl))
  mnase.mids <- apply(mnase.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(mnase.ttl))
  points(pos, log2(hmgd.mids / mnase.mids), type="l", col=clrs[i])
}

segments(x0=-1000, x1=1000,
         y0=0, y1=0, col="black", lty=2)

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)


ylim <- c(-2, 2)
plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     las=1,
     ylim=ylim,
     main="H1 normalized by total nucleosomes",
     ylab="log2(H1 rate / total nucleosome rate)",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  h1.mids <- apply(h1.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(h1.ttl))
  mnase.mids <- apply(mnase.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(mnase.ttl))
  points(pos, log2(h1.mids / mnase.mids), type="l", col=clrs[i])
}

segments(x0=-1000, x1=1000,
         y0=0, y1=0, col="black", lty=2)

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)

dev.off()














##################################################
# aggregate plots of high / low expression genes
##################################################



max.mismatch <- 1


tata.tab <- read.table(p("~/data/Dmel/tss_motif/tatabox.txt.gz"),
                       header=T)

tata.rev.tab <- read.table(p("~/data/Dmel/tss_motif/tatabox.rev.txt.gz"),
                       header=T)

ccaat.tab <- read.table(p("~/data/Dmel/tss_motif/ccaat.txt.gz"),
                        header=T)

mte.tab <- read.table(p("~/data/Dmel/tss_motif/mte.txt.gz"),
                      header=T)

bre.u.tab <- read.table(p("~/data/Dmel/tss_motif/bre.u.txt.gz"),
                       header=T)

bre.d.tab <- read.table(p("~/data/Dmel/tss_motif/bre.d.txt.gz"),
                       header=T)

inr.tab <- read.table(p("~/data/Dmel/tss_motif/inr.txt.gz"),
                       header=T)

dpe.tab <- read.table(p("~/data/Dmel/tss_motif/dpe.txt.gz"),
                       header=T)

max.mismatch <- 1
has.ccaat <- !is.na(ccaat.tab$MOTIF.MISMATCH) & ccaat.tab$MOTIF.MISMATCH <= max.mismatch

has.tata <- !is.na(tata.tab$MOTIF.MISMATCH) & tata.tab$MOTIF.MISMATCH <= max.mismatch

has.tata.rev <- !is.na(tata.rev.tab$MOTIF.MISMATCH) & tata.rev.tab$MOTIF.MISMATCH <= max.mismatch

has.bre.u <- !is.na(bre.u.tab$MOTIF.MISMATCH) & bre.u.tab$MOTIF.MISMATCH <= max.mismatch
has.bre.d <- !is.na(bre.d.tab$MOTIF.MISMATCH) & bre.d.tab$MOTIF.MISMATCH <= max.mismatch
has.inr <- !is.na(inr.tab$MOTIF.MISMATCH) & inr.tab$MOTIF.MISMATCH <= max.mismatch
has.dpe <- !is.na(dpe.tab$MOTIF.MISMATCH) & dpe.tab$MOTIF.MISMATCH <= max.mismatch


tata.promoter <- (has.tata)
dpe.promoter <- (has.dpe & has.inr) & !has.tata
other.promoter <- !has.tata & !has.inr

motif.filters <-
  data.frame(TATA.PROMOTER=tata.promoter,
             TATA.REV.PROMOTER=has.tata.rev,
             INR.DPE.PROMOTER=dpe.promoter,
             OTHER.PROMOTER=other.promoter)





xlim <- c(-750, 750)

## val.matrix <- h1.matrix
## val.ttl <- h1.ttl
## ylim <- c(0, 0.02)
## data.type <- "H1"

## val.matrix <- hmgd.matrix
## val.ttl <- hmgd.ttl
## ylim <- c(0, 0.06)
## data.type <- "HMGD"

val.matrix <- mnase.matrix
val.ttl <- mnase.ttl
ylim <- c(0, 0.02)
data.type <- "MNase"


sliding.window <- function(vals, win.sz=10) {
  window <- rep(1, win.sz)
  smoothed <- filter(vals, window) / sum(window)
  smoothed[is.na(smoothed)] <- 0
  return(smoothed)
}


for(expr.level in names(expr.filters)) {
  pdf(p(data.type, "_tss.by_expr.", expr.level, ".pdf"), width=7,
      height=7)

  par(mfrow=c(1, 1))

  expr.filter <- expr.filters[[expr.level]]
  
  all.mids <- apply(val.matrix[expr.filter,], 2, sum)*1e6 /
    (val.ttl * as.numeric(sum(expr.filter)))

  all.mids.smooth <- sliding.window(all.mids)

  plot(c(0), c(0),
       las=1,
       xlab="distance from TSS (bp)",
       ylab="midpoints per site / million mapped reads",
       xlim=xlim, type="n",
       main=data.type, ylim=ylim)
      
    segments(x0=region.tab$START, x1=region.tab$START,
             y0=0, y1=ylim[2], col="grey70")

    segments(x0=region.tab$END, x1=region.tab$END,
             y0=0, y1=ylim[2], col="grey70")

    lines(x=c(0,0), y=ylim, col="blue", lty=2, lwd=2)    
    points(pos, all.mids.smooth, type="l", col="black")

  ns <- rep(0, ncol(motif.filters))
  
  for(i in 1:ncol(motif.filters)) {
    
    motif <- names(motif.filters)[i]

    filter <- motif.filters[[motif]] & expr.filters[[expr.level]]

    mids <- apply(val.matrix[filter,], 2, sum)*1e6 /
      (val.ttl * as.numeric(sum(filter)))

    mids.smooth <- sliding.window(mids)

    points(pos, mids.smooth, type="l", col=color.pal[i])

    ns[i] <- sum(filter)
  }

  legend("topleft", legend=c(p("all (n=", sum(expr.filter), ")"),
                      p(names(motif.filters), " (n=", ns, ")")),
         col=c("black", color.pal[1:ncol(motif.filters)]),
           lty=c(1,1,1),  bg="white")

  
  dev.off()
}






expr <- gene.tab[["EXPR.S2_R."]]

f <- !is.na(expr)

for(i in 1:nrow(region.tab)) {
  region.name <- as.character(region.tab$NAME[i])
  region.start <- region.tab$START[i]
  region.end <- region.tab$END[i]

  cat(region.name, "\n")
  
  in.region <- (pos >= region.start) & (pos <= region.end)
  n.site <- as.numeric(sum(in.region))

  region.h1.counts <- apply(h1.matrix[,in.region], 1, sum)
  region.hmgd.counts <- apply(hmgd.matrix[,in.region], 1, sum)
  region.mnase.counts <- apply(mnase.matrix[,in.region], 1, sum)
      
  region.h1.rate <- (region.h1.counts + 1) * 1e6 /
    ((n.site) * h1.ttl)

  region.hmgd.rate <- (region.hmgd.counts + 1) * 1e6 /
    ((n.site) * hmgd.ttl)

  region.mnase.rate <- (region.mnase.counts + 1) * 1e6 /
    ((n.site) * mnase.ttl)

  region.hmgd.h1.ratio <- log2(region.hmgd.rate / region.h1.rate)
  region.hmgd.mnase.ratio <- log2(region.hmgd.rate / region.mnase.rate)
  region.h1.mnase.ratio <- log2(region.h1.rate / region.mnase.rate)

  res <- cor.test(expr[f], log2(region.mnase.rate[f]))
  region.tab[i, "MNASE.COR"] <- res$estimate
  region.tab[i,"MNASE.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"MNASE.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"MNASE.COR.P.VALUE"] <- res$p.value

  res <- cor.test(expr[f], log2(region.hmgd.rate[f]))
  region.tab[i,"HMGD.COR"] <- res$estimate
  region.tab[i,"HMGD.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"HMGD.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"HMGD.COR.P.VALUE"] <- res$p.value

  res <- cor.test(expr[f], log2(region.h1.rate[f]))
  region.tab[i,"H1.COR"] <- res$estimate
  region.tab[i,"H1.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"H1.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"H1.COR.P.VALUE"] <- res$p.value

  res <- cor.test(expr[f], region.hmgd.h1.ratio[f])
  region.tab[i,"HMGD.H1.RATIO.COR"] <- res$estimate
  region.tab[i,"HMGD.H1.RATIO.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"HMGD.H1.RATIO.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"HMGD.H1.RATIO.P.VALUE"] <- res$p.value
  
  res <- cor.test(expr[f], region.hmgd.mnase.ratio[f])
  region.tab[i,"HMGD.MNASE.RATIO.COR"] <- res$estimate
  region.tab[i,"HMGD.MNASE.RATIO.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"HMGD.MNASE.RATIO.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"HMGD.MNASE.RATIO.P.VALUE"] <- res$p.value

  res <- cor.test(expr[f], region.h1.mnase.ratio[f])
  region.tab[i,"H1.MNASE.RATIO.COR"] <- res$estimate
  region.tab[i,"H1.MNASE.RATIO.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"H1.MNASE.RATIO.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"H1.MNASE.RATIO.P.VALUE"] <- res$p.value
  
  mdl <- lm(expr[f] ~ log2(region.h1.rate[f]) + log2(region.hmgd.rate[f]))
  res <- cor.test(expr[f], mdl$fitted.values)
  region.tab[i, "H1.HMGD.MDL.COR"] <- res$estimate
  region.tab[i, "H1.HMGD.MDL.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i, "H1.HMGD.MDL.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i, "H1.HMGD.MDL.COR.P.VALUE"] <- res$p.value

  
  # is it significantly better to use both H1 and HMGD in linear
  # model than just H1 or just HMGD
  mdl.h1 <- lm(expr[f] ~ log2(region.h1.rate[f]))
  mdl.hmgd <- lm(expr[f] ~ log2(region.hmgd.rate[f]))
  mdl.hmgd.h1 <- lm(expr[f] ~ log2(region.h1.rate[f]) +
                    log2(region.hmgd.rate[f]))
  
  res1 <- anova(mdl.h1, mdl.hmgd.h1)
  res2 <- anova(mdl.hmgd, mdl.hmgd.h1)  
  region.tab[i, "H1.vs.BOTH.F.P.VAL"] <- res1[,6][2]
  region.tab[i, "HMGD.vs.BOTH.F.P.VAL"] <- res2[,6][2]
  
}




pdf("hmgd_h1_region_cors.pdf", width=7, height=8)

idx <- 1:(nrow(region.tab)-2)

ylim=range(region.tab[idx,4:ncol(region.tab)], -0.5, 0.75)

plot(idx, region.tab[idx,"MNASE.COR"], xaxt="n", type="o",
     ylim=ylim, las=1, cex=0.5, pch=21, xlab="region",
     ylab="pearson correlation", col="black", bg="black")

axis(1, at=idx, labels=region.tab$NAME[idx], cex.axis=0.75)

segments(x0=idx, x1=idx, y0=region.tab[idx,'MNASE.COR.CI.LOW'],
         y1=region.tab[idx,"MNASE.COR.CI.UP"])

lines(x=c(min(idx)-1,max(idx)+1), y=c(0,0), col="grey50", lty=2)

segments(x0=idx, x1=idx, y0=region.tab[idx,'HMGD.COR.CI.LOW'],
         y1=region.tab[idx,"HMGD.COR.CI.UP"], col=color.pal[1])
         
points(idx, region.tab[idx,"HMGD.COR"], type="o",
       col=color.pal[1], cex=0.5, pch=21, bg=color.pal[1])

segments(x0=idx, x1=idx, y0=region.tab[idx,'H1.COR.CI.LOW'],
         y1=region.tab[idx,"H1.COR.CI.UP"], col=color.pal[2])

points(idx, region.tab[idx,"H1.COR"], type="o",
       col=color.pal[2], cex=0.5, pch=21, bg=color.pal[2])

points(idx, region.tab[idx,"HMGD.H1.RATIO.COR"], type='o',
       col=color.pal[3], cex=0.5, pch=21, bg=color.pal[3])

segments(x0=idx, x1=idx, y0=region.tab[idx,'HMGD.H1.RATIO.COR.CI.LOW'],
         y1=region.tab[idx,"HMGD.H1.RATIO.COR.CI.UP"], col=color.pal[3])



legend("bottomleft", legend=c("MNase", "HMGD", "H1",
                       "HMGD/H1 ratio"),
       pt.cex=0.5, lty=1, col=c("black", color.pal[1:4]),
       pch=21, pt.bg=c("black", color.pal[1:4]))


dev.off()








#
# Scatter plots of H1 / HMGD ratio
#

library(ggplot2)



expr <- gene.tab[["EXPR.S2_R."]]
f <- !is.na(expr)

dnase.region.name <- "NFR"
dnase.region.start <- region.tab$START[region.tab$NAME == dnase.region.name]
dnase.region.end <- region.tab$END[region.tab$NAME == dnase.region.name]
in.region <- (pos >= dnase.region.start) & (pos <= dnase.region.end)
n.site <- as.numeric(sum(in.region))
region.dnase.counts <- apply(dnase.matrix[,in.region], 1, sum)



region.name <- "plus1"
region.start <- region.tab$START[region.tab$NAME == region.name]
region.end <- region.tab$END[region.tab$NAME == region.name]

in.region <- (pos >= region.start) & (pos <= region.end)
n.site <- as.numeric(sum(in.region))

region.hmgd.counts <- apply(hmgd.matrix[,in.region], 1, sum)
region.h1.counts <- apply(h1.matrix[,in.region], 1, sum)

ratio <- log2((region.hmgd.counts+1) / (region.h1.counts+1))

expr.level <- factor(rep("low", sum(f)), c("low", "medium", "high"))
expr.level[expr[f] <= -5] <- "low"
expr.level[expr[f] > -5 & expr[f] <= -2] <- "medium"
expr.level[expr[f] > -2] <- "high"

dat <- data.frame(dnase=log10(region.dnase.counts[f]+1),
                  expr=expr[f], hmgd.h1.ratio=ratio[f],
                  expr.level=expr.level)

library(ggplot2)

pdf(p("HMGD_H1_ratio_vs_expr.", region.name, ".pdf"), width=6, height=6)

ggplot(dat, aes(x=expr, y=hmgd.h1.ratio, color=expr.level)) +
  geom_point(shape=19,  alpha=1/10) +
  geom_smooth(method=loess)

dev.off()


ggplot(dat, aes(x=dnase, y=hmgd.h1.ratio, color=expr.level)) +
  geom_point(shape=19,  alpha=1/3)



mdl <- lm(dat$hmgd.h1.ratio ~ dat$expr)
summary(mdl)

mdl <- lm(dat$hmgd.h1.ratio ~ dat$dnase)
summary(mdl)

mdl <- lm(dat$hmgd.h1.ratio ~ dat$expr + dat$dnase)
summary(mdl)

ggplot(dat, aes(x=dnase, y=hmgd.h1.ratio)) +
  geom_point(shape=19,  alpha=1/10) +
  geom_smooth(method=loess)



# ggplot(dat, aes(expr)) + geom_density()

cor.test(dat$expr[expr.level == "low"],
         dat$hmgd.h1.ratio[expr.level == "low"])

cor.test(dat$expr[expr.level == "high"],
         dat$hmgd.h1.ratio[expr.level == "high"])

cor.test(dat$expr,
         dat$hmgd.h1.ratio)
