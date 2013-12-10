######################################################################
# Plots correlations between HMGD, H1, MNase and expression
# at each nucleosome position
######################################################################

source("~/proj/script/R/lib/util.R")

library(RColorBrewer)
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



#
# Define nucleosome regions
# 
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

flanking <- 1000
pos <- seq(-flanking, flanking)



# add small pseudocount to avoid zeros when taking log
pseudo <- 0.5 * min(gene.tab["EXPR.S2_R."][gene.tab["EXPR.S2_R."] > 0])
expr.rpkm <- gene.tab[["EXPR.S2_R."]] + pseudo

f <- !is.na(expr.rpkm)

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

  pseudo <- 2
  
  region.h1.rpkm <- (region.h1.counts + pseudo) * 1e9 / (n.site * h1.ttl)

  region.hmgd.rpkm <- (region.hmgd.counts + 1) * 1e9 / (n.site * hmgd.ttl)

  region.mnase.rpkm <- (region.mnase.counts + 1) * 1e9 / (n.site * mnase.ttl)

  region.hmgd.h1.ratio <- log10(region.hmgd.rpkm / region.h1.rpkm)
  region.hmgd.mnase.ratio <- log10(region.hmgd.rpkm / region.mnase.rpkm)
  region.h1.mnase.ratio <- log10(region.h1.rpkm / region.mnase.rpkm)

  res <- cor.test(log10(expr.rpkm[f]), log10(region.mnase.rpkm[f]))
  region.tab[i, "MNASE.COR"] <- res$estimate
  region.tab[i,"MNASE.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"MNASE.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"MNASE.COR.P.VALUE"] <- res$p.value

  res <- cor.test(log10(expr.rpkm[f]), log10(region.hmgd.rpkm[f]))
  region.tab[i,"HMGD.COR"] <- res$estimate
  region.tab[i,"HMGD.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"HMGD.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"HMGD.COR.P.VALUE"] <- res$p.value

  res <- cor.test(log10(expr.rpkm[f]), log10(region.h1.rpkm[f]))
  region.tab[i,"H1.COR"] <- res$estimate
  region.tab[i,"H1.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"H1.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"H1.COR.P.VALUE"] <- res$p.value

  res <- cor.test(log10(expr.rpkm[f]), region.hmgd.h1.ratio[f])
  region.tab[i,"HMGD.H1.RATIO.COR"] <- res$estimate
  region.tab[i,"HMGD.H1.RATIO.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"HMGD.H1.RATIO.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"HMGD.H1.RATIO.P.VALUE"] <- res$p.value
  
  res <- cor.test(log10(expr.rpkm[f]), region.hmgd.mnase.ratio[f])
  region.tab[i,"HMGD.MNASE.RATIO.COR"] <- res$estimate
  region.tab[i,"HMGD.MNASE.RATIO.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"HMGD.MNASE.RATIO.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"HMGD.MNASE.RATIO.P.VALUE"] <- res$p.value

  res <- cor.test(log10(expr.rpkm[f]), region.h1.mnase.ratio[f])
  region.tab[i,"H1.MNASE.RATIO.COR"] <- res$estimate
  region.tab[i,"H1.MNASE.RATIO.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i,"H1.MNASE.RATIO.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i,"H1.MNASE.RATIO.P.VALUE"] <- res$p.value
  
  mdl <- lm(log10(expr.rpkm[f]) ~ log10(region.h1.rpkm[f]) + log10(region.hmgd.rpkm[f]))
  res <- cor.test(log10(expr.rpkm[f]), mdl$fitted.values)
  region.tab[i, "H1.HMGD.MDL.COR"] <- res$estimate
  region.tab[i, "H1.HMGD.MDL.COR.CI.LOW"] <- res$conf.int[1]
  region.tab[i, "H1.HMGD.MDL.COR.CI.UP"] <- res$conf.int[2]
  region.tab[i, "H1.HMGD.MDL.COR.P.VALUE"] <- res$p.value

  
  # is it significantly better to use both H1 and HMGD in linear
  # model than just H1 or just HMGD
  mdl.h1 <- lm(log10(expr.rpkm[f]) ~ log10(region.h1.rpkm[f]))
  mdl.hmgd <- lm(log10(expr.rpkm[f]) ~ log10(region.hmgd.rpkm[f]))
  mdl.hmgd.h1 <- lm(log10(expr.rpkm[f]) ~ log10(region.h1.rpkm[f]) +
                    log10(region.hmgd.rpkm[f]))
  
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


