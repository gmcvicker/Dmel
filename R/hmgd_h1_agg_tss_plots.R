######################################################################
# Draws aggregate HMGD, H1, Total nucleosomes at TSSs
# stratified by expression level
######################################################################


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



expr.pct <- gene.tab[["EXPR.PCT.S2_R."]]

expr.filters <-
  data.frame(high.expr = !is.na(expr.pct) & (expr.pct >= 0.75),
             mid.expr = !is.na(expr.pct) & (expr.pct >= 0.25) &
                         (expr.pct <= 0.75),
             low.expr = !is.na(expr.pct) & (expr.pct < 0.25))



#####
# plot by datatype (HMGD1, H1, HMGD1/MNase, H1/MNase) with different lines
# for each expression level

pdf("mnase_h1_hmgd_tss_by_expr.pdf", width=6, height=12)

par(mfrow=c(3,1))


#
# Total nucleosomes
#
ylim <- c(0, 25)
data.matrix <- mnase.matrix
data.ttl <- mnase.ttl

plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     las=1,
     ylim=ylim,
     main="Total nucleosomes",
     ylab="midpoint density (RPKM)",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  mean.rpkm <- apply(data.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(data.ttl))  
  points(pos, mean.rpkm, type="l", col=clrs[i])
}

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)



#
# HMGD1 Plot
#

clrs <- c(color.pal[1], color.pal[4], color.pal[2])
  
ylim <- c(0, 60)
data.matrix <- hmgd.matrix
data.ttl <- hmgd.ttl

plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     ylim=ylim, las=1,
     main="HMGD1",
     ylab="midpoint density (RPKM)",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  mean.rpkm <- apply(data.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(data.ttl))  
  points(pos, mean.rpkm, type="l", col=clrs[i])
}

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)

legend("topright", col=clrs, lty=1,
       legend=names(expr.filters))


#
# H1 Plot
#
ylim <- c(0, 25)
data.matrix <- h1.matrix
data.ttl <- h1.ttl

plot(c(0), c(0), type="n",
     xlim=c(-1000, 1000),
     las=1,
     ylim=ylim,
     main="H1",
     ylab="midpoint density (RPKM)",
     xlab="distance from TSS (bp)")

for(i in 1:ncol(expr.filters)) {
  f <- expr.filters[,i]
  mean.rpkm <- apply(data.matrix[f,], 2, sum)*1e9 /
    (as.numeric(sum(f))*as.numeric(data.ttl))  
  points(pos, mean.rpkm, type="l", col=clrs[i])
}

segments(x0=0, x1=0, y0=ylim[1], y1=ylim[2],
         col="black", lty=2)



dev.off()
