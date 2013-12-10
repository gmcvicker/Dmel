

source("~/proj/script/R/lib/util.R")

data.dir <- "~/data/Dmel/ss_mnase_mids/hmgd_h1"

read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")

ss.info.tab <- read.table(p(data.dir, "/info.txt.gz"), header=T)

h1.matrix <- read.matrix(p(data.dir, "/h1.txt.gz"))
hmgd.matrix <- read.matrix(p(data.dir, "/hmgd.txt.gz"))
mnase.matrix <- read.matrix(p(data.dir, "/mnase.txt.gz"))


threshold <- 5
h1.matrix[h1.matrix > threshold] <- threshold
hmgd.matrix[hmgd.matrix > threshold] <- threshold
mnase.matrix[mnase.matrix > threshold] <- threshold

h1.ttl.mapped <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]
hmgd.ttl.mapped <- read.totals.tab$V2[read.totals.tab$V1 == "hmgd"]
mnase.ttl.mapped <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]


flanking <- 200

n.h1 <- apply(h1.matrix, 1, sum)

n.site <- flanking + flanking + 1

filter <- (ss.info.tab$TRANSCRIPT.EXPR > 0) & (ss.info.tab$NEAR.TSS == 0)

hist(log10(n.h1[filter]+1) / n.site, col="grey80", breaks=50,
     xlab="log10(number of H11 midpoints per base)")
q <- log10(quantile(n.site, 0.99)/n.site.0)

lines(x=c(q, q), y=c(0, 50000), col="red")

# filter <- (ss.info.tab$TRANSCRIPT.EXPR > 0)


f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p"  
h1.5p <- apply(h1.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(h1.ttl.mapped))
f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p"  
h1.3p <- apply(h1.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(h1.ttl.mapped))
f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p" & ss.info.tab$IS.CANONICAL
h1.canon.5p <- apply(h1.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(h1.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p" & ss.info.tab$IS.CANONICAL
h1.canon.3p <- apply(h1.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(h1.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p" & !ss.info.tab$IS.CANONICAL
h1.noncanon.5p <- apply(h1.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(h1.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p" & !ss.info.tab$IS.CANONICAL
h1.noncanon.3p <- apply(h1.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(h1.ttl.mapped))



f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p"  
hmgd.5p <- apply(hmgd.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(hmgd.ttl.mapped))
f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p"  
hmgd.3p <- apply(hmgd.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(hmgd.ttl.mapped))
f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p" & ss.info.tab$IS.CANONICAL
hmgd.canon.5p <- apply(hmgd.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(hmgd.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p" & ss.info.tab$IS.CANONICAL
hmgd.canon.3p <- apply(hmgd.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(hmgd.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p" & !ss.info.tab$IS.CANONICAL
hmgd.noncanon.5p <- apply(hmgd.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(hmgd.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p" & !ss.info.tab$IS.CANONICAL
hmgd.noncanon.3p <- apply(hmgd.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(hmgd.ttl.mapped))



f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p"  
mnase.5p <- apply(mnase.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(mnase.ttl.mapped))
f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p"  
mnase.3p <- apply(mnase.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(mnase.ttl.mapped))
f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p" & ss.info.tab$IS.CANONICAL
mnase.canon.5p <- apply(mnase.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(mnase.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p" & ss.info.tab$IS.CANONICAL
mnase.canon.3p <- apply(mnase.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(mnase.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "5p" & !ss.info.tab$IS.CANONICAL
mnase.noncanon.5p <- apply(mnase.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(mnase.ttl.mapped))

f <- filter & ss.info.tab$SPLICE.SITE.TYPE == "3p" & !ss.info.tab$IS.CANONICAL
mnase.noncanon.3p <- apply(mnase.matrix[f,], 2, sum)*1e9 / (sum(f) * as.numeric(mnase.ttl.mapped))





win.sz <- 20
smooth.h1.5p <- filter(h1.5p, win.sz) / win.sz
smooth.h1.3p <- filter(h1.3p, win.sz) / win.sz
smooth.h1.canon.5p <- filter(h1.canon.5p, win.sz) / win.sz
smooth.h1.canon.3p <- filter(h1.canon.3p, win.sz) / win.sz
smooth.h1.noncanon.5p <- filter(h1.noncanon.5p, win.sz) / win.sz
smooth.h1.noncanon.3p <- filter(h1.noncanon.3p, win.sz) / win.sz

smooth.hmgd.5p <- filter(hmgd.5p, win.sz) / win.sz
smooth.hmgd.3p <- filter(hmgd.3p, win.sz) / win.sz
smooth.hmgd.canon.5p <- filter(hmgd.canon.5p, win.sz) / win.sz
smooth.hmgd.canon.3p <- filter(hmgd.canon.3p, win.sz) / win.sz
smooth.hmgd.noncanon.5p <- filter(hmgd.noncanon.5p, win.sz) / win.sz
smooth.hmgd.noncanon.3p <- filter(hmgd.noncanon.3p, win.sz) / win.sz

smooth.mnase.5p <- filter(mnase.5p, win.sz) / win.sz
smooth.mnase.3p <- filter(mnase.3p, win.sz) / win.sz
smooth.mnase.canon.5p <- filter(mnase.canon.5p, win.sz) / win.sz
smooth.mnase.canon.3p <- filter(mnase.canon.3p, win.sz) / win.sz
smooth.mnase.noncanon.5p <- filter(mnase.noncanon.5p, win.sz) / win.sz
smooth.mnase.noncanon.3p <- filter(mnase.noncanon.3p, win.sz) / win.sz


pos <- seq(-200, 200)

pdf("h1_hmgd_SS_aggregate.canon.noncanon.pdf", width=7.5, height=10)

par(mfrow=c(3,2))

ylim <- range(smooth.h1.3p, smooth.h1.5p, 0)

plot(pos, smooth.h1.canon.5p, type="l", col="blue", las=1,
     ylim=ylim, main="5' splice site: H1",
     ylab="midpoints per billion mapped reads",
     xlab="distance from 5' splice site (bp)")

points(pos, smooth.h1.noncanon.5p, type="l", col="red")

legend("topright", legend=c("canonical splice sites",
                     "non-canonical splice sites"),
       lty=1, col=c("blue", "red"))


plot(pos, smooth.h1.canon.3p, type="l", col="blue", las=1,
     ylim=ylim,
     ylab="midpoints per billion mapped reads",
     xlab="distance from 3' splice site (bp)",
     main="3' splice site: H1")

points(pos, smooth.h1.noncanon.3p, type="l", col="red")


ylim <- range(smooth.hmgd.3p, smooth.hmgd.5p, 0)

plot(pos, smooth.hmgd.canon.5p, type="l", col="blue", las=1,
     ylim=ylim, main="5' splice site: HMGD",
     ylab="midpoints per billion mapped reads",
     xlab="distance from 5' splice site (bp)")

points(pos, smooth.hmgd.noncanon.5p, type="l", col="red")

legend("topright", legend=c("canonical splice sites",
                     "non-canonical splice sites"),
       lty=1, col=c("blue", "red"))


plot(pos, smooth.hmgd.canon.3p, type="l", col="blue", las=1,
     ylim=ylim,
     ylab="midpoints per billion mapped reads",
     xlab="distance from 3' splice site (bp)",
     main="3' splice site: HMGD")

points(pos, smooth.hmgd.noncanon.3p, type="l", col="red")


ylim <- range(smooth.mnase.3p, smooth.mnase.5p, 0)

plot(pos, smooth.mnase.canon.5p, type="l", col="blue", las=1,
     ylim=ylim, main="5' splice site: MNase",
     ylab="midpoints per billion mapped reads",
     xlab="distance from 5' splice site (bp)")

points(pos, smooth.mnase.noncanon.5p, type="l", col="red")

legend("topright", legend=c("canonical splice sites",
                     "non-canonical splice sites"),
       lty=1, col=c("blue", "red"))


plot(pos, smooth.mnase.canon.3p, type="l", col="blue", las=1,
     ylim=ylim,
     ylab="midpoints per billion mapped reads",
     xlab="distance from 3' splice site (bp)",
     main="3' splice site: MNase")

points(pos, smooth.mnase.noncanon.3p, type="l", col="red")

dev.off()





pdf("h1_hmgd_SS_aggregate.normalized.pdf", width=7.5, height=7.5)
par(mfrow=c(2,2))

ylim <- range(log2(smooth.h1.3p / smooth.mnase.3p),
              log2(smooth.h1.5p / smooth.mnase.5p), 0)

h1.ratio.5p <- log2(smooth.h1.5p / smooth.mnase.5p)
h1.ratio.3p <- log2(smooth.h1.3p / smooth.mnase.3p)

plot(pos, h1.ratio.5p, type="l", col="blue", las=1,
     ylim=ylim, main="5' splice site: H1/MNase ratio",
     ylab="log2 ratio",
     xlab="distance from 5' splice site (bp)")

plot(pos, h1.ratio.3p, type="l", col="blue", las=1,
     ylim=ylim,
     ylab="log2 ratio",
     xlab="distance from 3' splice site (bp)",
     main="3' splice site: H1/MNase ratio")



ylim <- range(log2(smooth.hmgd.3p / smooth.mnase.3p),
              log2(smooth.hmgd.5p / smooth.mnase.5p), 0)

hmgd.ratio.5p <- log2(smooth.hmgd.5p / smooth.mnase.5p)
hmgd.ratio.3p <- log2(smooth.hmgd.3p / smooth.mnase.3p)

plot(pos, hmgd.ratio.5p, type="l", col="blue", las=1,
     ylim=ylim, main="5' splice site: HMGD/MNase ratio",
     ylab="log2 ratio",
     xlab="distance from 5' splice site (bp)")

plot(pos, hmgd.ratio.3p, type="l", col="blue", las=1,
     ylim=ylim,
     ylab="log2 ratio",
     xlab="distance from 3' splice site (bp)",
     main="3' splice site: HMGD/MNase ratio")


dev.off()










