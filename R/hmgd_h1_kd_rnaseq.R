


#
# Read RNA-seq expression
#
tab <- read.table("~/data/Dmel/KD_RNA_seq/tr_rna_seq_counts.txt", header=T)
expr <- tab[,c("RNASEQ.RPKM.S2_KD.Wt", "RNASEQ.RPKM.S2_KD.HMGD1KD", "RNASEQ.RPKM.S2_KD.H1KD")]

#
# read total number of mapped HMGD1, H1, and total nucleosome reads
#
read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")
h1.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]
hmgd.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "hmgd"]
mnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]

#
# Read number of mapped HMGD1, H1, and total nucleosome reads at each promoter
#
hmgd.h1.tab <- read.table("~/data/Dmel/KD_RNA_seq/hmgd_h1_at_rnaseq_promoters.txt", header=T)
tss.region.len <- hmgd.h1.tab$END - hmgd.h1.tab$START + 1
hmgd.rpkm <- ((hmgd.h1.tab$HMGD.COUNT+1) * 1e9) / (tss.region.len * hmgd.ttl)
h1.rpkm <- ((hmgd.h1.tab$H1.COUNT+1) * 1e9) / (tss.region.len * h1.ttl)
mnase.rpkm <- ((hmgd.h1.tab$MNASE.COUNT+1) * 1e9) / (tss.region.len * mnase.ttl)

hmgd.h1.ratio <- log2(hmgd.rpkm / h1.rpkm)


s2.expr <- expr[["RNASEQ.RPKM.S2_KD.Wt"]]
h1.kd.expr <- expr[["RNASEQ.RPKM.S2_KD.H1KD"]]
hmgd.kd.expr <- expr[["RNASEQ.RPKM.S2_KD.HMGD1KD"]]

# require above a minimum expression threshold in at least one experiment
min.expr <- 5

pdf("hmgd_h1_kd_expr_histograms.pdf", width=15, height=5)
par(mfrow=c(1,3))

hist(log2(s2.expr), breaks=50, main="S2 expr")
segments(x0=log2(min.expr), x1=log2(min.expr), y0=0, y1=2500, col="red")
         
hist(log2(hmgd.kd.expr), breaks=50, main="HMGD KD expr",
     xlab="expression (log2 RPKM)")
segments(x0=log2(min.expr), x1=log2(min.expr), y0=0, y1=2500, col="red")

hist(log2(h1.kd.expr), breaks=50, main="H1 KD expr",
     xlab="expression (log2 RPKM)")
segments(x0=log2(min.expr), x1=log2(min.expr), y0=0, y1=2500, col="red")

dev.off()

is.expressed <- !is.na(s2.expr) & !is.na(hmgd.kd.expr) & !is.na(h1.kd.expr) & 
                ((s2.expr > min.expr) | (hmgd.kd.expr > min.expr) | (h1.kd.expr > min.expr)) &
                tab[["RNASEQ.COUNT.S2_KD.Wt"]] > 25 &
                tab[["RNASEQ.COUNT.S2_KD.H1KD"]] > 25 &
                tab[["RNASEQ.COUNT.S2_KD.HMGD1KD"]] > 25
                

# there are a handful of outlier regions (probably collapsed repeats) with far more read counts
# than all other regions
filter <- hmgd.rpkm < 500 & h1.rpkm < 500 & mnase.rpkm < 500 &
          !is.na(hmgd.rpkm) & !is.na(h1.rpkm) & !is.na(mnase.rpkm) & is.expressed



boxplot(h1.rpkm[filter], hmgd.rpkm[filter], mnase.rpkm[filter])



test.expr.diff <- function(tab) {
  hmgd.pvals <- rep(NA, nrow(tab))
  h1.pvals <- rep(NA, nrow(tab))
                    
  for(i in 1:nrow(tab)) {
    wt.count <- tab[i, "RNASEQ.COUNT.S2_KD.Wt"]
    wt.total <- tab[i, "RNASEQ.TOTAL.MAPPED.S2_KD.Wt"]
    
    h1.count <- tab[i, "RNASEQ.COUNT.S2_KD.H1KD"]
    h1.total <- tab[i, "RNASEQ.TOTAL.MAPPED.S2_KD.H1KD"]

    hmgd.count <- tab[i, "RNASEQ.COUNT.S2_KD.HMGD1KD"]
    hmgd.total <- tab[i, "RNASEQ.TOTAL.MAPPED.S2_KD.HMGD1KD"]

    if(is.na(hmgd.count) | is.na(h1.count)) {
      h1.pvals[i] <- NA
      hmgd.pvals[i] <- NA
    } else {
    
      m <- matrix(c(wt.count, wt.total, h1.count, h1.total), nrow=2)
      res <- chisq.test(m)
      h1.pvals[i] <- res$p.val

      m <- matrix(c(wt.count, wt.total, hmgd.count, hmgd.total), nrow=2)
      res <- chisq.test(m)
      hmgd.pvals[i] <- res$p.val
    }
  }

  return(cbind(hmgd.pvals, h1.pvals))
}



hmgd.idx <- which(tab$GENE.SYMBOL == "HmgD")
h1.idx <- grep("^His1", tab$GENE.SYMBOL)

#
# Identify genes with change in expression in response to H1 or HMGD1 knockdown
#
hmgd.mdl <- lm(log2(hmgd.kd.expr[filter]) ~ log2(s2.expr[filter]))
hmgd.resid <- hmgd.mdl$residuals

h1.mdl <- lm(log2(h1.kd.expr[filter]) ~ log2(s2.expr[filter]))
h1.resid <- h1.mdl$residuals

h1.expr.change <- log2(h1.kd.expr) - log2(s2.expr)
hmgd.expr.change <- log2(hmgd.kd.expr) - log2(s2.expr)

expr.change.p.vals <- test.expr.diff(tab)
h1.p.vals <- expr.change.p.vals[,2]
hmgd.p.vals <- expr.change.p.vals[,1]

thresh <- 0.75
h1.upreg <- (h1.resid > thresh) & (h1.p.vals[filter] < 0.01)
h1.downreg <- (h1.resid < -thresh) & (h1.p.vals[filter] < 0.01)
h1.nochange <- !h1.upreg & !h1.downreg


pdf("h1_kd_expr_plots.pdf", width=10, height=5)

par(mfrow=c(1,2))

plot(log2(s2.expr[filter][h1.nochange]), log2(h1.kd.expr[filter][h1.nochange]), pch=20,
     col="grey70", xlim=range(log2(s2.expr[filter]), na.rm=T),
     ylab="H1 KD expr (log2 RPKM)", xlab="Wt expr (log2 RPKM)")
     ylim=range(log2(h1.kd.expr[filter]), log2(hmgd.kd.expr[filter]), na.rm=T))

points(log2(s2.expr[filter][h1.upreg]), log2(h1.kd.expr[filter][h1.upreg]),
       col="orange", pch=20)
abline(h1.mdl$coef, col="red")
points(log2(s2.expr[filter][h1.downreg]), log2(h1.kd.expr[filter][h1.downreg]),
       col="blue", pch=20)

points(log2(s2.expr[h1.idx]), log2(h1.kd.expr[h1.idx]), col="green",
       cex=2.0, pch=21)


#
# compare HMGD1/H1 ratio in genes that changed vs. did not change in expression
#
boxplot(hmgd.h1.ratio[filter][h1.upreg],
        hmgd.h1.ratio[filter][h1.nochange],
        hmgd.h1.ratio[filter][h1.downreg],
        col=c("orange", "grey70", "blue"),
        ylab="HMGD1/H1 ratio (log2 RPKM ratio)",
        las=1, main="H1 knockdown", outline=F)

dev.off()



pdf("hmgd_kd_expr_plots.pdf", width=10, height=5)

par(mfrow=c(1,2))

thresh <- 0.75
hmgd.upreg <- (hmgd.resid > thresh) & (hmgd.p.vals[filter] < 0.01)
hmgd.downreg <- (hmgd.resid < -thresh) & (hmgd.p.vals[filter] < 0.01)
hmgd.nochange <- !hmgd.upreg & !hmgd.downreg & abs(hmgd.resid) < thresh

plot(log2(s2.expr[filter][hmgd.nochange]), log2(hmgd.kd.expr[filter][hmgd.nochange]), pch=20,
     col="grey70", xlim=range(log2(s2.expr[filter]), na.rm=T),
     ylim=range(log2(h1.kd.expr[filter]), log2(hmgd.kd.expr[filter]), na.rm=T),
     ylab="HMGD KD expr (log2 RPKM)", xlab="Wt expr (log2 RPKM)")

points(log2(s2.expr[filter][hmgd.upreg]), log2(hmgd.kd.expr[filter][hmgd.upreg]),
       col="orange", pch=20)
abline(hmgd.mdl$coef, col="red")
points(log2(s2.expr[filter][hmgd.downreg]), log2(hmgd.kd.expr[filter][hmgd.downreg]),
       col="blue", pch=20)

points(log2(s2.expr[hmgd.idx]), log2(hmgd.kd.expr[hmgd.idx]), col="green",
       pch=21, cex=2.0)

boxplot(hmgd.h1.ratio[filter][hmgd.upreg],
        hmgd.h1.ratio[filter][hmgd.nochange],
        hmgd.h1.ratio[filter][hmgd.downreg],
        col=c("orange", "grey70", "blue"),
        ylab="HMGD1/H1 ratio (log2 RPKM ratio)",
        las=1)

dev.off()











################################

source("~/proj/script/R/lib/util.R")


pdf("hmgd_h1_kd_expression_change.pdf", width=6, height=6)

n.bin <- 5

h1.bins <- bin.data(hmgd.h1.ratio[filter],
                    y=log2(h1.kd.expr[filter]) - log2(s2.expr[filter]),
                    n.bin=n.bin)

hmgd.bins <- bin.data(hmgd.h1.ratio[filter],
                      y=log2(hmgd.kd.expr[filter])  - log2(s2.expr[filter]),
                      n.bin=n.bin)

hmgd.ci.low <- hmgd.bins$y.mean - 2*hmgd.bins$y.sem
hmgd.ci.up  <- hmgd.bins$y.mean + 2*hmgd.bins$y.sem

h1.ci.low <- h1.bins$y.mean - 2*h1.bins$y.sem
h1.ci.up  <- h1.bins$y.mean + 2*h1.bins$y.sem

plot(h1.bins$x.mean, h1.bins$y.mean, ylim=range(hmgd.ci.up, hmgd.ci.low, h1.ci.up, h1.ci.low),
     xlab="HMGD / H1 ratio (log2 RPKM ratio)",
     ylab="expression change in response to knockdown (log2 RPKM ratio)")

points(hmgd.bins$x.mean, hmgd.bins$y.mean, col="red")

segments(x0=-2, x1=3.5, y0=0, y1=0, col="grey")

segments(x0=h1.bins$x.mean, x1=h1.bins$x.mean,
         y0=h1.ci.up, y1=h1.ci.low)

segments(x0=hmgd.bins$x.mean, x1=hmgd.bins$x.mean,
         y0=hmgd.ci.up, y1=hmgd.ci.low, col="red")

legend("topright", legend=c("H1 knockdown", "HMGD knockdown"),
       pch=21, col=c("black", "red"))
       

dev.off()








n.bin <- 5

h1.bins <- bin.data(hmgd.h1.ratio[filter],
                    y=log2(h1.kd.expr[filter]) - log2(s2.expr[filter]),
                    n.bin=n.bin)

hmgd.bins <- bin.data(hmgd.h1.ratio[filter],
                      y=log2(hmgd.kd.expr[filter])  - log2(s2.expr[filter]),
                      n.bin=n.bin)

hmgd.ci.low <- hmgd.bins$y.mean - 2*hmgd.bins$y.sem
hmgd.ci.up  <- hmgd.bins$y.mean + 2*hmgd.bins$y.sem

h1.ci.low <- h1.bins$y.mean - 2*h1.bins$y.sem
h1.ci.up  <- h1.bins$y.mean + 2*h1.bins$y.sem

plot(h1.bins$x.mean, h1.bins$y.mean, ylim=range(h1.ci.up, h1.ci.low),
     xlab="HMGD / H1 ratio (log2 RPKM ratio)",
     ylab="expression change in response to knockdown (log2 RPKM ratio)")

points(hmgd.bins$x.mean, hmgd.bins$y.mean, col="#008F00")

segments(x0=-2, x1=3.5, y0=0, y1=0, col="grey")

segments(x0=h1.bins$x.mean, x1=h1.bins$x.mean,
         y0=h1.ci.up, y1=h1.ci.low)

segments(x0=hmgd.bins$x.mean, x1=hmgd.bins$x.mean,
         y0=hmgd.ci.up, y1=hmgd.ci.low, col="red")

legend("topright", legend=c("H1 knockdown"),
       pch=21, col=c("black"))

dev.off()





pdf("hmgd_h1_kd_expression_change.pdf", width=10, height=4)

par(mfrow=c(1,3))

n.bin <- 50

filter <- hmgd.rpkm < 500 & h1.rpkm < 500 & mnase.rpkm < 500 &
          !is.na(hmgd.rpkm) & !is.na(h1.rpkm) & !is.na(mnase.rpkm) & is.expressed


gc.content <- tab$EXON.GC.COUNT / tab$TOTAL.EXON.LEN

h1.mdl <- lm(log2(h1.kd.expr[filter]) ~ gc.content[filter] + log2(s2.expr[filter]))
h1.resid <- h1.mdl$resid

h1.bins <- bin.data(h1.resid, y=h1.rpkm[filter], n.bin=n.bin)
h1.ci.low <- h1.bins$y.mean - 2*h1.bins$y.sem
h1.ci.up  <- h1.bins$y.mean + 2*h1.bins$y.sem


plot(h1.bins$x.mean, h1.bins$y.mean, ylim=range(h1.ci.up, h1.ci.low),
     xlab="relative expression change (log2 RPKM)",
     ylab="H1 RPKM")

segments(x0=h1.bins$x.mean, x1=h1.bins$x.mean,
         y0=h1.ci.up, y1=h1.ci.low)


hmgd.mdl <- lm(log2(hmgd.kd.expr[filter]) ~ gc.content[filter] + log2(s2.expr[filter]))
hmgd.resid <- hmgd.mdl$resid

hmgd.bins <- bin.data(hmgd.resid, y=hmgd.rpkm[filter], n.bin=n.bin)
hmgd.ci.low <- hmgd.bins$y.mean - 2*hmgd.bins$y.sem
hmgd.ci.up  <- hmgd.bins$y.mean + 2*hmgd.bins$y.sem


plot(hmgd.bins$x.mean, hmgd.bins$y.mean, ylim=range(hmgd.ci.up, hmgd.ci.low),
     xlab="relative expression change (log2 RPKM)",
     ylab="HMGD RPKM")

segments(x0=hmgd.bins$x.mean, x1=hmgd.bins$x.mean,
         y0=hmgd.ci.up, y1=hmgd.ci.low)


ratio.h1.bins <- bin.data(h1.resid, y=hmgd.h1.ratio[filter], n.bin=n.bin)
ratio.h1.ci.low <- ratio.h1.bins$y.mean - 2*ratio.h1.bins$y.sem
ratio.h1.ci.up  <- ratio.h1.bins$y.mean + 2*ratio.h1.bins$y.sem

ratio.hmgd.bins <- bin.data(hmgd.resid, y=hmgd.h1.ratio[filter], n.bin=n.bin)
ratio.hmgd.ci.low <- ratio.hmgd.bins$y.mean - 2*ratio.hmgd.bins$y.sem
ratio.hmgd.ci.up  <- ratio.hmgd.bins$y.mean + 2*ratio.hmgd.bins$y.sem


plot(ratio.h1.bins$x.mean, ratio.h1.bins$y.mean,
     ylim=range(ratio.h1.ci.up, ratio.h1.ci.low, ratio.hmgd.ci.up, ratio.hmgd.ci.low),
     xlab="relative expression change (log2 RPKM)",
     ylab="HMGD1 / H1 ratio (log2 RPKM)", las=1)

segments(x0=ratio.h1.bins$x.mean, x1=ratio.h1.bins$x.mean,
         y0=ratio.h1.ci.up, y1=ratio.h1.ci.low)


points(ratio.hmgd.bins$x.mean, ratio.hmgd.bins$y.mean,
       col="red")

segments(x0=ratio.hmgd.bins$x.mean, x1=ratio.hmgd.bins$x.mean,
         y0=ratio.hmgd.ci.up, y1=ratio.hmgd.ci.low,
         col="red")

legend("topright", legend=c("H1 KD", "HMGD KD"), pch=21, col=c("black", "red"))


dev.off()

