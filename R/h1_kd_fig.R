
#
# Read RNA-seq expression
#
tab <- read.table("~/data/Dmel/KD_RNA_seq/tr_rna_seq_counts.txt", header=T)
expr <- tab[,c("RNASEQ.RPKM.S2_KD.Wt", "RNASEQ.RPKM.S2_KD.H1KD")]

#
# read total number of mapped H1, and total nucleosome reads
#
read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")
h1.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]
mnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]

#
# Read number of mapped H1, and total nucleosome reads at each promoter
#
hmgd.h1.tab <- read.table("~/data/Dmel/KD_RNA_seq/hmgd_h1_at_rnaseq_promoters.txt", header=T)
tss.region.len <- hmgd.h1.tab$END - hmgd.h1.tab$START + 1
hmgd.rpkm <- ((hmgd.h1.tab$HMGD.COUNT+1) * 1e9) / (tss.region.len * hmgd.ttl)
h1.rpkm <- ((hmgd.h1.tab$H1.COUNT+1) * 1e9) / (tss.region.len * h1.ttl)
mnase.rpkm <- ((hmgd.h1.tab$MNASE.COUNT+1) * 1e9) / (tss.region.len * mnase.ttl)

hmgd.h1.ratio <- log2(hmgd.rpkm / h1.rpkm)

s2.expr <- expr[["RNASEQ.RPKM.S2_KD.Wt"]]
h1.kd.expr <- expr[["RNASEQ.RPKM.S2_KD.H1KD"]]

# require above a minimum expression threshold in at least one experiment
min.expr <- 5

pdf("h1_kd_expr_histograms.pdf", width=15, height=5)

par(mfrow=c(1,2))

hist(log2(s2.expr), breaks=50, main="S2 expr")
segments(x0=log2(min.expr), x1=log2(min.expr), y0=0, y1=2500, col="red")

hist(log2(h1.kd.expr), breaks=50, main="H1 KD expr",
     xlab="expression (log2 RPKM)")
segments(x0=log2(min.expr), x1=log2(min.expr), y0=0, y1=2500, col="red")

dev.off()

is.expressed <- !is.na(s2.expr) & !is.na(h1.kd.expr) & 
                ((s2.expr > min.expr) | (h1.kd.expr > min.expr))

                

# there are a handful of outlier regions (probably collapsed repeats) with far more read counts
# than all other regions
filter <- h1.rpkm < 500 & mnase.rpkm < 500 &
          !is.na(h1.rpkm) & !is.na(mnase.rpkm) & is.expressed


################################

h1.kd.expr.diff <- log2(h1.kd.expr) - log2(s2.expr)

source("~/proj/script/R/lib/util.R")

n.bin <- 5

h1.bins <- bin.data(h1.rpkm[filter],
                    y=h1.kd.expr.diff[filter],
                    n.bin=n.bin)

h1.ci.low <- h1.bins$y.mean - 2*h1.bins$y.sem
h1.ci.up  <- h1.bins$y.mean + 2*h1.bins$y.sem

plot(h1.bins$x.mean, h1.bins$y.mean, ylim=range(h1.ci.up, h1.ci.low),
     xlab="H1 promoter occupancy (RPKM)", col="#008F00",
     bg="#008F00", xlim=c(0, 12),
     ylab="mean expression change after H1 knockdown (log2 RPKM ratio)",
     pch=21)

segments(x0=-2, x1=3.5, y0=0, y1=0, col="grey")

segments(x0=h1.bins$x.mean, x1=h1.bins$x.mean,
         y0=h1.ci.up, y1=h1.ci.low, col="#008F00",)

dev.off()

cor.test(h1.rpkm[filter], h1.kd.expr.diff[filter])
