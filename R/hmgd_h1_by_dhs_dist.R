
library(ggplot2)

tab <- read.table("/home/gmcvicker/data/Dmel/DNase/hmgd_h1_by_dhs_dist.txt",
                  header=T)

read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")
h1.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]
hmgd.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "hmgd"]
mnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]
dnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "dnase"]


has.data <- tab$HMGD.COUNT > 0 &
             tab$DNASE.COUNT > 0 &
             tab$H1.COUNT > 0 &
             tab$MNASE.COUNT > 0

is.outlier <- tab$HMGD.COUNT > quantile(tab$HMGD.COUNT, 0.9999) |
              tab$DNASE.COUNT > quantile(tab$HMGD.COUNT, 0.9999) |
              tab$H1.COUNT > quantile(tab$H1.COUNT, 0.9999) |
              tab$MNASE.COUNT > quantile(tab$MNASE.COUNT, 0.9999)

het.chrom.names <- grep("Het", unique(tab$CHROM), value=T)

is.het.chrom <- tab$CHROM %in% het.chrom.names

is.sex.chrom <- tab$CHROM %in% c("chrX", "chrY")

win.sz <- tab$END[1]-tab$START[1]+1

dat <- data.frame(MEAN.DHS.DIST=tab$MEAN.DHS.DIST,
                  MEAN.TSS.DIST=tab$MEAN.TSS.DIST,
                  MNASE.RPKM=tab$MNASE.COUNT*1e9 / (win.sz*mnase.ttl),
                  DNASE.RPKM=tab$DNASE.COUNT*1e9 / (win.sz*dnase.ttl),
                  HMGD.RPKM=tab$HMGD.COUNT*1e9 / (win.sz*hmgd.ttl),
                  H1.RPKM=tab$H1.COUNT*1e9 / (win.sz*hmgd.ttl))
                  
f <- has.data & !is.outlier & !is.sex.chrom & !is.het.chrom
dat <- dat[f,]

dat["HMGD.H1.RATIO"] <- log2(dat$HMGD.RPKM / dat$H1.RPKM)

dat["DHS.DIST.CLASS"] <- cut(dat$MEAN.DHS.DIST,
                             breaks=c(0,2 * 10**seq(2, 8, by=1)),
                             include.lowest=T)

dat["TSS.DIST.CLASS"] <- cut(dat$MEAN.TSS.DIST,
                             breaks=c(0,2000, 2e6),
                             include.lowest=T)

dat["DNASE.CLASS"] <- cut(dat$DNASE.RATE, 5)


region.sz <- 1000

dat["ANNOTATION"] <- factor(levels=c("intergenic", "intron", "tss.upstream",
                              "tss.downstream", "exon"))

dat$ANNOTATION[tab$INTRON.COUNT[f] == region.sz] <- "intron"
dat$ANNOTATION[tab$INTERGENIC.COUNT[f] == region.sz] <- "intergenic"
dat$ANNOTATION[tab$EXON.COUNT[f] > 0] <- "exon"
dat$ANNOTATION[tab[f,"TSS_UPSTREAM.COUNT"] > 0] <- "tss.upstream"
dat$ANNOTATION[tab[f,"TSS_DOWNSTREAM.COUNT"] > 0] <- "tss.downstream"


dat["OVERLAP.INTERGENIC"] <- tab$INTERGENIC.COUNT[f] > 0
dat["OVERLAP.EXON"] <- tab$EXON.COUNT[f] > 0
dat["OVERLAP.TSS.UP"] <- tab[,"TSS_UPSTREAM.COUNT"][f] > 0
dat["OVERLAP.TSS.DOWN"] <- tab[,"TSS_DOWNSTREAM.COUNT"][f] > 0
dat["OVERLAP.INTRON"] <- tab[,"INTRON.COUNT"][f] > 0
dat["ALL.INTRON"] <- tab[,"INTRON.COUNT"][f] == 1000
dat["ALL.INTERGENIC"] <- tab[,"INTERGENIC.COUNT"][f] == 1000


f1 <- as.numeric(dat$DHS.DIST.CLASS) == 1
f2 <- as.numeric(dat$DHS.DIST.CLASS) == 2
f3 <- as.numeric(dat$DHS.DIST.CLASS) == 3
f4 <- as.numeric(dat$DHS.DIST.CLASS) == 4
f5 <- as.numeric(dat$DHS.DIST.CLASS) == 5

types <- c("HMGD.RATE", "H1.RATE", "HMGD.H1.RATIO")

pdf("DHS_dist_boxplots.pdf", width=8, height=4)

par(mfrow=c(1,3))
for(i in 1:length(types)) {  
  boxplot(dat[f1, types[i]],
          dat[f2, types[i]],
          dat[f3, types[i]],
          dat[f4, types[i]],
          dat[f5, types[i]],
          main=types[i],
          outline=F, las=1, col="grey")
}

dev.off()



pdf("DHS_dist_and_TSS_dist_boxplots.pdf")

par(mfrow=c(1,2))

tss.near <- as.numeric(dat$TSS.DIST.CLASS) == 1
tss.far <- as.numeric(dat$TSS.DIST.CLASS) == 2

boxplot(dat[f1 & tss.near, "HMGD.H1.RATIO"],
        dat[f2 & tss.near, types[i]],
        dat[f3 & tss.near, types[i]],
        dat[f4 & tss.near, types[i]],
        dat[f5 & tss.near, types[i]],
        main="near TSS",
        outline=F, las=1, col="grey")

boxplot(dat[f1 & tss.far, "HMGD.H1.RATIO"],
        dat[f2 & tss.far, types[i]],
        dat[f3 & tss.far, types[i]],
        dat[f4 & tss.far, types[i]],
        dat[f5 & tss.far, types[i]],
        main="far TSS",
        outline=F, las=1, col="grey")

dev.off()









pdf("DHS_dist_HMGD_boxplot.pdf", width=6, height=6)
ggplot(dat, aes(DHS.DIST.CLASS, HMGD.RATE)) + geom_boxplot() +
  ylab("HMGD rate (midpoints per million mapped reads)") +
  xlab("DNase hypersensitive site distance (bp)")
dev.off()

pdf("DHS_dist_H1_boxplot.pdf", width=6, height=6)
ggplot(dat, aes(DHS.DIST.CLASS, H1.RATE)) + geom_boxplot() +
  ylab("H1 rate (midpoints per million mapped reads)") +
  xlab("DNase hypersensitive site distance (bp)")
dev.off()

pdf("DHS_dist_HMGD_H1_ratio_boxplot.pdf", width=6, height=6)
ggplot(dat, aes(DHS.DIST.CLASS, HMGD.H1.RATIO)) +
  geom_boxplot() +
  ylab("HMGD/H1 rate ratio (log2)") +
  xlab("DNase hypersensitive site distance (bp)")
dev.off()

pdf("DHS_dist_HMGD_H1_ratio_boxplot_by_TSS_dist.pdf", width=8, height=6)
ggplot(dat, aes(DHS.DIST.CLASS, HMGD.H1.RATIO)) +
  geom_boxplot() +
  facet_grid(~ TSS.DIST.CLASS) +
  ylab("HMGD/H1 rate ratio (log2)") +
  xlab("DNase hypersensitive site distance (bp)")
dev.off()


pdf("DHS_dist_HMGD_H1_ratio_boxplot_by_annotation.pdf", width=8, height=6)
ggplot(dat, aes(ANNOTATION, HMGD.H1.RATIO)) +
  geom_boxplot() +
  ylab("HMGD/H1 rate ratio (log2)") +
  xlab("Annotation type")
#dev.off()


# pdf("DHS_dist_HMGD_H1_ratio_boxplot_by_TSS_dist.pdf", width=8, height=6)
ggplot(dat, aes(ANNOTATION, HMGD.RATE)) +
  geom_boxplot() +
  ylab("HMGD/H1 rate ratio (log2)") +
  xlab("Annotation type")
#dev.off()

# pdf("DHS_dist_HMGD_H1_ratio_boxplot_by_TSS_dist.pdf", width=8, height=6)
ggplot(dat, aes(ANNOTATION, H1.RATE)) +
  geom_boxplot() +
  ylab("HMGD/H1 rate ratio (log2)") +
  xlab("Annotation type")
#dev.off()





cor.test(dat$MEAN.DHS.DIST, dat$HMGD.H1.RATIO)

cor.test(dat$MEAN.DHS.DIST, log2(dat$HMGD.RATE))
cor.test(dat$MEAN.DHS.DIST, log2(dat$H1.RATE))
cor.test(dat$MEAN.DHS.DIST, dat$HMGD.H1.RATIO)



cor.test(dat$DNASE.RATE, dat$HMGD.H1.RATIO)


