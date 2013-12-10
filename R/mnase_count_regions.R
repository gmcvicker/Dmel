
chisq.pval <- function(count.row) {  
  count.matrix <- matrix(count.row, nrow=2, ncol=2)
  result <- chisq.test(count.matrix)
  return(result$p.value)
}


p <- function(...) {
  paste(..., sep='')
}


lineage <- "haltere"

tab <- read.table(p("~/data/Dmel/MNase/counts_per_region.", lineage, ".txt"),
                    header=F)

combined.count <- tab$V4 + tab$V6

filter <- combined.count >= 20
tab <- tab[filter,]


lineage.counts <- tab$V4
total.lineage.counts <- tab$V5
other.counts <- tab$V6
other.total.counts <- tab$V7

count.tab <- tab[,4:7]


p.vals <- apply(count.tab, 1, chisq.pval)


ratio <- log2((lineage.counts / total.lineage.counts) / (other.counts / other.total.counts))


hist(ratio[p.vals < 1e-4])

high.mnase <- (p.vals < 1e-4) & (ratio > 0)
low.mnase <- (p.vals < 1e-4) & (ratio < 0)

high.mnase.tab <- data.frame(CHROM=tab$V1[high.mnase],
                        START=tab$V2[high.mnase],
                        END=tab$V3[high.mnase],
                        LIN.MNASE.COUNT=lineage.counts[high.mnase],
                        TOT.LIN.MNASE.COUNT=total.lineage.counts[high.mnase],
                        OTHER.MNASE.COUNT=other.counts[high.mnase],
                        TOT.OTHER.MNASE.COUNT=other.total.counts[high.mnase],
                        LOG2.RATE.RATIO=ratio[high.mnase],
                        P.VAL=p.vals[high.mnase])

write.table(high.mnase.tab,
            file=p("~/data/Dmel/MNase/high_mnase_regions.", lineage, ".txt"),
            quote=FALSE,
            sep=" ",
            col.names=TRUE,
            row.names=FALSE)


low.mnase.tab <- data.frame(CHROM=tab$V1[low.mnase],
                        START=tab$V2[low.mnase],
                        END=tab$V3[low.mnase],
                        LIN.MNASE.COUNT=lineage.counts[low.mnase],
                        TOT.LIN.MNASE.COUNT=total.lineage.counts[low.mnase],
                        OTHER.MNASE.COUNT=other.counts[low.mnase],
                        TOT.OTHER.MNASE.COUNT=other.total.counts[low.mnase],
                        LOG2.RATE.RATIO=ratio[low.mnase],
                        P.VAL=p.vals[low.mnase])

write.table(low.mnase.tab,
            file=p("~/data/Dmel/MNase/low_mnase_regions.",
              lineage, ".txt"),
            quote=FALSE,
            sep=" ",
            col.names=TRUE,
            row.names=FALSE)



############################




kmer.tab <- read.table(p("~/data/Dmel/MNase/low_mnase_region_kmer_counts.", lineage, ".txt"), header=T)

tot.kmers <- sum(kmer.tab$ALL.COUNT)
tot.region.kmers <- sum(kmer.tab$REGION.COUNT)


pseudo.count <- 1

kmer.count.tab <- cbind(kmer.tab$REGION.COUNT + pseudo.count,
                        tot.region.kmers + pseudo.count,
                        kmer.tab$ALL.COUNT + pseudo.count,
                        tot.kmers + pseudo.count)

p.vals <- apply(kmer.count.tab, 1, chisq.pval)
log.ratio <- log2(((kmer.tab$REGION.COUNT + pseudo.count) /
                    (tot.region.kmers + pseudo.count)) /
                  ((kmer.tab$ALL.COUNT + pseudo.count) /
                   (tot.kmers + pseudo.count)))


lineage.rate <- (kmer.tab$REGION.COUNT + pseudo.count) /
  (tot.region.kmers + pseudo.count)

genome.wide.rate <- (kmer.tab$ALL.COUNT + pseudo.count) /
  (tot.kmers + pseudo.count)

plot(lineage.rate, genome.wide.rate)

ord <- order(p.vals)
kmer.tab[['P.VAL']] <- p.vals
kmer.tab[['LOG.OBS.EXP.RATIO']] <- log.ratio

enriched.tab <- kmer.tab[log.ratio > 0.0,]
depleted.tab <- kmer.tab[log.ratio < 0.0,]

head(enriched.tab[order(enriched.tab$P.VAL),], 50)
head(depleted.tab[order(depleted.tab$P.VAL),], 50)

write.table(head(enriched.tab[order(enriched.tab$P.VAL),], 50),
            file=p("kmers_enriched_in_low_mnase.",lineage,".txt"),
            sep="\t", col.names=T, row.names=F, quote=F)

write.table(head(depleted.tab[order(depleted.tab$P.VAL),], 50),
            file=p("kmers_depleted_in_low_mnase.",lineage,".txt"),
            sep="\t", col.names=T, row.names=F, quote=F)


f <- p.vals < 1e-100
text(lineage.rate[f], genome.wide.rate[f],
     labels=kmer.tab$KMER[f], pos=4, col="red")
