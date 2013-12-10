
library(RColorBrewer)
color.pal <- c("grey50", brewer.pal(9, "Set1"))

h1.tab <- read.table("~/data/Dmel/reads_per_chrom/H1_combined.txt",
                     header=T)

hmgd.tab <- read.table("~/data/Dmel/reads_per_chrom/HMGD_combined.txt",
                     header=T)

mnase.tab <- read.table("~/data/Dmel/reads_per_chrom/MNase_S2_combined.txt",header=T)

get.rate <- function(tab) {
  (tab$MIDPOINTS.PER.SITE * 1e6) / sum(tab$N.MIDPOINTS) 
}

mnase.tab['RATE'] <- get.rate(mnase.tab)
h1.tab['RATE'] <- get.rate(h1.tab)
hmgd.tab['RATE'] <- get.rate(hmgd.tab)


h1.tab["H1.REL.RATE"] <- h1.tab$RATE / mnase.tab$RATE
hmgd.tab["HMGD.REL.RATE"] <- hmgd.tab$RATE / mnase.tab$RATE


pdf("hmgd_h1_mids_by_chrom.pdf", width=10, height=4)

barplot(rbind(mnase.tab$RATE, h1.tab$RATE, hmgd.tab$RATE),
        names=h1.tab$CHROM, beside=T, las=1,
        legend=c("nucleosomes", "H1", "HMGD"),
        col=color.pal[1:3], ylab="mids per base per 10^6 mapped reads")

dev.off()


is.het <- rep(FALSE, nrow(mnase.tab))
is.het[grep("Het", mnase.tab$CHROM)] <- TRUE

is.euc <- grep("Het", mnase.tab$CHROM, invert=T)


ttl.het.reads <- sum(mnase.tab[is.het,]$N.MIDPOINTS)
ttl.euc.reads <- sum(mnase.tab[is.euc,]$N.MIDPOINTS)

het.h1.rate <- h1.tab[is.het,]$N.MIDPOINTS / sum(h1.tab[is.het,]

                                                 
