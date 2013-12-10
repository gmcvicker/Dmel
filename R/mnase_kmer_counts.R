


read.kmer.counts <- function(lineages, kmer.size=7) {
  n.kmer <- 4**kmer.size / 2

  kmer.tab <- data.frame(KMER.NUM=1:n.kmer)

  for(lin in lineages) {
    filename <- paste("~/data/Dmel/MNase/kmer_counts/", lin, ".",
                      kmer.size, "mer.txt",
                      sep="")

    cat(lin, "\n")
    tab <- read.table(filename, header=T)
    kmer.tab["KMER"] <- tab$KMER
    kmer.tab["REVCOMP.KMER"] <- tab$REVCOMP.KMER
    kmer.tab["KMER.N.OCCURANCES"] <- tab$ALL.COUNT
    kmer.tab[lin] = tab$MNASE.COUNT
  }
  return(kmer.tab)
}

lineage.tab <- read.table("~/data/Dmel/lineage.txt",
                          col.names=c("ID","CELL.LINE","LINEAGE","BARCODE"))

lineages <- as.character(lineage.tab$LINEAGE)
kmer.tab <- read.kmer.counts(lineages)

n.lineage <- length(lineages)
lineage.totals <- rep(NA, n.lineage)

for(i in 1:n.lineage) {  
  lin <- lineages[i]
  lineage.totals[i] <- sum(as.numeric(kmer.tab[,lin]))

  rate.label <- paste(lin, ".RATE", sep="")
  
  kmer.tab[rate.label] <- kmer.tab[,lin] /
    (lineage.totals[i] * kmer.tab$KMER.N.OCCURANCES)
}


combined.count <- apply(kmer.tab[,lineages], 1, sum)


mean.rate <- apply(kmer.tab[,9:12], 1, mean)

# combined.rate <- combined.count / (kmer.tab$KMER.N.OCCURANCES * sum(lineage.totals))

pdf("mnase_kmer_rate_scatter.pdf", width=5, height=5)

for(lin in lineages) {  
  label <- paste(lin, ".RATE", sep="")

  lineage.rate <- kmer.tab[[label]]
  
  plot(log10(mean.rate), log10(lineage.rate), col="grey50",
       ylab="log10(lineage rate)", xlab="log10(combined rate)",
       main=lin)

  mdl <- lm(log10(lineage.rate) ~ log10(mean.rate))

  abline(a=mdl$coef[1], b=mdl$coef[2])
  # abline(a=0, b=1)

  rank.stat <- -abs(mdl$resid)
  f <- rank(rank.stat) <= 20
  
  text(log10(mean.rate[f]),
       log10(lineage.rate[f]),
       labels=kmer.tab$KMER[f], pos=4, col="red",
       cex=0.5)
  
}

dev.off()

