


read.kmer.counts <- function(lineages, kmer.size=7) {
  n.kmer <- 4**kmer.size / 2

  kmer.tab <- data.frame(KMER.NUM=1:n.kmer)

  for(lin in lineages) {
    filename <- paste("~/data/Dmel/MNase/kmer_counts/", lin, ".",
                      kmer.size, "mer.txt.gz",
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


experiment.tab <- read.table("~/data/Dmel/mnase_midpoint_tracks.txt",
                             col.names=c("NAME","TRACK"))

experiments <- as.character(experiment.tab$NAME)
kmer.tab <- read.kmer.counts(experiment)

n.exper <- length(experiments)
exper.totals <- rep(NA, n.exper)

for(i in 1:n.exper) {  
  exper <- experiments[i]
  exper.totals[i] <- sum(as.numeric(kmer.tab[,exper]))

  rate.label <- paste(exper, ".RATE", sep="")
  
  kmer.tab[rate.label] <- kmer.tab[,exper] /
    (exper.totals[i] * kmer.tab$KMER.N.OCCURANCES)
}


combined.count <- apply(kmer.tab[,experiments], 1, sum)

in.vitro.rate <- kmer.tab[,"S2_in_vitro_combined.RATE"]

in.vivo.exper <- c("S2", "antenna", "eye", "haltere", "leg")

# combined.rate <- combined.count / (kmer.tab$KMER.N.OCCURANCES * sum(lineage.totals))

pdf("mnase_kmer_rate_scatter.pdf", width=5, height=5)

for(exper in in.vivo.exper) {  
  label <- paste(exper, ".RATE", sep="")

  in.vivo.rate <- kmer.tab[[label]]

  par(xpd=FALSE)
  
  plot(log10(in.vitro.rate), log10(in.vivo.rate), col="grey50",
       ylab="log10(in vivo rate)", xlab="log10(in vitro rate)",
       main=exper)

  mdl <- lm(log10(in.vivo.rate) ~ log10(in.vitro.rate))

  abline(a=mdl$coef[1], b=mdl$coef[2])
  # abline(a=0, b=1)

  rank.stat <- -abs(mdl$resid)
  f <- rank(rank.stat) <= 20


  text(log10(in.vitro.rate[f]),
       log10(in.vivo.rate[f]),
       labels=kmer.tab$KMER[f], pos=4, col="red",
       cex=0.5)
  
}

dev.off()












pdf("mnase_kmer_rate_scatter.pdf", width=5, height=5)

cor.matrix <- matrix(nrow=length(experiments), ncol=length(experiments))
diag(cor.matrix) <- 1.0

for(i in 1:length(experiments)) {
  for(j in 1:i) {
    if(i == j) {
      next
    }
    exper1 <- experiments[i]
    exper2 <- experiments[j]
    label1 <- paste(exper1, ".RATE", sep="")
    label2 <- paste(exper2, ".RATE", sep="")

    cat(exper1, "|", exper2, "\n")
      
    rate1 <- kmer.tab[[label1]]
    rate2 <- kmer.tab[[label2]]

    par(xpd=FALSE)

    res <- cor.test(log10(rate1), log10(rate2))
    main <- paste(exper1, " vs. ", exper2, " (R=",
                  signif(res$estimate,3), ")", sep="")

    cor.matrix[i,j] <- res$est
    cor.matrix[j,i] <- res$est
    
    plot(log10(rate1), log10(rate2), col="grey50",
         ylab=paste(exper1, "log10(rate)"),
         xlab=paste(exper2, "log10(in vitro rate)"),
         main=main)

    mdl <- lm(log10(rate2) ~ log10(rate1))
    abline(a=mdl$coef[1], b=mdl$coef[2])    
    rank.stat <- -abs(mdl$resid)
    f <- rank(rank.stat) <= 20
    
    text(log10(rate1[f]),
         log10(rate2[f]),
         labels=kmer.tab$KMER[f], pos=4, col="red",
         cex=0.5)

    
  }
}

dev.off()

dist.matrix <- 1-cor.matrix^2
rownames(dist.matrix) <- experiments
colnames(dist.matrix) <- experiments
clust <- hclust(dist(dist.matrix))

pdf("mnase_kmer_corr_clustering.pdf", width=5, height=5)
plot(clust)
dev.off()
