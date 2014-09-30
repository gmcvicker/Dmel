

run.tests <- function(tab1, tab2) { 
  new.tab <- cbind(tab1$COUNT, tab1$TOTAL-tab1$COUNT,
                   tab2$COUNT, tab2$TOTAL-tab2$COUNT)
  p.values <- apply(new.tab, 1, test.chisq)

  return(p.values)
}


test.chisq <- function(x) {
  m <- matrix(x, nrow=2, ncol=2)

  res <- chisq.test(m)
  
  return(res$p.value)
}



cell.lines <- c("eye", "haltere", "antenna", "leg")

d <- "~/data/Dmel/MNase/mnase_region_counts/signif_diff/"

for(cell.line in cell.lines) {

    cat(cell.line, "\n")

    high.kmer.tab <- read.table(paste(d, cell.line, "_fdr0.01_higher_kmers.txt", sep=""), header=T)
    low.kmer.tab <- read.table(paste(d, cell.line, "_fdr0.01_lower_kmers.txt", sep=""), header=T)

    high.ttl.region.kmer <- sum(as.numeric(high.kmer.tab$REGION.COUNT))
    low.ttl.region.kmer <- sum(as.numeric(low.kmer.tab$REGION.COUNT))

    tab1 <- data.frame(COUNT=high.kmer.tab$REGION.COUNT,
                       TOTAL=rep(high.ttl.region.kmer, nrow(high.kmer.tab)))

    tab2 <- data.frame(COUNT=low.kmer.tab$REGION.COUNT,
                       TOTAL=rep(low.ttl.region.kmer, nrow(low.kmer.tab)))

    p.val <- run.tests(tab1, tab2)

    high.rate <- high.kmer.tab$REGION.COUNT / high.ttl.region.kmer
    low.rate <- low.kmer.tab$REGION.COUNT / low.ttl.region.kmer

    log2.ratio <- log2(high.rate / low.rate)


    plot(high.rate, low.rate)

# make table of kmers that are enriched in regions with
# low occupancy (compared to S2) in this cell line
    f <- log2.ratio < 0
    decrease.tab <- data.frame(kmer=high.kmer.tab$KMER[f], high.rate=signif(high.rate[f], 3),
                               low.rate=signif(low.rate[f], 3),
                               log2.ratio=signif(log2.ratio[f], 3),
                               chisq.p.val=signif(p.val[f], 3))
    ord <- order(decrease.tab$log2.ratio)
    decrease.tab <- decrease.tab[ord,]

    out.file <- paste(d, cell.line, "_lower_occ_enriched_kmers.txt", sep="")
    write.table(decrease.tab, file=out.file, row.names=FALSE, quote=FALSE, col.names=TRUE)

# make table of kmers that are enriched in regions with
# high occupancy (compared to S2) in this cell line
    f <- log2.ratio > 0
    increase.tab <- data.frame(kmer=high.kmer.tab$KMER[f], high.rate=signif(high.rate[f], 3),
                               low.rate=signif(low.rate[f], 3), log2.ratio=signif(log2.ratio[f], 3),
                               chisq.p.val=signif(p.val[f], 3))
    ord <- order(increase.tab$log2.ratio, decreasing=T)
    increase.tab <- increase.tab[ord,]

    out.file <- paste(d, cell.line, "_higher_occ_enriched_kmers.txt", sep="")
    write.table(increase.tab, file=out.file, row.names=FALSE, quote=FALSE, col.names=TRUE)
    
    ord <- order(increase.tab$log2.ratio, decreasing=T)
    increase.tab <- increase.tab[ord,]
}


