
run.tests <- function(tab1, tab2) { 
  f <- (tab1$COUNT + tab2$COUNT) >= 20

  new.tab <- cbind(tab1$COUNT, tab1$TOTAL-tab1$COUNT,
                   tab2$COUNT, tab2$TOTAL-tab2$COUNT)
  p.values <- rep(NA, nrow(new.tab))
  p.values[f] <- apply(new.tab[f,], 1, test.chisq)

  return(p.values)
}


test.chisq <- function(x) {
  m <- matrix(x, nrow=2, ncol=2)

  res <- chisq.test(m)
  
  return(res$p.value)
}
           


d <- "~/data/Dmel/MNase/mnase_region_counts"
out.d <- "~/data/Dmel/MNase/mnase_region_counts/signif_diff"

s2.tab <- read.table(paste(d, "/S2.region_counts.txt.gz", sep=""), header=T)

# in.vitro.tab <- read.table(paste(d, "/S2_in_vitro_combined.region_counts.txt.gz", sep=""), header=T)

tab1 <- s2.tab

s2.fpkm <- as.numeric(s2.tab$COUNT+1)*1e9 /
  ((s2.tab$END - s2.tab$START)*as.numeric(s2.tab$TOTAL.COUNTS))


# there are a few regions in S2 and in S2 in vitro that have really extreme counts
# compared to other cell lines. Likely some sort of mapping issue
# due to collapsed repeat or massive amplification
# filter out regions with RPKM > 150
max.fpkm <- 150

# pull out regions with doubling or halving of MNase density (log10(2) = 0.3)
log10.change.threshold <- 0.3


for(lineage in c("S2_in_vitro_combined", "eye", "leg", "haltere", "antenna")) {
  cat(lineage, "\n")
  tab2 <- read.table(paste(d, "/", lineage, ".region_counts.txt.gz", sep=""),
                     header=T)

  fpkm <- as.numeric(tab2$COUNT+1)*1e9 /
    ((tab2$END - tab2$START) * as.numeric(tab2$TOTAL.COUNTS))

  exclude <- s2.fpkm > max.fpkm
  
  p.values <- run.tests(tab1, tab2)

  png(paste(out.d, "/", "S2_vs_", lineage, "_200bp_fpkm.png", sep=""),
      width=500, height=500)

  f <- !is.na(p.values) & !exclude
  fdr <- rep(NA, length(p.values))
  fdr[f] <- p.adjust(p.values[f], method="BH")

  f <- TRUE
  plot(log10(s2.fpkm[f]), log10(fpkm[f]), pch=21,
       bg="grey", col="grey", xlab="S2 log10(FPKM)",
       ylab=paste(lineage, "log10(FPKM))"))

  log.ratio <- log10(fpkm) - log10(s2.fpkm)
  
  f <- !is.na(fdr) & (fdr <= 0.01) & (log.ratio > log10.change.threshold) & !exclude
  
  new.tab <- cbind(tab2, S2.COUNT=s2.tab$COUNT,
                   S2.TOTAl=s2.tab$TOTAL.COUNTS,
                   P.VALUE=signif(p.values,4), FDR=signif(fdr, 4))
  
  output.name <- paste(out.d, "/", lineage, "_fdr0.01_higher.txt", sep="")
  write.table(new.tab[f,], file=output.name, quote=F, row.names=F,
              col.names=T, sep=" ")
  
  points(log10(s2.fpkm[f]), log10(fpkm[f]), pch=21,
         bg="red", col="red")
  
  f <- !is.na(fdr) & (fdr <= 0.01) & (log.ratio < -log10.change.threshold) & !exclude

  output.name <- paste(out.d, "/", lineage, "_fdr0.01_lower.txt", sep="")
  write.table(new.tab[f,], file=output.name, quote=F,
              row.names=F, col.names=T, sep=" ")

  points(log10(s2.fpkm[f]), log10(fpkm[f]), pch=21,
         bg="blue", col="blue")

  abline(a=0, b=1)
  
  dev.off()
}





