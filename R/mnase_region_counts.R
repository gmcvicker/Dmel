
run.tests <- function(tab1, tab2) { 
  f <- (tab1$COUNT + tab2$COUNT) >= 20

  new.tab <- cbind(tab1$COUNT, tab1$TOTAL, tab2$COUNT, tab2$TOTAL)
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

# s2.tab <- read.table(paste(d, "/S2.region_counts.txt.gz", sep=""), header=T)

in.vitro.tab <- read.table(paste(d, "/S2_in_vitro_combined.region_counts.txt.gz", sep=""), header=T)


tab1 <- in.vitro.tab

in.vitro.fpkm <- as.numeric(in.vitro.tab$COUNT+1)*1e9 /
  ((in.vitro.tab$END - in.vitro.tab$START)*as.numeric(in.vitro.tab$TOTAL.COUNTS))


for(lineage in c("S2", "eye", "leg", "haltere", "antenna")) {
  cat(lineage, "\n")
  tab2 <- read.table(paste(d, "/", lineage, ".region_counts.txt.gz", sep=""),
                     header=T)

  fpkm <- as.numeric(tab2$COUNT+1)*1e9 /
    ((tab2$END - tab2$START) * as.numeric(tab2$TOTAL.COUNTS))

  p.values <- run.tests(tab1, tab2)

  png(paste(out.d, "/", "in_vitro_vs_", lineage, "_200bp_fpkm.png", sep=""),
      width=500, height=500)

  f <- !is.na(p.values)
  fdr <- rep(NA, length(p.values))
  fdr[f] <- p.adjust(p.values[f], method="BH")

  f <- TRUE
  plot(log10(in.vitro.fpkm[f]), log10(fpkm[f]), pch=21,
       bg="grey", col="grey", xlab="in vitro log10(FPKM)",
       ylab=paste(lineage, "log10(FPKM))"))

  log.ratio <- log10(fpkm) - log10(in.vitro.fpkm)
  
  f <- !is.na(fdr) & (fdr <= 0.01) & (log.ratio > 0.0)
  
  new.tab <- cbind(tab2, IN.VITRO.COUNT=in.vitro.tab$COUNT,
                   IN.VITRO.TOTAl=in.vitro.tab$TOTAL.COUNTS,
                   P.VALUE=p.values, FDR=fdr)
  
  output.name <- paste(out.d, "/", lineage, "_fdr0.01_higher.txt", sep="")
  write.table(new.tab[f,], file=output.name, quote=F, row.names=F,
              col.names=T, sep=" ")
  
  points(log10(in.vitro.fpkm[f]), log10(fpkm[f]), pch=21,
         bg="red", col="red")
  
  f <- !is.na(fdr) & (fdr <= 0.01) & (log.ratio < 0.0)

  output.name <- paste(out.d, "/", lineage, "_fdr0.01_lower.txt", sep="")
  write.table(new.tab[f,], file=output.name, quote=F,
              row.names=F, col.names=T, sep=" ")

  points(log10(in.vitro.fpkm[f]), log10(fpkm[f]), pch=21,
         bg="blue", col="blue")

  abline(a=0, b=1)
  
  dev.off()
}





