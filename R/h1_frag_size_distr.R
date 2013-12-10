
get.combined.hist <- function(tab) {
  hist.tab <- NULL
  
  for(i in 1:nrow(tab)) {
    sample <- tab$V3[i]
    path <- paste(dir, "/frag_size_distr/", sample, ".txt", sep="")

    sample.hist.tab <- read.table(path, header=F, #skip=2,
                                  col.names=c("FRAG.LEN", "COUNT", "PRP"))

    if(is.null(hist.tab)) {
      hist.tab <- sample.hist.tab[,1:2]
    } else {
      hist.tab$COUNT <- as.numeric(sample.hist.tab$COUNT) +
        as.numeric(hist.tab$COUNT)
    }
  }

  return(hist.tab)
}


library(RColorBrewer)

dir <- "~/data/Dmel/H1"
summary.tab <- read.table(paste(dir, "/tasklist.txt", sep=""), header=F)

sample.names <- summary.tab$V3

pdf("frag_size_hist.h1.pdf", width=10, height=8)

colors <- brewer.pal(length(sample.names), "Set1")

plot(x=c(0), y=c(0), type="n", ylim=c(0, 0.03), xlim=c(90, 220),
     ylab="fraction of fragments", xlab="fragment size", las=1)

sample.medians <- rep(NA, length(sample.names))
sample.means <- rep(NA, length(sample.names))
sample.sds <- rep(NA, length(sample.names))

for(i in 1:length(sample.names)) {
  
  filename <- paste(dir, "/frag_size_distr/", sample.names[i], ".txt", sep="")
  hist.tab <- read.table(filename, header=F)
  frag.len <- hist.tab$V1
  count <- hist.tab$V2
  
  ttl.reads <- sum(count)
  frac <- count / ttl.reads
  cm.prp <- cumsum(frac)
  sample.medians[i] <- frag.len[cm.prp >= 0.5][1]
  sample.means[i] <- sum(frac * frag.len)
  sample.sds[i] <- sqrt(sum(count * (frag.len - sample.means[i])^2) / sum(count))

  lines(frag.len, frac , type="l", col=colors[i], lwd=2)
}

hist.tab <- get.combined.hist(summary.tab)
count <- hist.tab$COUNT
frag.len <- hist.tab$FRAG.LEN
ttl.reads <- sum(count)
frac <- count / ttl.reads

cm.prp <- cumsum(frac)
# central 95%
x1 <- frag.len[cm.prp >= 0.025][1]
x2 <- rev(frag.len[cm.prp <= 0.975])[1]
med.frag.len <- frag.len[cm.prp >= 0.5][1]
mean.frag.len <- sum(frac * frag.len)

sd.frag.len <- sqrt(sum(count*(frag.len - mean.frag.len)^2) / sum(count))
                   

lines(frag.len, frac , type="l", col="black", lwd=2)

labels <- c(paste(sample.names,
                  " (mean: ", round(sample.means, 1),
                  ", med: ", sample.medians, ")",
                  sep=""),
            paste("combined (mean: ", round(mean.frag.len, 1),
                  ", med: ", med.frag.len, ")", sep=""))

legend("topleft",
       legend=labels,
       col=c(colors, "black"), lty=1, lwd=2,
       title="sample")

dev.off()


