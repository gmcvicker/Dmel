

source("~/proj/script/R/lib/util.R")

source("~/proj/script/R/10_IND/heatmap.R")


summary.tab <- read.table("~/data/Dmel/nuc_repeat_len_by_region/region_summary.txt.gz", header=T)

scores.tab <- read.matrix("~/data/Dmel/nuc_repeat_len_by_region/scores.txt.gz")

mnase.tab <- read.matrix("~/data/Dmel/nuc_repeat_len_by_region/mnase.txt.gz")



hmgd.ttl <- summary.tab$HMGD.TTL
h1.ttl <- summary.tab$H1.TTL

hmgd.h1.ratio <- log2((hmgd.ttl+1) / (h1.ttl+1))



q <- quantile(hmgd.h1.ratio, 0.9)
f <- hmgd.h1.ratio > q

# f <- rep(TRUE, nrow(scores.tab))
mid <- 1001
score.sum <- apply(scores.tab[f,], 2, sum)
mnase.sum <- apply(mnase.tab[f,], 2, sum)

# fold counts about midpoint
score.sum <- score.sum[mid:length(score.sum)] + rev(score.sum[1:mid])
mean.score <- score.sum / sum(f)  / 2

mnase.sum <- mnase.sum[mid:length(mnase.sum)] + rev(mnase.sum[1:mid])
mean.mnase <- mnase.sum / sum(f)

plot(mean.score, type="l")

repeat.len <- 180
offset <- 100
peaks <- seq(from=offset, to=length(scores), by=repeat.len)




f <- rep(TRUE, nrow(mnase.tab))
mnase.sum <- apply(mnase.tab[f,], 2, sum)
mean.mnase <- mnase.sum / sum(f)

plot(mean.mnase)



plot(scores.tab[2,])

scores <- scores.tab[3,]

threshold <- 0.2
peak.regions <- rle(scores > threshold)

n.peaks <- 0
start <- 1

peaks <- rep(NA, sum(peak.regions$values))

for(i in 1:length(peak.regions$lengths)) {
  end <- start + peak.regions$lengths[i] - 1

  is.peak <- peak.regions$values[i]

  if(is.peak) {
    max.score <- max(scores[start:end])
    n.peaks <- n.peaks + 1
    peaks[n.peaks] <- start + which(scores[start:end] == max.score)[1] - 1
  }
  
  start <- end+1
}

plot(scores, type="l")
points(peaks, scores[peaks], col="red")



repeat.len <- 180
offset <- 100
peaks <- seq(from=offset, to=length(scores), by=repeat.len)


plot(scores, type="l")
points(peaks, scores[peaks], col="red")


plot(mnase.tab[3,])
