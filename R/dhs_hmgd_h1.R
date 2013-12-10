
library(RColorBrewer)
source("~/proj/script/R/lib/util.R")



color.pal <-  brewer.pal(9, "Set1")


data.dir <- "~/data/Dmel/DNase/dhs_counts/"


read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")

dhs.tab <- read.table(p(data.dir, "/dhs_info.txt.gz"), header=T)

h1.matrix <- read.matrix(p(data.dir, "/h1.txt.gz"))
hmgd.matrix <- read.matrix(p(data.dir, "/hmgd.txt.gz"))
mnase.matrix <- read.matrix(p(data.dir, "/mnase.txt.gz"))


dnase.matrix <- read.matrix(p(data.dir, "/dnase.txt.gz"))

h1.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]
hmgd.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "hmgd"]
mnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]
dnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "dnase"]

read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")



max.count <- 20
mnase.matrix[mnase.matrix > max.count] <- max.count
hmgd.matrix[hmgd.matrix > max.count] <- max.count
h1.matrix[h1.matrix > max.count] <- max.count
dnase.matrix[dnase.matrix > max.count] <- max.count

tss.dist <- dhs.tab$TSS.DIST


pos <- seq(-1000, 1000)

in.region <- (pos > -200) & (pos < 200)
dnase.rate <- apply(dnase.matrix[,in.region], 1, sum)*1e6 /
  (as.numeric(sum(in.region)) * as.numeric(dnase.ttl))

h1.rate <- apply(h1.matrix[,in.region], 1, sum)*1e6 /
  (as.numeric(sum(in.region)) * as.numeric(h1.ttl))

hmgd.rate <- apply(hmgd.matrix[,in.region], 1, sum)*1e6 /
  (as.numeric(sum(in.region)) * as.numeric(hmgd.ttl))


tss.label <- factor(rep("near", length(dnase.rate)),
                    levels=c("near", "med", "far"))

tss.label[tss.dist < 100] <- "near"
tss.label[tss.dist >= 100] <- "med"
tss.label[tss.dist >= 2000] <- "far"



par(mfrow=c(3,1))

for(label in c("near", "med", "far")) {
  f <- tss.label == label
  
  mean.h1.mids <- apply(h1.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(h1.ttl))

  mean.hmgd.mids <- apply(hmgd.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(hmgd.ttl))
  
  mean.mnase.mids <- apply(mnase.matrix[f,], 2, sum)*1e6 /
    (as.numeric(sum(f))*as.numeric(mnase.ttl))
  
  
  mean.dnase.mids <- apply(dnase.matrix[f,], 2, sum)*1e6 /
  (as.numeric(sum(f))*as.numeric(dnase.ttl))
  
  
  plot(pos, mean.dnase.mids, ylim=range(mean.dnase.mids,0),
       type="l", col=color.pal[4])

  points(pos, mean.mnase.mids, type="l", col="black")
  
  points(pos, mean.hmgd.mids, type="l", col=color.pal[1])
  
  points(pos, mean.h1.mids, type="l", col=color.pal[2])
}








cor.test(dnase.rate, h1.rate)

cor.test(dnase.rate, hmgd.rate)


dat <- data.frame(dnase=log2(dnase.rate+1),
                  h1=log2(h1.rate+1),
                  hmgd=log2(hmgd.rate+1),
                  hmgd.h1.ratio=log2((hmgd.rate+1)/(h1.rate+1)),
                  tss.dist=tss.dist,
                  tss.label=tss.label)



library(ggplot2)

ggplot(dat, aes(x=dnase, y=hmgd.h1.ratio, color=tss.label)) +
  geom_point(shape=19, alpha=1/10)


cor.test(dat$dnase[near.tss], dat$hmgd[near.tss])
cor.test(dat$dnase[!near.tss], dat$hmgd[!near.tss])

plot(log2(dnase.rate+1), log2((hmgd.rate+1) / (h1.rate+1)))

