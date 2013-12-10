
calc.rate <- function(tab) {
  total.reads <- sum(as.numeric(tab$READ.COUNT))
  rate.per.site <- tab$READ.COUNT*1e6 / (total.reads * tab$TOTAL.SITES)
  return(rate.per.site)
}

data.types <- c("HMGD", "H1", "MNase_S2")

data.dir <- "~/data/Dmel/reads_by_region_type"

path <- paste(data.dir, "/", data.types, ".txt", sep="")


mnase.tab <- read.table(paste(data.dir, "/MNase_S2.txt", sep=""), header=T)
h1.tab <- read.table(paste(data.dir, "/H1.txt", sep=""), header=T)
hmgd.tab <- read.table(paste(data.dir, "/HMGD.txt", sep=""), header=T)

mnase.tab['RATE'] <- calc.rate(mnase.tab)
h1.tab['RATE'] <- calc.rate(h1.tab)
hmgd.tab['RATE'] <- calc.rate(hmgd.tab)


combined.tab <- rbind(mnase.tab, h1.tab, hmgd.tab)
combined.tab['type'] <- factor(levels=c("MNase", "H1", "HMGD"))

combined.tab['type'] <- c(rep("MNase", nrow(mnase.tab)),
                          rep("H1", nrow(mnase.tab)),
                          rep("HMGD", nrow(mnase.tab)))
                                   

rate.matrix <- as.matrix(data.frame(MNase=combined.tab$RATE[combined.tab$type=="MNase"],
                                    H1=combined.tab$RATE[combined.tab$type=="H1"],
                                    HMGD=combined.tab$RATE[combined.tab$type=="HMGD"]),
                         row.names=tab$REGION.TYPE)

library(RColorBrewer)
color.pal <- brewer.pal(9,"Set1")

pdf("hmgd_h1_rate_by_region.pdf", width=7, height=5)

barplot(t(rate.matrix), beside=T, las=1,
        #col=color.pal[3:5],
        names=tab$REGION.TYPE,
        legend=c("MNase", "H1", "HMGD"),
        ylab="reads per million mapped reads",
        xlab="genomic region")

dev.off()

