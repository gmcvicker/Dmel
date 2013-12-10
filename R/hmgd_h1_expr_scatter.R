######################################################################
# Draws HMGD and H1 vs expression scatter plots used for paper figure
######################################################################

source("~/proj/script/R/lib/util.R")

data.dir <- "~/data/Dmel/MNase/tss_mnase_midpoints_by_gene"

read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")

gene.tab <- read.table(p(data.dir, "/gene_summary.fb_gene_names.txt.gz"), header=T)

h1.matrix <- read.matrix(p(data.dir, "/h1.txt.gz"))
h1.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "h1"]

hmgd.matrix <- read.matrix(p(data.dir, "/hmgd.txt.gz"))
hmgd.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "hmgd"]

mnase.matrix <- read.matrix(p(data.dir, "/S2.txt.gz"))
mnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "S2"]

dnase.matrix <- read.matrix(p(data.dir, "/dnase.txt.gz"))
dnase.ttl <- read.totals.tab$V2[read.totals.tab$V1 == "dnase"]

max.count <- 20
mnase.matrix[mnase.matrix > max.count] <- max.count
hmgd.matrix[hmgd.matrix > max.count] <- max.count
h1.matrix[h1.matrix > max.count] <- max.count




#
# Scatter plots of H1 / HMGD ratio
#

# add small pseudocount to avoid zeros when taking log
pseudo <- 0.5 * min(gene.tab["EXPR.S2_R."][gene.tab["EXPR.S2_R."] > 0])
expr.rpkm <- gene.tab[["EXPR.S2_R."]] + pseudo

f <- !is.na(expr.rpkm)


flanking <- 1000
pos <- seq(-flanking, flanking)

region.start <- min(pos)
region.end <- max(pos)

in.region <- (pos >= region.start) & (pos <= region.end)

n.site <- as.numeric(sum(in.region))

region.hmgd.counts <- apply(hmgd.matrix[,in.region], 1, sum)
region.h1.counts <- apply(h1.matrix[,in.region], 1, sum)
region.mnase.counts <- apply(mnase.matrix[,in.region], 1, sum)

pseudo <- 2
hmgd.rpkm <- (region.hmgd.counts+pseudo)*1e9 / ((region.end - region.start) * hmgd.ttl)
h1.rpkm <- (region.h1.counts+pseudo)*1e9 / ((region.end - region.start) * h1.ttl)
mnase.rpkm <- (region.mnase.counts+pseudo)*1e9 / ((region.end - region.start) * mnase.ttl)

hmgd.h1.ratio <- log2(hmgd.rpkm / h1.rpkm)
          
pdf("HMGD1_vs_expr_scatter.pdf", width=6, height=6)

mdl <- lm(log10(hmgd.rpkm[f]) ~ log10(expr.rpkm[f]) )
res <- cor.test(log10(expr.rpkm[f]), log10(hmgd.rpkm[f]))
r <- res$est
plot(log10(expr.rpkm[f]), log10(hmgd.rpkm[f]), col="grey50",
     xlab="S2 RNA-seq expression (log10 RPKM)", ylab="HMGD1 (log10 RPKM)",
     main=paste("HMGD1, R=", signif(r,2), sep=""))

abline(mdl$coef[1], mdl$coef[2], col="red")
          
dev.off()


pdf("H1_vs_expr_scatter.pdf", width=6, height=6)

mdl <- lm(log10(h1.rpkm[f]) ~ log10(expr.rpkm[f]) )
res <- cor.test(log10(expr.rpkm[f]), log10(h1.rpkm[f]))
r <- res$est
plot(log10(expr.rpkm[f]), log10(h1.rpkm[f]), col="grey50",
     xlab="S2 RNA-seq expression (log10 RPKM)", ylab="H1 (log10 RPKM)",
     main=paste("HMGD1, R=", signif(r,2), sep=""))
abline(mdl$coef[1], mdl$coef[2], col="red")
          
dev.off()


pdf("MNase_vs_expr_scatter.pdf", width=6, height=6)
mdl <- lm(log10(mnase.rpkm[f]) ~ log10(expr.rpkm[f]) )
res <- cor.test(log10(expr.rpkm[f]), log10(mnase.rpkm[f]))
r <- res$est
plot(log10(expr.rpkm[f]), log10(mnase.rpkm[f]), col="grey50",
     xlab="S2 RNA-seq expression (log10 RPKM)", ylab="Total nucleosomes (log10 RPKM)",
     main=paste("Total nucleosomes, R=", signif(r,2), sep=""))
abline(mdl$coef[1], mdl$coef[2], col="red")
          
dev.off()

