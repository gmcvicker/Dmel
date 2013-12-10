

library(RColorBrewer)
pal <- brewer.pal(9, "Set1")

antenna.color <- pal[1]
eye.color <- pal[2]
haltere.color <- pal[3]
leg.color <- pal[4]


tfbs.tab <- read.table("~/data/Dmel/redfly/redfly_tfbs_mnase_midpoints.txt",
                       header=T)



#
# Identify TFs that have most significant diffs between lineages using
# chisq test
#

tf.names <- unique(tfbs.tab$TF.NAME)


# total number of mapped reads
antenna.total <- 8516055
eye.total <- 4205727
haltere.total <- 9757530
leg.total <- 12055247

combined.total <- antenna.total + eye.total + haltere.total + leg.total

tf.pval <- rep(NA, length(tf.names))
tf.n.tfbs <- rep(NA, length(tf.names))

for(i in 1:length(tf.names)) {
  tf <- tf.names[i]
  idx <- which(tfbs.tab$TF.NAME == tf)

  tf.n.tfbs[i] <- length(idx)
  
  eye.count <- tfbs.tab$eye.MNASE.COUNT[idx]
  haltere.count <- tfbs.tab$haltere.MNASE.COUNT[idx]
  antenna.count <- tfbs.tab$antenna.MNASE.COUNT[idx]
  leg.count <- tfbs.tab$leg.MNASE.COUNT[idx]

  combined.count <- eye.count + haltere.count + antenna.count + leg.count
  
  eye.rate <- eye.count / eye.total
  haltere.rate <- haltere.count / haltere.total
  antenna.rate <- antenna.count / antenna.total
  leg.rate <- leg.count / leg.total

  combined.rate <- combined.count / combined.total

  
  model <- lm(eye.rate ~ haltere.rate + antenna.rate + leg.rate)

  m <- matrix(c(sum(eye.count), sum(haltere.count),
                sum(antenna.count), sum(leg.count),
                eye.total, haltere.total, antenna.total, leg.total),
              dimnames=list(c("TFBS counts", "total counts"),
                c("eye", "haltere", "antenna", "leg")),
              ncol=4, nrow=2, byrow=T)

  stat <- chisq.test(m)
  
  tf.pval[i] <- stat$p.value
}


hist(tf.n.tfbs)


f <- tf.n.tfbs >= 5
sum(f)




#
# Make boxplots of relative MNase rates for most significant TFs
#

profile.tab <- read.table("~/data/Dmel/redfly/redfly_tfbs_mnase_profile.txt",
                          header=F)
n.flanking <- 1000


filename <- "mnase_lineage_rates.pdf"

pdf(filename, width=10, height=5)

par(mfrow=c(1,2))

f <- (tf.n.tfbs >= 5) & (tf.pval < 1e-4)

pseudo.count <- 0.5

for(i in 1:length(tf.names[f])) {  
  tf <- tf.names[f][i]

  ########## make profile plots
  
  idx <- which((profile.tab$V1 == tf) & (profile.tab$V2 == "eye"))
  n.sites <- length(idx)  
  eye.counts <- as.matrix(profile.tab[idx, 3:ncol(profile.tab)])
  
  idx <- which((profile.tab$V1 == tf) & (profile.tab$V2 == "haltere"))
  haltere.counts <- as.matrix(profile.tab[idx, 3:ncol(profile.tab)])

  idx <- which((profile.tab$V1 == tf) & (profile.tab$V2 == "antenna"))
  antenna.counts <- as.matrix(profile.tab[idx, 3:ncol(profile.tab)])

  idx <- which((profile.tab$V1 == tf) & (profile.tab$V2 == "leg"))
  leg.counts <- as.matrix(profile.tab[idx, 3:ncol(profile.tab)])

  eye.rate <- apply(eye.counts, 2, sum) / (eye.total * n.sites)
  haltere.rate <- apply(haltere.counts, 2, sum) / (haltere.total * n.sites)
  antenna.rate <- apply(antenna.counts, 2, sum) / (antenna.total * n.sites)
  leg.rate <- apply(leg.counts, 2, sum) / (leg.total * n.sites)
  
  win <- rep(1, 150)
  smooth.eye.rate <- filter(eye.rate, win) / sum(win)             
  smooth.haltere.rate <- filter(haltere.rate, win) / sum(win)
  smooth.antenna.rate <- filter(antenna.rate, win) / sum(win)
  smooth.leg.rate <- filter(leg.rate, win) / sum(win)
  
  pos <- 1:ncol(eye.counts) - n.flanking
  
  plot(pos, smooth.eye.rate * 1e9, type="l",
       ylim=range(0, 22), col=eye.color,
       main=tf, ylab="MNase midpoints per billion mapped reads",
       xlab="distance from TFBS (bp)")
    
  points(pos, smooth.haltere.rate*1e9, type="l", col=haltere.color)
  points(pos, smooth.antenna.rate*1e9, type="l", col=antenna.color)
  points(pos, smooth.leg.rate*1e9, type="l", col=leg.color)


  ##################### make boxplots

  idx <- which(tfbs.tab$TF.NAME == tf)

  n <- length(idx)
  
  eye.count <- tfbs.tab$eye.MNASE.COUNT[idx] + pseudo.count
  haltere.count <- tfbs.tab$haltere.MNASE.COUNT[idx] + pseudo.count
  antenna.count <- tfbs.tab$antenna.MNASE.COUNT[idx] + pseudo.count
  leg.count <- tfbs.tab$leg.MNASE.COUNT[idx] + pseudo.count

  combined.count <- eye.count + haltere.count + antenna.count + leg.count
  
  eye.rate <- eye.count / eye.total
  haltere.rate <- haltere.count  / haltere.total
  antenna.rate <- antenna.count / antenna.total
  leg.rate <- leg.count / leg.total

  combined.rate <- combined.count / combined.total

  eye.rel.rate <- log2(eye.rate / combined.rate)
  haltere.rel.rate <- log2(haltere.rate / combined.rate)
  antenna.rel.rate <- log2(antenna.rate / combined.rate)
  leg.rel.rate <- log2(leg.rate / combined.rate)
  
  boxplot(cbind(antenna.rel.rate, eye.rel.rate, haltere.rel.rate,
                leg.rel.rate),
          col=c(antenna.color, eye.color, haltere.color, leg.color),
          ylab="relative nucleosome occupancy",
          main=paste(tf, "\nn sites=", n,
            ", p=", signif(tf.pval[f][i],3), sep=""),
          names=c("antenna", "eye", "haltere", "leg"))

}

dev.off()






##############################


gene.expr.tab <- read.table("~/data/Dmel/mod_encode/rna_seq/gene_rna_seq_expr.txt", header=T)


tf.expr.tab <- data.frame(TF.NAME=tf.names,
                          EYE.EXPR=rep(NA, length(tf)),
                          ANTENNA.EXPR=rep(NA, length(tf)),
                          HALTERE.EXPR=rep(NA, length(tf)),
                          LEG.EXPR=rep(NA, length(tf)))

# make table of TF gene expression

for(i in 1:length(tf.names)) {
  tf <- tf.names[i]  

  gene.idx <- which(gene.expr.tab$GENE.SYMBOL == as.character(tf))

  if(length(gene.idx) == 1) {
    
    expr <- c(gene.expr.tab[["RNASEQ.RPKM.ML_DmD11"]][gene.idx],
              gene.expr.tab[["RNASEQ.RPKM.ML_DmD20_c5"]][gene.idx],
              gene.expr.tab[["RNASEQ.RPKM.ML_DmD17_c3"]][gene.idx],
              gene.expr.tab[["RNASEQ.RPKM.CME_L1"]][gene.idx])

    tf.expr.tab[i,2:ncol(tf.expr.tab)] <- expr
  }  
}


tf.expr.sd <- apply(tf.expr.tab[,2:ncol(tf.expr.tab)], 1, sd)
tf.expr.mean <- apply(tf.expr.tab[,2:ncol(tf.expr.tab)], 1, mean)
tf.expr.coef.var <- tf.expr.sd / tf.expr.mean

# only consider TFs that have some expression variation
f <- (tf.n.tfbs >= 5) & !is.na(tf.expr.sd) & (tf.expr.sd > 1.0)

pseudo.count <- 0.5

expr.cor <- rep(NA, sum(f))
cor.p.val <- rep(NA, sum(f))

lm.beta <- rep(NA, sum(f))
lm.intercept <- rep(NA, sum(f))
lm.p.val <- rep(NA, sum(f))


pdf("dmel_gene_expr_mnase.pdf", width=5, height=5)

for(i in 1:length(tf.names[f])) {
  tf <- tf.names[f][i]

  idx <- which(tfbs.tab$TF.NAME == tf)

  n <- length(idx)
  
  eye.count <- tfbs.tab$eye.MNASE.COUNT[idx] + pseudo.count
  haltere.count <- tfbs.tab$haltere.MNASE.COUNT[idx] + pseudo.count
  antenna.count <- tfbs.tab$antenna.MNASE.COUNT[idx] + pseudo.count
  leg.count <- tfbs.tab$leg.MNASE.COUNT[idx] + pseudo.count

  combined.count <- eye.count + haltere.count + antenna.count + leg.count
  
  eye.rate <- eye.count / eye.total
  haltere.rate <- haltere.count  / haltere.total
  antenna.rate <- antenna.count / antenna.total
  leg.rate <- leg.count / leg.total

  combined.rate <- combined.count / combined.total

  eye.rel.rate <- log2(eye.rate / combined.rate)
  haltere.rel.rate <- log2(haltere.rate / combined.rate)
  antenna.rel.rate <- log2(antenna.rate / combined.rate)
  leg.rel.rate <- log2(leg.rate / combined.rate)
  
  expr.row <- tf.expr.tab[f,][i,2:ncol(tf.expr.tab)]

  expr <- c(rep(expr.row$EYE.EXPR, n),
                 rep(expr.row$ANTENNA.EXPR, n),
                 rep(expr.row$HALTERE.EXPR, n),
                 rep(expr.row$LEG.EXPR, n))
                  
  colors <- c(rep(eye.color, n), rep(antenna.color, n),
              rep(haltere.color, n), rep(leg.color, n))

  mnase.rate <- c(eye.rel.rate, antenna.rel.rate,
                      haltere.rel.rate, leg.rel.rate)

  #mnase.rate <- c(eye.rate, antenna.rate,
  #                haltere.rate, leg.rate)

  
  mdl <- lm(mnase.rate ~ expr)
  lm.intercept[i] <- mdl$coef[1]
  lm.beta[i] <- mdl$coef[2]
  smry <- summary(mdl)
  lm.p.val[i] <- pf(smry$fstat[1], df1=smry$df[1], df2=smry$df[2], lower.tail=F)
  
  stat <- cor.test(expr, mnase.rate, method="pearson")
  r.val <- stat$estimate
  p.val <- stat$p.value
  expr.cor[i] <- r.val
  cor.p.val[i] <- p.val

  main <- paste(tf, "\nR=", signif(r.val, 2), ", p=", signif(p.val, 2),
                sep="")
    
  plot(expr, mnase.rate, pch=21, col="black", bg=colors,
       main=main, ylab="relative nucleosome occupancy",
       xlab="gene expression (RPKM)",
       xlim=range(0, expr),
       ylim=range(mnase.rate, 0, -3, 3))

  legend("topright", pch=21,
         col="black",
         pt.bg=c(eye.color, antenna.color, haltere.color, leg.color),
         legend=c("ML-DmD11 (eye)", "ML-DmD20 (antenna)",
           "ML-DmD17 (haltere)",  "CMEL1 (leg)"))
  
}

dev.off()


# bonferroni correction
threshold <- 0.05 / sum(f)

plot(lm.beta, -log10(lm.p.val))
f <- lm.p.val < threshold
points(lm.beta[f], -log10(lm.p.val[f]), col="red", bg="red", pch=21)
lines(x=c(0,0), y=c(0, 10), col="grey50")

       
