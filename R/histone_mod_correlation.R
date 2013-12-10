
source("~/proj/script/R/lib/util.R")
source("~/proj/script/R/10_IND/heatmap.R")



draw.mark.heatmap <- function(mark.matrix, mark.names,
                              labels=NA, scale=0.8) {
  n.mark <- length(mark.names)

  # make correlation matrix
  cor.matrix <- matrix(nrow=n.mark, ncol=n.mark,
                       dimnames=list(mark.names, mark.names))

  for(i in 1:n.mark) {
    mark1 <- mark.names[i]
    
    for(j in 1:i) {
      mark2 <- mark.names[j]
      r <- cor(mark.matrix[[mark1]], mark.matrix[[mark2]])
      cor.matrix[i,j] <- r
      cor.matrix[j,i] <- r
    }
  }

  dist <- as.dist(1-cor.matrix)

  clust <- hclust(dist)

  draw.heatmap(cor.matrix[clust$order,clust$order], scale=0.8,
               labels=labels[clust$order])

  return(clust)
}



tab <- read.big.table("~/data/Dmel/histone_mod/mod_encode_histone_mod_vals.2kb.TSS.txt", header=T)

hmgd.h1.tab <- read.big.table("~/data/Dmel/histone_mod/hmgd_h1_vals.2kb.TSS.txt", header=T)

# normalize HMGD and H1 using plain MNase-seq
n.h1.rep <- 4
n.hmgd.rep <- 5
hmgd.names <- paste("HMGD_rep", seq(n.hmgd.rep), sep="")
h1.names <- paste("H1_rep", seq(n.h1.rep), sep="")

hmgd.h1.norm <- log2(hmgd.h1.tab[,c(hmgd.names, h1.names)] / hmgd.h1.tab[["mnase"]])

combined.tab <- as.data.frame(cbind(tab, hmgd.h1.norm))

ptm.tab <- read.table("~/data/Dmel/histone_ptms.txt", header=T)

filter <- (combined.tab$N.DEF > 1000) & (hmgd.h1.tab$mnase > 1.0)

mark.names <- colnames(combined.tab)[5:ncol(combined.tab)]
keep <- ptm.tab$KEEP==1


pdf("histone_mod_cor_matrix.TSS.pdf", width=20, height=20)

clust <- draw.mark.heatmap(combined.tab[filter,], mark.names[keep],
                           labels=ptm.tab$MARK.NAME[keep])

dev.off()


pdf("histone_mod_hclust.TSS.pdf", width=20, height=20)
plot(clust)
dev.off()







mark.names <- colnames(hmgd.h1.tab)[4:ncol(hmgd.h1.tab)]
                      
pdf("hmgd_h1_cor_matrix.pdf", width=8, height=8)
clust <- draw.mark.heatmap(hmgd.h1.tab[filter,], mark.names)
dev.off()

pdf("hmgd_h1_hclust.pdf", width=8, height=8)
plot(clust)
dev.off()

mark.names <- names(hmgd.h1.norm[filter,])
                      
pdf("hmgd_h1_norm_cor_matrix.pdf", width=8, height=8)
clust <- draw.mark.heatmap(hmgd.h1.norm[filter,], mark.names)
dev.off()

pdf("hmgd_h1_norm_hclust.pdf", width=8, height=8)
plot(clust)
dev.off()



