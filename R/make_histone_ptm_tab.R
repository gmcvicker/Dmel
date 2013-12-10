
source("~/proj/script/R/lib/util.R")

tab <- read.big.table("~/data/Dmel/histone_mod/mod_encode_histone_mod_vals.2kb.TSS.txt", header=T)

hmgd.h1.tab <- read.big.table("~/data/Dmel/histone_mod/hmgd_h1_vals.2kb.TSS.txt", header=T)

# normalize HMGD and H1 using plain MNase-seq
n.h1.rep <- 4
n.hmgd.rep <- 5
hmgd.names <- paste("HMGD_rep", seq(n.hmgd.rep), sep="")
h1.names <- paste("H1_rep", seq(n.h1.rep), sep="")

hmgd.h1.norm <- log2(hmgd.h1.tab[,c(hmgd.names, h1.names)] / hmgd.h1.tab[["mnase"]])

combined.tab <- as.data.frame(cbind(tab, hmgd.h1.norm))


mark.names <- colnames(combined.tab)[5:ncol(combined.tab)]

write.table(mark.names, file="mark_names.txt", quote=F,
            col.names=F, row.names=F)

