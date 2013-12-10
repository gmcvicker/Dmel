library(RColorBrewer)
source("~/proj/script/R/lib/util.R")


color.pal <-  brewer.pal(9, "Set1")

data.dir <- "~/data/Dmel/MNase/tss_mnase_midpoints_by_gene"

read.totals.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt")

gene.tab <- read.table(p(data.dir, "/gene_summary.fb_gene_names.txt.gz"), header=T)


expr.pct <- gene.tab[["EXPR.PCT.S2_R."]]

expr.filters <-
  data.frame(high.expr = !is.na(expr.pct) & (expr.pct >= 0.75),
             mid.expr = !is.na(expr.pct) & (expr.pct >= 0.25) &
                         (expr.pct <= 0.75),
             low.expr = !is.na(expr.pct) & (expr.pct < 0.25))



data.types <- c("S2", "S2_in_vitro", "S2_in_vitro_031212", "S2_in_vitro_080110",
                "leg", "haltere", "eye", "antenna")
                

flanking <- 1000
n.site <- flanking + flanking + 1
df <- data.frame(POS=seq(-flanking, flanking))

for(type in data.types) {
  count.matrix <- read.matrix(p(data.dir, "/", type, ".txt.gz"))
  ttl <- read.totals.tab$V2[read.totals.tab$V1 == type]

  max.count <- 20
  count.matrix[mnase.matrix > max.count] <- max.count
  
  for(i in 1:ncol(expr.filters)) {
    f <- expr.filters[,i]
    
    mean.fpkm <- apply(count.matrix[f,], 2, sum)*1e9 / (as.numeric(sum(f))*as.numeric(ttl))

    label <- p(type, ".", names(expr.filters)[i], ".FPKM")

    df[label] <- mean.fpkm
    cat(label, "\n")
  }
  
}



pdf("tss_mnase_mids.pdf", width=6, height=10)

par(mfrow=c(3,1))

win.sz <- 10

for(j in 1:ncol(expr.filters)) {

  expr.type <- names(expr.filters)[j]
  
  plot(type="n", x=c(0), y=c(0), ylim=range(0, df[,2:ncol(df)]), xlim=range(df$POS),
       xlab="position relative to TSS (bp)", ylab="FPKM", las=1,
       main=expr.type)

  legend("topleft", legend=data.types,
         lty=1, col=color.pal[1:length(data.types)])
  
  for(i in 1:length(data.types)) {
    type <- data.types[i]
    col <- color.pal[i]

    label <- p(type, ".", expr.type, ".FPKM")    
    cat(label, "\n")

    vals <- filter(df[,"label"], win.sz) / win.sz
    
    points(df[,"POS"], vals, lwd=1, type="l", col=color.pal[i])
  }
}

dev.off()



