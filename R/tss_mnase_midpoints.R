

library(RColorBrewer)

color.pal <- brewer.pal(9, "Set1")
color.pal[6] <- "black"


# count.tab <- read.table("~/data/Dmel/MNase/mapped_read_counts.txt",
#                        header=F)

count.tab <- read.table("~/data/Dmel/mapped_read_counts.txt.gz",
                        header=F)

lineages <- c("eye", "haltere", "leg", "antenna", "S2", "S2_in_vitro_031212")
n.lineage <- length(lineages)

smooth.win <- 20


for(expr.type in c("all", "low_expr", "mid_expr", "high_expr")) {

  pdf(paste("~/data/Dmel/tss_mnase_midpoints.", expr.type, ".pdf", sep=""), width=8, height=5)
  
  filename <- paste("~/data/Dmel/tss_mnase_midpoints.",
                    expr.type, ".txt.gz", sep="")

  tab <- read.table(filename, header=T)


  out.filename <- paste("mnase_density_tss_dm3.", expr.type, ".pdf", sep="")
  # pdf(out.filename, width=8, height=5)


  # ylim <- c(0, 0.00160)
  ylim <- range(0, 25)

  
  plot(x=c(0), y=c(0), type="n", ylim=ylim, xlim=c(-1000, 1000),
       ylab="FPKM",
       xlab="distance from TSS (bp)", las=1)

  for(i in 1:n.lineage) {
    lin <- lineages[i]
    counts <- tab[[lin]]

    total.mapped <- count.tab$V2[count.tab$V1 == lin]
    n.site <- tab$N.SITE

    # freq <- counts / sum(counts)
    rate <- (counts / n.site/ total.mapped) * 1e9
    
    smooth <- filter(rate, rep(1, smooth.win)) / smooth.win
    
    lines(tab$POS, smooth, type="l", col=color.pal[i])

                                        # draw positions of peaks
    pos.label <- paste("PEAK.POS.", lin, sep="")
    height.label <- paste("PEAK.HEIGHT.", lin, sep="")
    
                                        # points(peak.tab[[pos.label]], peak.tab[[height.label]],
                                        #       col=color.pal[i])  
  }

  legend("topleft", col=color.pal[1:n.lineage],
         legend=lineages, lwd=rep(1, n.lineage))

  dev.off()
}







###############




# estimate nucleosome repeat length from separation of peaks
region.size <- 180
left.bounds <- seq(20, 1000-region.size, by=region.size)
right.bounds <- left.bounds + region.size

segments(x0=left.bounds, x1=left.bounds,
         y0=0, y1=1)

segments(x0=right.bounds, x1=right.bounds,
         y0=0, y1=1)

n.region <- length(left.bounds)


peak.tab <- data.frame(PEAK.NUM=1:n.region)

for(lin.idx in 1:length(lineages)) {  
    lin <- lineages[lin.idx]
    counts <- tab[[lin]]

    freq <- counts / sum(counts)
    smooth <- filter(freq, rep(1, smooth.win)) / smooth.win

    peak.height <- rep(NA, n.region)
    peak.pos <- rep(NA, n.region) 
                        
    for(i in 1:n.region) {
      f <- (tab$POS >= left.bounds[i]) & (tab$POS < right.bounds[i])
      peak.height[i] <- max(smooth[f])
      peak.pos[i] <- which(smooth[f] >= max(smooth[f]))[1] + left.bounds[i]    
    }
    
    pos.label <- paste("PEAK.POS.", lin, sep="")
    height.label <- paste("PEAK.HEIGHT.", lin, sep="")
    
    peak.tab[[pos.label]] <- peak.pos
    peak.tab[[height.label]] <- peak.height
}

legend("topleft", legend=lineages,
       pch=21, col=color.pal[1:length(lineages)])



pdf("nucleosome_repeat_length_tss_dm3.pdf", width=6, height=6)

nucleosome.num <- 1:n.region

plot(x=c(0), y=c(0), type="n",
     xlim=c(1, 5), ylim=c(150, 1000),
     xlab="nucleosome number",
     ylab="distance from TSS (bp)")

nrl.estimate <- rep(NA, length(lineages))
  
for(lin.idx in 1:length(lineages)) {
  lin <- lineages[lin.idx]
  
  pos.label <- paste("PEAK.POS.", lin, sep="")
    
  peak.pos <- peak.tab[[pos.label]]
  
  mdl <- lm(peak.pos ~ nucleosome.num)
  points(nucleosome.num, peak.pos, col=color.pal[lin.idx])
  abline(a=mdl$coef[1], b=mdl$coef[2], col=color.pal[lin.idx])

  nrl.estimate[lin.idx] <-  mdl$coef[2]
}

legend("topleft", legend=paste(lineages,signif(nrl.estimate, 4), "bp"),
       pch=21, col=color.pal[1:length(lineages)])
       
dev.off()


