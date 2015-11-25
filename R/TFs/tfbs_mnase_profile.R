


cell.types <- c("S2", "S2_in_vitro_031212", "S2_in_vitro_080110",
                "antenna", "eye", "haltere", "leg")

colors <- c("blue", "black", "grey80",
            "red", "orange", "green", "purple")


out.dir <- "~/data/Dmel/TFs"

#tf.filename <- paste(out.dir, "/pwm_score_thresholds.most_variables.0.0001.txt",
#                     sep="")


#tf.filename <- paste(out.dir, "/kmers.txt",
#                     sep="")

tf.filename <- paste(out.dir, "/kmer_outliers.txt",
                     sep="")


tf.tab <- read.table(tf.filename, header=F)

tf.names <- as.character(tf.tab$V1)

for(flank in c(2000, 500, 200)) {
    
    pdf(paste("tf_mnase_profile.", flank, ".pdf", sep=""), width=8, height=5)

    for(tf in tf.names) {

        plot(x=c(), y=c(), xlim=c(-flank, flank), ylim=c(0,0.0004), type="n",
             ylab="MNase midpoint density", xlab="distance from TFBS",
             main=tf)

        for(i in 1:length(cell.types)) {
            color <- colors[i]
            cell.type <- cell.types[i]
            fwd.file <- paste(out.dir, "/", tf, "/fwd_", cell.type,
                              ".txt", sep="")
            tab <- read.table(fwd.file)

            pos <- tab$V1
            mid.density <- tab$V2 / sum(as.numeric(tab$V2))
            
            points(pos, mid.density, col=color, type="l")

            
            legend("bottomright", legend=cell.types, col=colors, lty=1)
        }

    }
    dev.off()
}



##########################


cell.types <- c("S2", "antenna", "eye", "haltere", "leg")

in.vitro.cell.type <- "S2_in_vitro_031212"

colors <- c("blue", "red", "orange", "green", "purple")

out.dir <- "~/data/Dmel/TFs"

tf.filename <- paste(out.dir, "/pwm_score_thresholds.most_variables.0.0001.txt",
                     sep="")

tf.tab <- read.table(tf.filename, header=F)

tf.names <- as.character(tf.tab$V1)

for(flank in c(2000, 500, 200)) {
    
    pdf(paste("tf_mnase_rel_profile.", flank, ".pdf", sep=""), width=8, height=5)


    
    for(tf in tf.names) {

        plot(x=c(), y=c(), xlim=c(-flank, flank), ylim=c(-1.0, 1.0), type="n",
             ylab="log2 MNase midpoint in vivo / in vitro ratio",
             xlab="distance from TFBS",
             main=tf)


        segments(x0=c(-2000), x1=c(2000),
                 y0=c(0), y1=c(0), col="black")

        in.vitro.file <- paste(out.dir, "/", tf, "/fwd_", in.vitro.cell.type,
                               ".txt", sep="")

        in.vitro.tab <- read.table(in.vitro.file)
        
        in.vitro.density <- in.vitro.tab$V2 / sum(as.numeric(in.vitro.tab$V2))

        
        for(i in 1:length(cell.types)) {
            color <- colors[i]
            cell.type <- cell.types[i]
            fwd.file <- paste(out.dir, "/", tf, "/fwd_", cell.type,
                              ".txt", sep="")
            tab <- read.table(fwd.file)
            
            pos <- tab$V1
            mid.density <- tab$V2 / sum(as.numeric(tab$V2))

            rel.density <- log2(mid.density / in.vitro.density)
            
            points(pos, rel.density, col=color, type="l")

            legend("bottomright", legend=cell.types, col=colors, lty=1)
        }
    }
    dev.off()
}

