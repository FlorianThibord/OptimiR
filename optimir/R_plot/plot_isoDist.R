#!/usr/bin/env Rscript
######################################
## OPTIMIR: Graphs pie Charts
######################################
## Draw pie charts representing isomiRs distribution (on 5'end, 3'end, and global:canonical/not_canonical)
## Take path to isomiR_dist.annot file as argument
## Rscript plot_isoDist.R /path/to/Results/isomiRs_dist.annot

args = commandArgs(trailingOnly=T)
isomiR_table_path <- args[1]
res_directory <- paste(dirname(isomiR_table_path), "graphs", sep="/")
if (!file.exists(res_directory)){
    dir.create(res_directory)
}
out_basename <- paste(res_directory, gsub(".annot", "", isomiR_table_path),sep="/")
iso_d = read.table(isomiR_table_path, header=T)

iso_plot <- function(l, out_basename){
    entry <- l$Reference
    total <- l$Total
    
    colors <- c("grey90", "grey50", "white", "black", "grey30") ## B&W
    colors <- c("olivedrab2","sienna2","dodgerblue4","lightgoldenrod","mediumaquamarine") ## crazycols
    main_title_cex = 2
    legend_cex = 1.5
    title_cex = 1.8
    labels_cex = 1.7
    png(paste0(out_basename, "_", entry, ".png"), width = 9000, height = 3000, units = 'px', res = 500)
    par(mfrow=c(1,3), oma=c(0,0,4,0), mar=c(0,7,3,12), xpd=T)
    ## Cheatsheet for mar and oma: c(b,l,t,r)
    
    ## Graph global dist
    cano_pct <- l$Canonical_PCT
    iso_pct <- 100. - cano_pct
    is_canonical_pct = c(cano_pct, iso_pct)
    is_canonical_labels = c("Canonicals", "IsomiRs")
    pielabels<- paste(round(is_canonical_pct,1), "%", sep="")
    pie(is_canonical_pct, main="Canonicals vs isomiRs \nproportions" , col=colors, radius = 1, labels=pielabels, cex=labels_cex, border="grey30", cex.main = title_cex)
    legend("topright", inset=c(-0.3,0), is_canonical_labels, cex=legend_cex, fill=colors)
    
    ## Graph 3' dist
    isomiRs_3_pct = c(l$End3_Canonical_PCT, l$End3_Trim_PCT, l$End3_Tail_NT_PCT, l$End3_TrimTail_PCT, l$End3_Tail_TE_PCT)
    isomiRs_labels_3 = c("3’ Canonical", "3’ Trim", "3’ NTA", "3’ Trim+NTA", "3’ TA")
    pielabels<- paste(round(isomiRs_3_pct, 1), "%", sep="")
    pie(isomiRs_3_pct, main="3' end isoform distribution", radius=1, col=colors, labels=pielabels, cex=labels_cex, border="grey30", cex.main = title_cex)
    legend("topright", inset=c(-0.3,0), isomiRs_labels_3, cex=legend_cex, fill=colors)

    ## Graph 5' dist
    isomiRs_5_pct = c(l$End5_Canonical_PCT, l$End5_Trim_PCT, l$End5_Tail_NT_PCT, l$End5_TrimTail_PCT, l$End5_Tail_TE_PCT)
    isomiRs_labels_5 = c("5’ Canonical", "5’ Trim", "5’ NTA", "5’ Trim+NTA", "5’ TA")
    pielabels<- paste(round(isomiRs_5_pct,1), "%", sep="")
    pie(isomiRs_5_pct, main="5' end isoform distribution", radius=1, col=colors, labels=pielabels, cex=labels_cex, border="grey30", cex.main = title_cex)
    legend("topright", inset=c(-0.3,0), isomiRs_labels_5, cex=legend_cex, fill=colors)

    title(main=paste(entry, "  -  Read counts = ", total, sep=""), outer=T, cex.main = main_title_cex)
    dev.off()
}

## Plot them ALL:
## by(iso_d, 1:nrow(iso_d), function(row) iso_plot(row, out_basename))
## Just plot TOTAL
iso_plot(iso_d[which(iso_d$Reference == "TOTAL"),], out_basename)
