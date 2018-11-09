#!/usr/bin/env Rscript
######################################
## OPTIMIR: Graphs pie Charts
######################################
## Draw scatter plots for each polymiRs representing Allele Specific Expression
## Take path to isomiR_dist.annot file as argument
## Rscript plot_ASE.R /path/to/Results/polymiRs_table.annot
library(gridExtra)
library(ggplot2)

args = commandArgs(trailingOnly=T)
polymiRs_table_path <- args[1]
res_directory <- paste(dirname(polymiRs_table_path), "graphs", sep="/")
if (!file.exists(res_directory)){
    dir.create(res_directory)
}
out_basename <- paste(res_directory, gsub(".annot", "", polymiRs_table_path),sep="/")
table = read.table(polymiRs_table_path, header=T, sep="\t")
table$polymiR_id <- gsub(",", "_", paste(table$polymiR, table$rsID_list, sep="_"))
sub_table = table[which(table$Genotype_list == "0/1" | table$Genotype_list == "1/0"),]

plot_ASE = function(sub_table, title, poly_name) {
    max_scale = max(sub_table$Counts_Alternative, sub_table$Counts_Reference)
    ggplot(sub_table, aes(x=Counts_Alternative, y=Counts_Reference)) +
        ## annotate("rect", xmin=-0.1, xmax=5, ymin=-0.1, ymax=5, fill="red", alpha=0.1) + 
        geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=0.5) +
        geom_point(shape=0, size=0.4) +
        ## facet_wrap(~polymiR_id, scales="free") +
        ## labs(x="Reads aligned on sequence with alternative allele (log scale)", y="Reads aligned on sequence with reference allele (log scale)") +
        labs(x="Reads aligned on sequence with alternative allele", y="Reads aligned on sequence with reference allele") +
        ggtitle(title) +
        theme(
            plot.title = element_text(face="bold", size = 6, hjust=0.5),
            axis.line = element_line(colour = "darkgrey", size = 0.5, linetype = "solid"),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(face="bold", angle = 50, size = 4),
            axis.text.y = element_text(face="bold", angle = 0, size = 4),
            axis.title.x = element_text(face="bold", size = 6),
            axis.title.y = element_text(face="bold", angle = 90, size = 6),
            ) +
        ## scale_x_continuous(breaks=c(0,1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000,2000,3000,4000), trans="log1p", limits=c(-0.01, max(sub_table$Counts_Alternative) + 100), expand=c(0,0)) +
        ## scale_y_continuous(breaks=c(0,1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000,2000,3000,4000), trans="log1p", limits=c(-0.01, max(sub_table$Counts_Reference) + 100), expand=c(0,0))
        scale_x_continuous(limits=c(-(max_scale/100), max_scale + max_scale / 10), expand=c(0,0)) +
        scale_y_continuous(limits=c(-(max_scale/100), max_scale + max_scale / 10), expand=c(0,0))
    ggsave(paste(out_basename, l, "png", sep="."), dpi = 500, width = 6, height = 6, unit="in")
}

sub_table$polymiR_id <- as.factor(gsub(",", "_", paste(sub_table$polymiR, sub_table$rsID_list, sep="_")))
out_basename <- "graphs/ASE"
for (l in levels(sub_table$polymiR_id)){
    subsub_table = sub_table[which(sub_table$polymiR_id == l),]
    title = paste(gsub("_", " / ", l), " (", nrow(subsub_table), " samples)", sep="")
    plot_ASE(subsub_table, title, l)
}


