#!/usr/local/bin/Rscript

# script generates volcano plots for visualising the dmrs significance - utilises the ggplot2 package. highlights those with fdr < 5%
args <- commandArgs();
path2annotations<-args[length(args)]

setwd(path2annotations)
enr<-dir("./",pattern="EnrichedRegions_BinaryWindowData_FDR7.*peak_cg.txt")
enr<-read.table(enr,head=T,sep="\t")
red<-dir("./",pattern="ReducedRegions_BinaryWindowData_FDR7.*peak_cg.txt")
red<-read.table(red,head=T,sep="\t")
dmr<-rbind(red,enr)
dmr$threshold = as.factor(abs(dmr$BH_Log2) >= 1 & dmr$BH_FDR >= 13)
require(ggplot2)
g = ggplot(data=dmr, aes(x=BH_Log2, y=BH_FDR, colour=threshold))+geom_point(alpha=0.8, size=1.75)+opts(legend.position = "none")+xlab("log2 fold change") + ylab("FDR")
jpeg("dmr_volcano.jpg")
g
dev.off()
