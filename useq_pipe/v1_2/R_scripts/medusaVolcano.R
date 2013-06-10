#!/usr/local/bin/Rscript

# script generates volcano plots for visualising the dmrs significance - utilises the ggplot2 package. highlights those with fdr < 5%
args <- commandArgs();
pvalue<-as.numeric(args[length(args)-2])
path2dmrs<-args[length(args)-1]
path2annotations<-args[length(args)]

setwd(path2annotations)
enr<-read.table(paste(path2dmrs,"/medips_dmrs_hyperTreatment_p",pvalue,".txt",sep=""),head=T,sep="\t")
red<-read.table(paste(path2dmrs,"/medips_dmrs_hypoTreatment_p",pvalue,".txt",sep=""),,head=T,sep="\t")
dmr<-rbind(red,enr)
strictPvalue<-pvalue/10
dmr$threshold = as.factor(abs(dmr$edgeR.logFC) >= 1 & dmr$edgeR.adj.p.value >= strictPvalue)
require(ggplot2)
g = ggplot(data=dmr, aes(x=edgeR.logFC, y=-log10(edgeR.adj.p.value), colour=threshold))+geom_point(alpha=0.8, size=1.75)+opts(legend.position = "none")+xlab("log fold change") + ylab("-log10(adjusted pvalues)")
jpeg("dmr_volcano.jpg")
g
dev.off()
