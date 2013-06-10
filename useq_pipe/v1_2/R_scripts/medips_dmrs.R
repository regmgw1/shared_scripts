#!/usr/local/bin/Rscript

# script runs the various QC steps for the medusa pipeline - utilises the MEDIPS bioconductor package.
args <- commandArgs();
treat<-args[length(args)-7]
control<-args[length(args)-6]
species<-args[length(args)-5]
readDepth<-as.numeric(args[length(args)-4])
dmrSize<-as.numeric(args[length(args)-3])
pvalue<-as.numeric(args[length(args)-2])
path2input<-args[length(args)-1]
path2output<-args[length(args)]

#treat<-"HN105,HN125,HN39"
#control<-"HN29,HN32,HN96"
#species<-"Human"
#readDepth<-10
#dmrSize<-500
#path2output<-"/medical_genomics/tmpFiles/medipsTest_160513/"

library(MEDIPS)
if (species == "Mouse"){
	library(BSgenome.Mmusculus.UCSC.mm9)
	genome <- "BSgenome.Mmusculus.UCSC.mm9"
} else if (species == "Human"){
	library(BSgenome.Hsapiens.UCSC.hg19)
	genome <- "BSgenome.Hsapiens.UCSC.hg19"
} else if (species == "Dog"){
	library(BSgenome.Cfamiliaris.UCSC.canFam2)
	genome <- "BSgenome.Cfamiliaris.UCSC.canFam2"
} else if (species == "Chimp"){
	library(BSgenome.Ptroglodytes.UCSC.panTro2)
	genome<-"BSgenome.Ptroglodytes.UCSC.panTro2"
} else if (species == "Macaca"){
	library(BSgenome.Mmulatta.NCBI.mmul1)
	genome<-"BSgenome.Mmulatta.NCBI.mmul1"
} else {
        stop("Species must be Mouse, Human, Chimp, Macaca or Dog. Currently using mouse mm9, human hg19, chimp panTro2, macaca mmul1 and dog canFam2")
}

setwd(path2output)

treat.v<-as.character(do.call('rbind', strsplit(as.character(treat), ',', fixed=TRUE)))
control.v<-as.character(do.call('rbind', strsplit(as.character(control), ',', fixed=TRUE)))

treat.set<-MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",treat.v[1],".bed",sep=""),window_size=dmrSize)
control.set<-MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",control.v[1],".bed",sep=""),window_size=dmrSize)

for(i in 2:length(treat.v))
{
treat.set<-c(treat.set,MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",treat.v[i],".bed",sep=""),window_size=dmrSize))
}

for(i in 2:length(control.v))
{
control.set<-c(control.set,MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",control.v[i],".bed",sep=""),window_size=dmrSize))
}

cor.matrix<-MEDIPS.correlation(MSets=c(treat.set,control.set),plot=F)
write.table(cor.matrix,file="correlation_matrix.txt",append=F,quote=F,sep="\t",row.names=T,col.names=T)

CS<-MEDIPS.couplingVector(pattern="CG",refObj=treat.set[[1]])
mr.edgeR<-MEDIPS.meth(MSet1=treat.set,MSet2=control.set,CSet=CS,p.adj="bonferroni",diff.method="edgeR",prob.method="poisson",MeDIP=T,CNV=F,minRowSum=readDepth)
mr.edgeR.s<-MEDIPS.selectSig(results=mr.edgeR,p.value=pvalue,adj=T,ratio=NULL,bg.counts=NULL,CNV=F)
mr.edgeR.s.treatLoss<-mr.edgeR.s[which(mr.edgeR.s[,grep("logFC",colnames(mr.edgeR.s))] <0),]
mr.edgeR.s.treatGain<-mr.edgeR.s[which(mr.edgeR.s[,grep("logFC",colnames(mr.edgeR.s))] >0),]
mr.edgeR.s.treatLoss.m<-MEDIPS.mergeFrames(frames=mr.edgeR.s.treatLoss,distance=0)
mr.edgeR.s.treatGain.m<-MEDIPS.mergeFrames(frames=mr.edgeR.s.treatGain,distance=0)
write.table(mr.edgeR.s.treatGain,file=paste("medips_dmrs_hyperTreatment_p",pvalue,".txt",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=T)
write.table(mr.edgeR.s.treatLoss,file=paste("medips_dmrs_hypoTreatment_p",pvalue,".txt",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=T)
dmrFiles<-c(paste("medips_dmrs_hyperTreatment_p",pvalue,".txt",sep=""),paste("medips_dmrs_hypoTreatment_p",pvalue,".txt",sep=""))
write.table(dmrFiles,file="dmrFiles.txt",,append=F,quote=F,sep="\t",row.names=F,col.names=F)
mr.edgeR.s.treatGain.m<-data.frame(mr.edgeR.s.treatGain.m$chr,as.numeric(mr.edgeR.s.treatGain.m$start),as.numeric(mr.edgeR.s.treatGain.m$stop),mr.edgeR.s.treatGain.m$ID)
mr.edgeR.s.treatLoss.m<-data.frame(mr.edgeR.s.treatLoss.m$chr,as.numeric(mr.edgeR.s.treatLoss.m$start),as.numeric(mr.edgeR.s.treatLoss.m$stop),mr.edgeR.s.treatLoss.m$ID)
write.table(mr.edgeR.s.treatLoss.m,file=paste("medips_dmrs_hypoTreatment_p",pvalue,".bed",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=F)
write.table(mr.edgeR.s.treatGain.m,file=paste("medips_dmrs_hyperTreatment_p",pvalue,".bed",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=F)
