#!/usr/local/bin/Rscript

# script runs the various QC steps for the medusa pipeline - utilises the MEDIPS bioconductor package.
args <- commandArgs();
treat<-args[length(args)-8]
control<-args[length(args)-7]
species<-args[length(args)-6]
refName<-args[length(args)-5]
readDepth<-as.numeric(args[length(args)-4])
dmrSize<-as.numeric(args[length(args)-3])
pvalue<-as.numeric(args[length(args)-2])
path2input<-args[length(args)-1]
path2output<-args[length(args)]

#treat<-"AST7,AST8,AST9"
#control<-"NSC1,NSC2,NSC3"
#species<-"Human"
#readDepth<-10
#dmrSize<-500
#path2output<-"/medical_genomics/tmpFiles/medipsTest_160513/"

library(MEDIPS)
library(gtools)
if (species == "Mouse"){
	bsGenome<-paste("BSgenome.Mmusculus.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Human"){
	bsGenome<-paste("BSgenome.Hsapiens.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Dog"){
	bsGenome<-paste("BSgenome.Cfamiliaris.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Chimp"){
	bsGenome<-paste("BSgenome.Ptroglodytes.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Macaca"){
	bsGenome<-paste("BSgenome.Mmulatta.NCBI.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else {
        stop("Species must be Mouse, Human, Chimp, Macaca or Dog. Currently using mouse mm9, human hg19, chimp panTro2, macaca mmul1 and dog canFam2")
}

setwd(path2output)

treat.v<-as.character(do.call('rbind', strsplit(as.character(treat), ',', fixed=TRUE)))
control.v<-as.character(do.call('rbind', strsplit(as.character(control), ',', fixed=TRUE)))

treat.set<-MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",treat.v[1],".bed",sep=""),window_size=dmrSize)

fileName<-paste(treat.v[1],".bed",sep="")
GRange.Reads<-getGRange(fileName, path2input, 0, 0, NULL, T)

## Build dataframe of genomic coords
if (length(unique(seqlevels(GRange.Reads))) > 1) 
{
	chromosomes = mixedsort(unique(seqlevels(GRange.Reads)))
}
if (length(unique(seqlevels(GRange.Reads))) == 1)
{
	chromosomes = unique(seqlevels(GRange.Reads))
}
dataset<-get(ls(paste("package:", genome, sep = "")))
chr_lengths<-as.numeric(seqlengths(dataset)[chromosomes])
no_chr_windows<-ceiling(chr_lengths/dmrSize)
supersize_chr<-cumsum(no_chr_windows)
Granges.genomeVec<-MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, dmrSize)
overlap<-countOverlaps(Granges.genomeVec, GRange.Reads)
genomicCounts.df<-data.frame(as.data.frame(Granges.genomeVec),overlap)
colnames(genomicCounts.df)[dim(genomicCounts.df)[2]]<-treat.v[1]

for(i in 2:length(treat.v))
{
treat.set<-c(treat.set,MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",treat.v[i],".bed",sep=""),window_size=dmrSize))
fileName<-paste(treat.v[i],".bed",sep="")
GRange.Reads<-getGRange(fileName, path2input, 0, 0, NULL, T)
Granges.genomeVec<-MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, dmrSize)
overlap<-countOverlaps(Granges.genomeVec, GRange.Reads)
genomicCounts.df<-data.frame(genomicCounts.df,overlap)
colnames(genomicCounts.df)[dim(genomicCounts.df)[2]]<-treat.v[i]
}

control.set<-MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",control.v[1],".bed",sep=""),window_size=dmrSize)
fileName<-paste(control.v[1],".bed",sep="")
GRange.Reads<-getGRange(fileName, path2input, 0, 0, NULL, T)
Granges.genomeVec<-MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, dmrSize)
overlap<-countOverlaps(Granges.genomeVec, GRange.Reads)
genomicCounts.df<-data.frame(genomicCounts.df,overlap)
colnames(genomicCounts.df)[dim(genomicCounts.df)[2]]<-control.v[1]

for(i in 2:length(control.v))
{
control.set<-c(control.set,MEDIPS.createSet(BSgenome=genome,file=paste(path2input,"/",control.v[i],".bed",sep=""),window_size=dmrSize))
fileName<-paste(control.v[i],".bed",sep="")
GRange.Reads<-getGRange(fileName, path2input, 0, 0, NULL, T)
Granges.genomeVec<-MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, dmrSize)
overlap<-countOverlaps(Granges.genomeVec, GRange.Reads)
genomicCounts.df<-data.frame(genomicCounts.df,overlap)
colnames(genomicCounts.df)[dim(genomicCounts.df)[2]]<-control.v[i]
}
write.table(genomicCounts.df,file="readCounts_matrix.txt",append=F,quote=F,sep="\t",row.names=F,col.names=T)
genomicCountsThresh.df<-genomicCounts.df[rowSums(genomicCounts.df[,6:dim(genomicCounts.df)[2]])>=readDepth,]
genomicCountsThresh.bed<-genomicCountsThresh.df[,1:3]
write.table(genomicCountsThresh.bed,file="thresholdWindows.bed",append=F,quote=F,sep="\t",row.names=F,col.names=F)

cor.matrix<-MEDIPS.correlation(MSets=c(treat.set,control.set),plot=F)
write.table(cor.matrix,file="correlation_matrix.txt",append=F,quote=F,sep="\t",row.names=T,col.names=T)

CS<-MEDIPS.couplingVector(pattern="CG",refObj=treat.set[[1]])
mr.edgeR<-MEDIPS.meth(MSet1=treat.set,MSet2=control.set,CSet=CS,p.adj="bonferroni",diff.method="edgeR",prob.method="poisson",MeDIP=T,CNV=F,minRowSum=readDepth)
mr.edgeR.s<-MEDIPS.selectSig(results=mr.edgeR,p.value=pvalue,adj=T,ratio=NULL,bg.counts=NULL,CNV=F)
mr.edgeR.s.treatLoss<-mr.edgeR.s[which(mr.edgeR.s[,grep("logFC",colnames(mr.edgeR.s))] <0),]
mr.edgeR.s.treatGain<-mr.edgeR.s[which(mr.edgeR.s[,grep("logFC",colnames(mr.edgeR.s))] >0),]

dmrFiles<-character()
if (dim(mr.edgeR.s.treatLoss)[1] > 0) {
write.table(mr.edgeR.s.treatLoss,file=paste("medips_dmrs_hypoTreatment_p",pvalue,".txt",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=T)
mr.edgeR.s.treatLoss.m<-MEDIPS.mergeFrames(frames=mr.edgeR.s.treatLoss,distance=0)
mr.edgeR.s.treatLoss.m<-data.frame(mr.edgeR.s.treatLoss.m$chr,as.numeric(mr.edgeR.s.treatLoss.m$start),as.numeric(mr.edgeR.s.treatLoss.m$stop),mr.edgeR.s.treatLoss.m$ID)
write.table(mr.edgeR.s.treatLoss.m,file=paste("medips_dmrs_hypoTreatment_p",pvalue,".bed",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=F)
dmrFiles<-c(dmrFiles,paste("medips_dmrs_hypoTreatment_p",pvalue,".txt",sep=""))
}

if (dim(mr.edgeR.s.treatGain)[1] > 0) {
write.table(mr.edgeR.s.treatGain,file=paste("medips_dmrs_hyperTreatment_p",pvalue,".txt",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=T)
mr.edgeR.s.treatGain.m<-MEDIPS.mergeFrames(frames=mr.edgeR.s.treatGain,distance=0)
mr.edgeR.s.treatGain.m<-data.frame(mr.edgeR.s.treatGain.m$chr,as.numeric(mr.edgeR.s.treatGain.m$start),as.numeric(mr.edgeR.s.treatGain.m$stop),mr.edgeR.s.treatGain.m$ID)
write.table(mr.edgeR.s.treatGain.m,file=paste("medips_dmrs_hyperTreatment_p",pvalue,".bed",sep=""),append=F,quote=F,sep="\t",row.names=F,col.names=F)
dmrFiles<-c(dmrFiles,paste("medips_dmrs_hyperTreatment_p",pvalue,".txt",sep=""))
}
write.table(dmrFiles,file="dmrFiles.txt",,append=F,quote=F,sep="\t",row.names=F,col.names=F)

