#!/usr/local/bin/Rscript

# script runs the various QC steps for the medusa pipeline - utilises the MEDIPS bioconductor package.
args <- commandArgs();
sample<-args[length(args)-4]
species<-args[length(args)-3]
windowSize<-as.numeric(args[length(args)-2])
strands<-as.numeric(args[length(args)-1])
path2output<-args[length(args)]

#sample<-"HN105"
#species<-"Human"
#windowSize<-100
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

sampleWig.set<-MEDIPS.createSet(BSgenome=genome,file=paste(path2output,"/",sample,".bed",sep=""),window_size=windowSize)
MEDIPS.exportWIG(Set=sampleWig.set,file=paste(sample,"_rpm.wig",sep=""),format="rpkm")
if(strands == 0) {
sr<-MEDIPS.saturation(file=paste(path2output,"/",sample,".bed",sep=""),BSgenome=genome,uniq=T,extend=0,shift=0,window_size=windowSize,nit=10,rank=T)
png(file=paste(sample,"_saturation.png",sep=""),type="cairo")
MEDIPS.plotSaturation(sr)
dev.off()
cr<-MEDIPS.seqCoverage(file=paste(path2output,"/",sample,".bed",sep=""),pattern="CG",BSgenome=genome,uniq=T,extend=0,shift=0)
png(file=paste(sample,"_coveragePie.png",sep=""),type="cairo")
MEDIPS.plotSeqCoverage(cr,type="pie",cov.level=c(0,1,5,10,20,50))
dev.off()
er<-MEDIPS.CpGenrich(file=paste(path2output,"/",sample,".bed",sep=""),BSgenome=genome,uniq=T,extend=0,shift=0)
write.table(er,paste(sample,"_enrichment.txt",sep=""),sep="\t",quote=F)
}
