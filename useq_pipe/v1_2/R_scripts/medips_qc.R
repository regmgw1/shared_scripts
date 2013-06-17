#!/usr/local/bin/Rscript

# script runs the various QC steps for the medusa pipeline - utilises the MEDIPS bioconductor package.
args <- commandArgs();
sample<-args[length(args)-5]
species<-args[length(args)-4]
refName<-args[length(args)-3]
windowSize<-as.numeric(args[length(args)-2])
strands<-as.numeric(args[length(args)-1])
path2output<-args[length(args)]

library(MEDIPS)
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
