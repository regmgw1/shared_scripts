#!/usr/local/bin/Rscript

args <- commandArgs();
input<-args[length(args)-3]
species<-args[length(args)-2]
path2scripts<-args[length(args)-1]
path2output<-args[length(args)]
print (input)

library(MEDIPS)
if (species == "Mouse"){
	library(BSgenome.Mmusculus.UCSC.mm9)
	genome <- "BSgenome.Mmusculus.UCSC.mm9"
} else if (species == "Human"){
	library(BSgenome.Hsapiens.UCSC.hg19)
	genome <- "BSgenome.Hsapiens.UCSC.hg19"
} else {
        stop("Species must be Mouse or Human. Currently using mouse mm9 and human hg19")
}
source(paste(path2scripts,"/medips_strand_function.R",sep=""))
setwd(path2output)
data.set<-MEDIPS.readAlignedSequences(BSgenome=genome,file=paste(path2output,"/",input,".medips",sep=""))
data.set<-MEDIPS.genomeVector(data=data.set, extend=0, bin_size=50)
sr.data<-MEDIPS.saturationAnalysis(data=data.set, bin_size=50, extend=0, no_iterations=10, no_random_iterations=1,rank=T)
png(file=paste(input,"_saturation.png",sep=""),type="cairo")
MEDIPS.plotSaturation(sr.data)
dev.off()
data.set<-MEDIPS.getPositions(data=data.set, pattern="CG")
cr.data<-MEDIPS.coverageAnalysis(data=data.set, extend=0, no_iterations=10)
png(file=paste(input,"_coverage.png",sep=""), type="cairo")
MEDIPS.plotCoverage(cr.data)
dev.off()
MEDIPS.exportWIG(file=paste(input,"_rpm.wig",sep=""),data=data.set,raw=T,descr=paste(input,".rpm",sep=""))
er.data <- MEDIPS.CpGenrich(data=data.set)
write.table(er.data,paste(input,"_enrichment.txt",sep=""),sep="\t",quote=F)
dataPos.set<-genomeVectorStrand(data=data.set,extend=0,bin_size=50,strandSpec="+")
exportWigStrand(data=dataPos.set,file=paste(input,"_for_rpm.wig",sep=""),raw=T,descr=paste(input,"_for.rpm",sep=""),sample_name=paste(input,"_for",sep=""))
rm(dataPos.set)
gc()
dataNeg.set<-genomeVectorStrand(data=data.set,extend=0,bin_size=50,strandSpec="-")
exportWigStrand(data=dataNeg.set,file=paste(input,"_rev_rpm.wig",sep=""),raw=T,descr=paste(input,"_rev.rpm",sep=""),sample_name=paste(input,"_rev",sep=""),color="0.255.0")
