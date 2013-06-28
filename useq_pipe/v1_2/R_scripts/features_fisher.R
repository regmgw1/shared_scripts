#!/usr/local/bin/Rscript

args <- commandArgs();
path2input <-args[length(args)-1]
dmrSet<-args[length(args)]

##load data
data<-read.table(paste(path2input,dmrSet,"_enrichment.txt",sep=""),head=T,sep="\t")
LZ<-length(data[,1])
FirstB<-c()
SecondB<-c()
ThirdB<-c()
FourthB<-c()
Results<-c()
Totaldmr<-data[1,3]
Totalpeaks<-data[1,5]
Results.df<-data.frame()
for (z in 1:LZ){
KK<-c()
MM<-c()
FirstB<-data[z,2]
SecondB<-(data[z,4] - FirstB)
ThirdB<-Totaldmr-FirstB
FourthB<-(Totalpeaks-Totaldmr)-SecondB
KK<-c(FirstB,SecondB,ThirdB,FourthB)
MM<-matrix(KK,2,2)
fish<-fisher.test(MM, alternative="greater")
pval<-fish[1]
odds<-fish[3]
tmp.v<-c(as.character(data[z,1]),pval,odds)
Results.df<-rbind(Results.df,t(tmp.v))
print("______________________")
}
Results.df
P <- as.data.frame(lapply(Results.df,unlist))
P$p.value_FDRadjust = p.adjust(P$p.value, "fdr")
P$p.value_Bonferroniadjust = p.adjust(P$p.value, "bonferroni")
colnames(P)[1]<-"Feature"
write.table(P, paste(path2input,dmrSet,"_fisherOut.txt",sep=""), sep="\t", quote=F, row.names=F)
