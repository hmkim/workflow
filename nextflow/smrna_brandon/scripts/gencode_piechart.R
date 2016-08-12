
# USAGE: gencode_pichart.R <gencode gene counts> <sample_name> <output dir>

############################################################################
## Description:
## create a pie chart in PDF format from GENCODE classification table
##
## Author: Jared Evans
## Date: 5/22/14
##
## Parameters:
## <GENCODE counts> - tab delimited file from HTSeq output
## <sample name> - Name of Sample
## <output dir> - Directory where output PDF should be written
##
############################################################################



stdin <- commandArgs(TRUE)

if(length(stdin) != 3){
	stop("ERROR! Incorrect number of arguments. \nUSAGE: gencode_piechart.R <gencode_gene_counts> <sample_name> <output dir>")
}
gene_counts.path <- stdin[1]
sample <- stdin[2]
output.dir <- stdin[3]
keep<-c("lincRNA","miRNA","misc_RNA","protein_coding","rRNA","snoRNA","snRNA","Mt_rRNA","Mt_tRNA")
#keep<-c("lincRNA","miRNA","misc_RNA","protein_coding","rRNA","snoRNA","snRNA")

genecnt<-read.table(gene_counts.path,sep="\t")
cnt.bytype<-rowsum(genecnt[,2],group=genecnt$V4)
#cnt.bytype100<-cnt.bytype[rowSums(cnt.bytype>100)>0,]
cnt.bytypekeep<-cnt.bytype[rownames(cnt.bytype) %in% keep,]
# add all other counts into an other category
other <- sum(cnt.bytype[!rownames(cnt.bytype) %in% keep,])
newnames<-c(names(cnt.bytypekeep),"other")
cnt.bytypekeep<-c(cnt.bytypekeep,other)
names(cnt.bytypekeep) <- newnames

#pdf(paste(output.dir,"/",sample,"_gencode_piechart.pdf",sep=""),width=11,height=8.5)
png(paste(output.dir,"/",sample,"_gencode_piechart.png",sep=""),width=2400,height=2400,res=300)
#par(mar = c(10,6,5,3))
par(mar=c(0,0,4.1,0))
pie(cnt.bytypekeep,labels=names(cnt.bytypekeep),main=paste(sample," RNA Quantification",sep=""),cex=0.8)
dev.off()

write.table(cnt.bytype,paste(output.dir,"/",sample,"_gencode_counts.txt",sep=""),sep="\t",quote=F,col.names=F)

