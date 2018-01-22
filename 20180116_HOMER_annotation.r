# This code uses HOMER to annotate the delta 20 methylation regions and then analyses the results in R

cd risk

PATH=$PATH:/data/seqtools/homer/bin/
PATH=$PATH:/data/seqtools/weblogo/
PATH=$PATH:/data/seqtools/samtools-1.1/

annotatePeaks.pl OvCa_0_delta20_CI.txt hg19 > homer_OvCa_0_delta20_CI.txt
annotatePeaks.pl BrCa_0_delta20_CI.txt hg19 > homer_BrCa_0_delta20_CI.txt

nice R

# Read in HOMER annotation output files

annotations <- list.files()[grep("homer",list.files())]
for(i in 1:length(annotations)){
	reading <- read.table(annotations[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
	reading$Annotation[grep("intron",reading$Annotation)] <- "intron"
	reading$Annotation[grep("promoter-TSS",reading$Annotation)] <- "promoter-TSS"
	reading$Annotation[grep("exon",reading$Annotation)] <- "exon"
	reading$Annotation[grep("non-coding",reading$Annotation)] <- "non-coding"
	reading$Annotation[grep("TTS",reading$Annotation)] <- "TTS"
	cat(annotations[i]," ",table(reading$Annotation),"\n")
	assign(annotations[i],reading)
}

cancerType <- c("OvCa","BrCa")
min.diff <- 0

# Plot the genomic distribution of delta 20 meth windows

for(i in 1:length(cancerType)){
	
	plotting <- get(annotations[i])
	toName <- paste0(cancerType[i],"_",min.diff,"_barplot_annotation_numbers.png")

	png(toName,h=6,w=3,unit="in",res=300)
	par(mar=c(7,5,2,2))
	barplot(table(plotting$Annotation),las=2,ylab="Number of windows",cex.axis=0.6)
	dev.off()
	
	toName <- paste0(cancerType[i],"_",min.diff,"_barplot_annotation_pc.png")

	png(toName,h=6,w=3,unit="in",res=300)
	par(mar=c(7,5,2,2))
	barplot((table(plotting$Annotation)/sum(table(plotting$Annotation)))*100,las=2,ylab="% of windows",cex.axis=0.6)
	dev.off()

}