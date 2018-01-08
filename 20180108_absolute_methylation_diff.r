# This code plots the absolute differences in methylation values for cases and controls

bin_0_2 <- c()
bin_2_5 <- c()
bin_5_10 <- c()
bin_10_20 <- c()
bin_20_n <- c()

for(i in 1:length(gr.all)){
chr <- gr.all[[i]]
diffs <- chr$Casemean-chr$Controlmean
abs.diffs <- abs(diffs)
bin_0_2 <- length(which(abs.diffs<2))
bin_2_5 <- c(bin_2_5,length(which(abs.diffs>=2 & abs.diffs<5)))
bin_5_10 <- c(bin_5_10,length(which(abs.diffs>=5 & abs.diffs<10)))
bin_10_20 <- c(bin_10_20,length(which(abs.diffs>=10 & abs.diffs<20)))
bin_20_n <- c(bin_20_n,length(which(abs.diffs>20)))
}

bins <- rbind(bin_0_2,bin_2_5,bin_5_10,bin_10_20,bin_20_n)
colnames(bins) <- c(1:22,"X")
write.table(bins,"abs_meth_diff_chr.txt",sep="\t",quote=F)

n_bins <- colSums(bins)
bins.prop <- bins

for(i in 1:23){
bins.prop[,i] <- bins.prop[,i]/n_bins[i]
}

png("abs_meth_diff_chr.png",h=3,w=6,unit="in",res=300)
par(mar=c(3,5,2,2))
barplot(bins.prop,xlab="Chr",ylab="Proportion of windows",cex.axis=0.6,cex.names=0.5)
dev.off()

bins.total <- rowSums(bins)
bins.total.prop <- bins.total/sum(bins.total)

png("abs_meth_diff_overall.png",h=6,w=3,unit="in",res=300)
par(mar=c(3,5,2,2))
barplot(cbind(bins.total.prop),xlab="Chr",ylab="Proportion of windows",cex.axis=0.6,cex.names=0.5)
dev.off()

# Do the CIs overlap?

chrX.cases <- IRanges(chr$CaseLCI,chr$CaseUCI)
chrX.controls <- IRanges(chr$ControlLCI,chr$ControlUCI)

overlaps <- c()
for(i in 1:length(chrX.cases)){
case <- chrX.cases[[i]]
control <- chrX.controls[[i]]
overlaps <- c(overlaps,length(which(case %in% control)))
}
png("CI_overlaps.png")
hist((overlaps))
dev.off()

# Write out genomic regions with >20% delta meth

bins <- as.list(rep(NA,23))
for(i in 1:length(gr.all)){
chr <- gr.all[[i]]
diffs <- chr$Casemean-chr$Controlmean
abs.diffs <- abs(diffs)
bins[[i]] <- which(abs.diffs>20)
}

gr.diff <- gr.all
lengths <- c()
for(i in 1:23){
gr.diff[[i]] <- gr.diff[[i]][bins[[i]]]
lengths <- c(lengths,length(bins[[i]]))
}

lengths.200 <- lengths*200

hg19 <- read.table("/data/emmabell42/resources/hg19.chrom.sizes.txt",stringsAsFactors=F)
hg19 <- hg19[c(1:7,9:19,22,20,24,23,8),]

chr.props <- lengths.200/hg19[,3]
chr.pc <- chr.props*100 

png("chr_pc_delta20.png",h=3,w=6,unit="in",res=300)
par(mar=c(3,5,1,1))
barplot(chr.pc,ylab="Percentage of chr",names=c(1:22,"X"),cex.names=0.5,cex.axis=0.6)
dev.off()

df.diff <- c() 
for(i in 1:23){
df.diff <- rbind(df.diff,data.frame(seqnames=seqnames(gr.diff[[i]]),starts=start(gr.diff[[i]])-1,ends=end(gr.diff[[i]])))
}
id <- paste(df.diff[,1],df.diff[,2],sep="_")
df.diff <- cbind(df.diff,id,"","0")
write.table(df.diff,"delta_meth_20_CI_0.txt",sep="\t",quote=F,row.names=F,col.names=F)

PATH=$PATH:/data/seqtools/homer/bin/
PATH=$PATH:/data/seqtools/weblogo/
PATH=$PATH:/data/seqtools/samtools-1.1/


annotatePeaks.pl delta_meth_20_CI_0.txt hg19 > homer_delta_meth_20_CI_0.txt

nice R
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

png("delta_meth_20_CI_0_annot.png",h=6,w=3,unit="in",res=300)
par(mar=c(7,5,2,2))
barplot(table(reading$Annotation),las=2,ylab="Number of peaks",cex.axis=0.6)
dev.off()
