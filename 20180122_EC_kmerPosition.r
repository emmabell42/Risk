library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
meth.seq <- getSeq(hg19,names=GRanges(seqnames=as.character(meth$chr),ranges=IRanges(start=meth$start,end=meth$end)))

this.kmer <- kmers[1]
posList <- sapply(seqs,function(x)start(matchPattern(pattern=this.kmer,subject=x)))

resolution <- 1000
posMatrix <- array(0,dim=c(length(seqs),resolution))
for(i in  which(sapply(posList,length)>0)){
	thisx <- round(resolution*(posList[[i]]/width(seqs[i])))
	posMatrix[i,thisx] <- 1
}
png(file="edtest.png")
plot(colMeans(posMatrix),xlab="relative position",ylab="enrichment (a.u.)",type="l",lwd=2)
dev.off()