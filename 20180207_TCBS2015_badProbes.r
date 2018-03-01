# This code analyses the TCBS targets from the 2015 experiment

report <- read.table("Resources/JF_v2_1x_moderate stringent_1_Report.txt",skip=45, head=T,sep="\t",stringsAsFactors=F)

# Load the required R packages
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE) 
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
library(Biostrings)
library(gplots)

toSplit <- c(":","-")
regions <- strsplit(report$Interval,paste0(toSplit,collapse="|"))
regions <- unlist(regions)
regions <- matrix(regions,ncol=3,byrow=T)
regions <- as.data.frame(regions,stringsAsFactors=F)
for(i in 2:3){
regions[,i] <- as.numeric(regions[,i])
}
regions[grep("chrx",regions[,1]),1] <- "chrX"

# Create a GenomicRanges object of your regions 
regions.gr <- GRanges(seqnames=regions[,1],IRanges(start=regions[,2],end=regions[,3]))

# Retrieve the DNA sequences of your regions 
seqs <- getSeq(Hsapiens, seqnames(regions.gr), start(regions.gr), end(regions.gr), as.character = T)
names(seqs) <- report$TargetID
seqs <- DNAStringSet(x=seqs,use.names=T)

# Count oligonucleotide frequency of the sequences
oligos <- as.list(rep(NA,5))
cors <- oligos
dendros <- oligos

for(i in 1:length(oligos)){
	width <- i
	for(j in 1:length(seqs)){
	oligos[[i]] <- rbind(oligos[[i]],oligonucleotideFrequency(seqs[j],width=width))	
	}

oligos[[i]] <- na.omit(oligos[[i]])
cors[[i]] <- cor(t(oligos[[i]]),method="spearman")
dendros[[i]] <- hclust(as.dist(1-cors[[i]]))

}

for(i in 1:4){

fileName <- paste0("Heatmap_TCBStargets_frequency_nucleotide_",i,".png")
png(fileName,h=6,w=12,unit="in",res=300)
heatmap.2(oligos[[i]],trace="none",scale="row",labRow="",col=bluered(100),RowSideCol=as.character(report$HighCoverage),Rowv=as.dendrogram(dendros[[i]]),cexCol=4/i)
dev.off()

fileName <- paste0("Heatmap_TCBStargets_proportion_nucleotide_",i,".png")
png(fileName,h=6,w=12,unit="in",res=300)
heatmap.2(oligos.prop[[i]],trace="none",scale="row",labRow="",col=bluered(100),RowSideCol=as.character(report$HighCoverage),Rowv=as.dendrogram(dendros.prop[[i]]),cexCol=4/i)
dev.off()

fileName <- paste0("Heatmap_TCBStargets_oe_nucleotide_",i,".png")
png(fileName,h=6,w=12,unit="in",res=300)
heatmap.2(oligos.oe[[i]],trace="none",scale="row",labRow="",col=bluered(100),RowSideCol=as.character(report$HighCoverage),Rowv=as.dendrogram(dendros.oe[[i]]),cexCol=4/i)
dev.off()

}

bases <- colnames(oligos[[1]])

for(i in 1:4){

base <- bases[i]
toGrep <-  paste0("^",base)
toPlot <- oligos[[4]][,grep(toGrep,colnames(oligos[[4]]))]

fileName <- paste0("Heatmap_TCBStargets_frequency_kmer4_",i,".png")
png(fileName,h=6,w=12,unit="in",res=300)

heatmap.2(toPlot,trace="none",scale="row",labRow="",col=bluered(100),RowSideCol=as.character(report$HighCoverage),cexCol=1)
dev.off()

}

cors <- array(NA,dim=c(ncol(oligos[[5]]),2))
colnames(cors) <- c("Cor","P")
rownames(cors) <- colnames(oligos[[5]])

for(i in 1:ncol(oligos[[5]])){

thisKmer <- oligos[[5]][,i]
cors[i,1] <- cor.test(thisKmer,report$Coverage)[[4]]
cors[i,2] <- cor.test(thisKmer,report$Coverage)[[3]]

}

cors[,2] <- p.adjust(cors[,2],method="bonferroni")
write.table(cors,"kmer_5_cors.txt",sep="\t",quote=F)

cors <- cors[order(cors[,1]),]
badKmers <- rownames(cors)[1:148]
toPlot <- oligos[[5]][,which(colnames(oligos[[5]]) %in% badKmers)]
rowSideCols <- as.character(report$HighCoverage)[which(rowSums(toPlot)>0)]
toPlot <- toPlot[which(rowSums(toPlot)>0),]
png("badKmers.png",h=6,w=12,unit="in",res=300)
heatmap.2(toPlot,trace="none",scale="row",labRow="",col=bluered(100),RowSideCol=rowSideCols,cexCol=0.5)
dev.off()

badBases <- c("AAA","TTT")
badKmers <- rownames(cors)[grep(paste0(badBases,collapse="|"),rownames(cors))]
toPlot <- oligos[[5]][,which(colnames(oligos[[5]]) %in% badKmers)]
rowSideCols <- as.character(report$HighCoverage)[which(rowSums(toPlot)>0)]
toPlot <- toPlot[which(rowSums(toPlot)>0),]
png("badKmers_runs.png",h=6,w=12,unit="in",res=300)
heatmap.2(toPlot,trace="none",scale="row",labRow="",col=bluered(100),RowSideCol=rowSideCols,cexCol=0.5)
dev.off()


# k-mer content

badKmers <- c("AAAAA","TTTTT")
regions <- c("high","low")

high <- seqs[which(report$HighCoverage==1)]
low <- seqs[which(report$HighCoverage==0)]

enrichment <- array(NA,c(1000,length(regions)*length(badKmers)))
colnames(enrichment) <- as.vector(sapply(badKmers,function(X) paste(X,regions,sep="_")))

for(i in 1:length(badKmers)){

	this.kmer <- badKmers[i]

	for(j in 1:length(regions)){
		
		regionsToCalc <- get(regions[j])
		
		posList <- sapply(regionsToCalc,function(x)start(matchPattern(pattern=this.kmer,subject=x)))
		resolution <- 1000
		posMatrix <- array(0,dim=c(length(regionsToCalc),resolution))
		for(k in  which(sapply(posList,length)>0)){
			thisx <- round(resolution*(posList[[k]]/width(regionsToCalc[k])))
			posMatrix[k,thisx] <- 1
			
			whichColumn <- paste(this.kmer,regions[j],sep="_")
			enrichment[,whichColumn] <- colMeans(posMatrix)	
		}
	}
}

ylim <- c(0,max(enrichment)*1.1)

for(i in 1:length(badKmers)){
kmer <- badKmers[i]
columns <- grep(kmer,colnames(enrichment))
fileName <- paste0("kmer_enrichment_",kmer,".png")
png(fileName,h=6,w=6,unit="in",res=300)
plot(enrichment[,columns[1]],xlab="Relative Position",ylab="Enrichment (A.U.)",type="l",lwd=2,ylim=ylim,main=kmer)
points(enrichment[,columns[2]],col="red",type="l",lwd=2)
legend("topright",regions,fill=c("black","red"))
dev.off()
}

baduns <- oligos[[5]][,"TTTTT"]+oligos[[5]][,"AAAAA"]
P <- paste("P =",formatC(t.test(baduns~report$HighCoverage)[[3]], format = "e", digits = 2))
png("boxplot_badKmers.png",,h=6,w=3,unit="in",res=300)
boxplot(baduns~report$HighCoverage,col="lightgrey",pch=20,notch=T,names=c("Low","High"),xlab="Probe coverage",ylab="Frequency of AAAAA or TTTTT")
legend("topright",legend=P,bty="n")
dev.off()

# GC content

GC <- c()

for(i in 1:length(seqs)){
nC <- oligos[[1]][i,2]
nG <- oligos[[1]][i,3]
bp <- length(seqs[[i]])
GC <- c(GC,(nC+nG)/bp)
}

cor.test(GC,report$Coverage)
P <- paste("P =",formatC(t.test(GC~report$HighCoverage)[[3]], format = "e", digits = 2))
png("boxplot_GC.png",,h=6,w=3,unit="in",res=300)
boxplot(GC~report$HighCoverage,col="lightgrey",pch=20,notch=T,names=c("Low","High"),xlab="Probe coverage",ylab="GC%")
legend("topright",legend=P,bty="n")
dev.off()


