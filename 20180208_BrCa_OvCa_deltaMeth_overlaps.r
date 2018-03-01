# This code compares the distribution of delta20 windows in BrCa
# vs OvCa

# Calculate the distribution of peaks
cancerType <- c("BrCa","OvCa")
annotations <- list.files()[grep("homer",list.files())]
annotations <- annotations[grep(paste0(cancerType,collapse="|"),annotations)]
for(i in 1:length(annotations)){
	toName <- gsub(".txt","",annotations[i])
	reading <- read.table(annotations[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
	if(length(grep("_annStats",toName))<1){
		colnames(reading)[1] <- "id"
		reading <- reading[order(reading[,"id"]),]
	}
	assign(toName,reading)
}

rm(annotations)
rm(i)
rm(reading)
rm(toName)

annotation <- OvCa_0_delta10_homer_annStats[,1]
annotation <- gsub("\\?","",annotation)
annotation <- unique(annotation)
annotation <- c(annotation,"LINE\\|L1","LINE\\|L2")

simple.annotation <- as.list(rep(NA,2))
names(simple.annotation) <- cancerType

for(i in 1:length(cancerType)){
	type <- cancerType[i]
	annotationFile <- paste0(type,"_0_delta20_homer")
	annotationFile <- get(annotationFile)
	ids <- annotationFile[,"id"]
	Detailed.Annotation <- annotationFile$Detailed.Annotation
	Simple.Annotation <- rep(NA,length(Detailed.Annotation))
	for(j in 1:length(annotation)){
		feature <- annotation[j]
		genomic.windows <- grep(feature,Detailed.Annotation,ignore.case=T)
		Simple.Annotation[genomic.windows] <- feature
		names(Simple.Annotation) <- ids
	}
	simple.annotation[[i]] <- Simple.Annotation
}

simple.annotation.toplot <- rbind(table(simple.annotation[[1]]),table(simple.annotation[[2]]))
rownames(simple.annotation.toplot) <- names(simple.annotation)
simple.annotation.toplot <- simple.annotation.toplot[,order(colSums(simple.annotation.toplot),decreasing=T)]
colnames(simple.annotation.toplot) <- gsub("\\\\\\|","\\|",colnames(simple.annotation.toplot))

png("BrCa_OvCa_0_delta20_genomic_distribution.png",w=6,h=6,unit="in",res=300)
barplot(simple.annotation.toplot,beside=T,las=2,ylab="Number of windows",cex.names=0.7,cex.axis=0.7,col=c("darkgrey","lightgrey"))
legend("topright",legend=rownames(simple.annotation.toplot),fill=c("darkgrey","lightgrey"))
dev.off()

# Assessing the overlap between BrCa and OvCa delta10 and delta20 windows

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(stringr)
require(venneuler)

rdsFiles <- list.files()[grep("20.RDS",list.files())]

for(i in 1:length(rdsFiles)){
rds <- readRDS(rdsFiles[i])
rds <- GRangesList(rds)
rds <- unlist(rds)
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

library(ChIPpeakAnno)
overlaps <- findOverlapsOfPeaks(BrCa_0_delta20.gr,OvCa_0_delta20.gr, maxgap=1000)
overlaps[[1]][1,3] <- 7729274-sum(overlaps[[1]][2:4,3])
overlaps.peaks <- overlaps$overlappingPeaks
table(overlaps.peaks[["BrCa_0_delta20.gr///OvCa_0_delta20.gr"]]$overlapFeature)

v <- overlaps[[1]][2:4,3]
conditions <- rev(unlist(lapply(str_split(colnames(overlaps[[1]]),"_"),"[[",1))[1:2])
names(v) <- c(conditions,paste(conditions,collapse="&"))
v <- venneuler(v)
png("BrCa_OvCa_0_delta20_venn.png",w=6,h=6,unit="in",res=300)
plot(v)
dev.off()

BrCa.ids <- paste(seqnames(BrCa_0_delta20.gr),start(BrCa_0_delta20.gr)-1,sep="_")
OvCa.ids <- paste(seqnames(OvCa_0_delta20.gr),start(OvCa_0_delta20.gr)-1,sep="_")

BrCa_0_delta20.gr <- BrCa_0_delta20.gr[which(BrCa.ids %in% BrCa_0_delta20_homer[,"id"])]
OvCa_0_delta20.gr <- OvCa_0_delta20.gr[which(OvCa.ids %in% OvCa_0_delta20_homer[,"id"])]
BrCa.ids <- BrCa.ids[which(BrCa.ids %in% BrCa_0_delta20_homer[,"id"])]
OvCa.ids <- OvCa.ids[which(OvCa.ids %in% OvCa_0_delta20_homer[,"id"])]

BrCa.annotation <- simple.annotation[[1]][match(BrCa.ids, names(simple.annotation[[1]]))]
OvCa.annotation <- simple.annotation[[2]][match(OvCa.ids, names(simple.annotation[[2]]))]

annotation <- sort(unique(BrCa.annotation))

overlaps.by.region <- as.list(rep(NA,length(annotation)))
names(overlaps.by.region) <- annotation

#fileName <- paste0("BrCa_OvCa_0_delta20_",region,"_venn.png")
png("BrCa_OvCa_0_delta20_byRegion.png",w=6,h=6,unit="in",res=300)
par(mfrow=c(6,4),mar=c(1,1,1,1))
for(i in 1:length(annotation)){
region <- annotation[i]
cat("\n",i," of ",length(annotation)," ",region,"\n")
brca.index <- which(BrCa.annotation==region)
brca <- BrCa_0_delta20.gr[brca.index]
ovca.index <- which(OvCa.annotation==region)
ovca <- OvCa_0_delta20.gr[ovca.index]
overlaps.tmp <- findOverlapsOfPeaks(brca,ovca, maxgap=1000)
overlaps.by.region[[i]] <- overlaps.tmp[[1]]
v <- overlaps.tmp[[1]][2:4,3]
names(v) <- c(conditions,paste(conditions,collapse="&"))
v <- venneuler(v)
plot(v,main=region)
}
dev.off()

hg19_repeat <- read.table(gzfile("Resources/hg19_rmsk.tsv.gz"),header=T,comment='$',stringsAsFactors=F)
l1s <- grep("L1",hg19_repeat$repFamily)
l2s <- grep("L2",hg19_repeat$repFamily)
l1.bp <- sum(hg19_repeat[l1s,"genoEnd"]-hg19_repeat[l1s,"genoStart"])
l2.bp <- sum(hg19_repeat[l2s,"genoEnd"]-hg19_repeat[l2s,"genoStart"])

q <- overlaps.by.region[[6]][4,3]*200
m <- (overlaps.by.region[[6]][3,3]*200)+q
n <- l1.bp-m
k <- overlaps.by.region[[6]][2,3]*200+q
phyper(q=q,m=m,n=n,k=k)



