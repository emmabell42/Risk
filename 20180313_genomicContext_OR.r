# This code generates odds ratios for the enrichment of delta20 windows 
# in different genomic context


annotations <- list.files()[grep("20_homer",list.files())]
for(i in 1:length(annotations)){
	toName <- gsub(".txt","",annotations[i])
	reading <- read.table(annotations[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
	if(length(grep("_annStats",toName))==1){
		reading[13,] <- NA
		reading <- na.omit(reading)
	}
	
	else {
	
	colnames(reading)[1] <- "id"
	reading <- reading[order(reading[,"id"]),]
	
	}
	
	assign(toName,reading)
}

rm(annotations)
rm(i)
rm(reading)
rm(toName)

cancerType <- c("BrCa","OvCa")

alpha <- 0.05

for(i in 1:length(cancerType)){

	cancer <- cancerType[i]
	annStats <- paste0(cancer,"_0_delta20_homer_annStats")
	annStats <- get(annStats)
	annStats <- annStats[!duplicated(annStats[,1]),]
	annStats[grep("\\?",annStats[,1]),1] <- NA
	annStats <- na.omit(annStats)
	
	annotationMatrix <- array(NA,dim=c(nrow(annStats),4))
	colnames(annotationMatrix) <- c("OR","LCI","UCI","P")
	rownames(annotationMatrix) <- annStats[,1]
	
	
	for(j in 1:nrow(annotationMatrix)){
		signifBPwAnnot <- as.numeric(annStats[j,2])
		signifBPwoAnnot <- sum(as.numeric(annStats[,2]))-signifBPwAnnot
		allBPwAnnot <- as.numeric(annStats[j,3])/200
		allBPwoAnnot <- sum(as.numeric(annStats[,3])/200)-allBPwAnnot
		
		OR <- (signifBPwAnnot/signifBPwoAnnot)/(allBPwAnnot/allBPwoAnnot)
		
		siglog <- sqrt((1/signifBPwAnnot) + (1/signifBPwoAnnot) + (1/allBPwAnnot) + (1/allBPwoAnnot))
		zalph <- qnorm(1 - alpha/2)
		logOR <- log(OR)
		loglo <- logOR - zalph * siglog
		loghi <- logOR + zalph * siglog
		
		ORlo <- exp(loglo)
		ORhi <- exp(loghi)

		annotationMatrix[j,1:3] <- c(OR,ORlo,ORhi)
		
		chiTable <- matrix(c(signifBPwAnnot,signifBPwoAnnot, allBPwAnnot, allBPwoAnnot),nrow=2,ncol=2,byrow=T)
		p <- chisq.test(chiTable)[[3]]
		annotationMatrix[j,4] <- p
	}
	
	annotationMatrix[,4] <- p.adjust(annotationMatrix[,4],method="BH")
	
	toName <- paste0(cancer,".ORtable")
	assign(toName,annotationMatrix)
	
	annotationToPlot <- c(which(annotationMatrix[,2]<1 & annotationMatrix[,3]<1),which(annotationMatrix[,2]>1 & annotationMatrix[,3]>1))
	annotationToPlot <- sort(annotationToPlot)
	
	annotationMatrix <- annotationMatrix[annotationToPlot,]
	annotationMatrix <- annotationMatrix[order(annotationMatrix[,"OR"],decreasing=T),]
	
	plotName <- paste0("Barplot_",cancer,"_annotation_OR.png")
	png(plotName,h=6,w=12,res=300,unit="in")
	par(mar=c(7,5,2,2))
	x <- barplot(annotationMatrix[,"OR"])
	ylim <- c(0,max(annotationMatrix[,"UCI"])*1.1)
	barplot(annotationMatrix[,"OR"],ylim=ylim,names=rownames(annotationMatrix),las=2,main=cancer,ylab="Odds Ratio")
	arrows(x, annotationMatrix[,"UCI"], x, annotationMatrix[,"LCI"], length=0.05, angle=90, code=3)
	abline(h=1,lty=2,lwd=2,col="red")
	
	dev.off()
	
}


library(GenomicRanges)
library(stringr)

rdsFiles <- list.files()[grep("allWindows.RDS",list.files())]
rdsFiles <- c(rdsFiles,list.files()[grep("_delta20.RDS",list.files())])

for(i in 1:length(rdsFiles)){
rds <- readRDS(rdsFiles[i])
rds <- GRangesList(rds)
rds <- unlist(rds)
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

# Get the number of bp of each annotation in the whole genome

simpleAnnotation <- OvCa_0_delta20_homer$Annotation
simpleAnnotation <- lapply(strsplit(simpleAnnotation," "), `[`, 1:2)
simpleAnnotation <- simpleAnnotation[!duplicated(simpleAnnotation)]
simpleAnnotation[grep("intron",simpleAnnotation)] <- "intron"
simpleAnnotation[grep("Intergenic",simpleAnnotation)] <- "intergenic"
simpleAnnotation[grep("TTS",simpleAnnotation)] <- "TTS"
simpleAnnotation[grep("3'",simpleAnnotation)] <- "3UTR"
simpleAnnotation[grep("non-coding",simpleAnnotation)] <- "non-coding"
simpleAnnotation[grep("promoter-TSS",simpleAnnotation)] <- "promoter-TSS"
simpleAnnotation[grep("exon",simpleAnnotation)] <- "exon"
simpleAnnotation[grep("5'",simpleAnnotation)] <- "5UTR"
simpleAnnotation <- simpleAnnotation[!duplicated(simpleAnnotation)]



directions <- c("hypo","hyper")

for(i in 1:length(cancerType)){

	cancer <- cancerType[i]
	delta20 <- paste0(cancer,"_0_delta20.gr")
	delta20.gr <- get(delta20)
	
	ids <- paste(seqnames(delta20.gr),start(delta20.gr)-1,sep="_")
	
	annStats <- paste0(cancer,"_0_delta20_homer_annStats")
	annStats <- get(annStats)
	annStats <- annStats[!duplicated(annStats[,1]),]
	annStats[grep("\\?",annStats[,1]),1] <- NA
	annStats <- na.omit(annStats)
	
	homer <- paste0(cancer,"_0_delta20_homer")
	homer <- get(homer)
	
	delta20.gr <- delta20.gr[which(ids %in% homer[,"id"])]
	ids <- ids[which(ids %in% homer[,"id"])]
	
	detailedAnnotation <- homer$Detailed.Annotation
	
	for(j in 1:nrow(annStats)){
		annotation <- annStats[j,1]
		if(length(grep("\\?",annotation))<1){
		detailedAnnotation[grep(annotation,detailedAnnotation,ignore.case=T)] <- annotation
		}
	}

	detailedAnnotation[grep("CpG",detailedAnnotation)] <- "CpG-Island"
	
	for(k in 1:length(directions)){
		
	direction <- directions[k]
	means <- mcols(delta20.gr)[,1:2]
	
	if(direction=="hypo"){
	
	index <- which(means[,2]<means[,1])
	direction.gr <- delta20.gr[index]
	
	}
	
	else{
	
	index <- which(means[,2]>means[,1])
	direction.gr <- delta20.gr[index]
	
	}

	directionAnnotation <- detailedAnnotation[index]
	directionAnnotationTable <- as.matrix(cbind(table(directionAnnotation),NA))
	
	for(l in 1:nrow(directionAnnotationTable)){
		
		toFind <- rownames(directionAnnotationTable)[l]
		directionAnnotationTable[l,2] <- annStats[which(annStats[,1]==toFind),3]
		
		}
		
	toName <- paste0(cancer,".",direction,".AnnotationTable")
	assign(toName,directionAnnotationTable)

	annotationMatrix <- array(NA,dim=c(nrow(directionAnnotationTable),4))
	colnames(annotationMatrix) <- c("OR","LCI","UCI","P")
	rownames(annotationMatrix) <- rownames(directionAnnotationTable)
	
	
	for(j in 1:nrow(annotationMatrix)){
		signifBPwAnnot <- as.numeric(directionAnnotationTable[j,1])
		signifBPwoAnnot <- sum(as.numeric(directionAnnotationTable[,1]))-signifBPwAnnot
		allBPwAnnot <- as.numeric(directionAnnotationTable[j,2])/200
		allBPwoAnnot <- sum(as.numeric(directionAnnotationTable[,2])/200)-allBPwAnnot
		
		OR <- (signifBPwAnnot/signifBPwoAnnot)/(allBPwAnnot/allBPwoAnnot)
		
		siglog <- sqrt((1/signifBPwAnnot) + (1/signifBPwoAnnot) + (1/allBPwAnnot) + (1/allBPwoAnnot))
		zalph <- qnorm(1 - alpha/2)
		logOR <- log(OR)
		loglo <- logOR - zalph * siglog
		loghi <- logOR + zalph * siglog
		
		ORlo <- exp(loglo)
		ORhi <- exp(loghi)

		annotationMatrix[j,1:3] <- c(OR,ORlo,ORhi)
		
		chiTable <- matrix(c(signifBPwAnnot,signifBPwoAnnot, allBPwAnnot, allBPwoAnnot),nrow=2,ncol=2,byrow=T)
		p <- chisq.test(chiTable)[[3]]
		annotationMatrix[j,4] <- p
	}
	
	annotationMatrix[,4] <- p.adjust(annotationMatrix[,4],method="BH")
	
	toName <- paste0(cancer,".",direction,".ORtable")
	assign(toName,annotationMatrix)
	
	annotationToPlot <- c(which(annotationMatrix[,2]<1 & annotationMatrix[,3]<1),which(annotationMatrix[,2]>1 & annotationMatrix[,3]>1))
	annotationToPlot <- sort(annotationToPlot)
	
	annotationMatrix <- annotationMatrix[annotationToPlot,]
	annotationMatrix <- annotationMatrix[order(annotationMatrix[,"OR"],decreasing=T),]
	
	plotName <- paste0("Barplot_",cancer,"_",direction,"_annotation_OR.png")
	png(plotName,h=6,w=12,res=300,unit="in")
	par(mar=c(7,5,2,2))
	x <- barplot(annotationMatrix[,"OR"])
	ylim <- c(0,max(annotationMatrix[,"UCI"])*1.1)
	barplot(annotationMatrix[,"OR"],ylim=ylim,names=rownames(annotationMatrix),las=2,main=paste(cancer,direction),ylab="Odds Ratio")
	arrows(x, annotationMatrix[,"UCI"], x, annotationMatrix[,"LCI"], length=0.05, angle=90, code=3)
	abline(h=1,lty=2,lwd=2,col="red")
	
	dev.off()

	
	}
}
