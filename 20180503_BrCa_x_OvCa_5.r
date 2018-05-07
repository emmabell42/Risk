# This code generates odds ratios for the enrichment of delta20 windows 
# in different genomic context

setwd("/data/emmabell42/risk/BrCa_x_OvCa")

# Load packages

library(GenomicRanges)

# Read in files
rdsFiles <- list.files("../")[grep("20.RDS",list.files("../"))]
rdsFiles <- c(rdsFiles,list.files("..")[grep("allWindows.RDS",list.files(".."))])
for(i in 1:length(rdsFiles)){
toRead <- paste0("../",rdsFiles[i])
rds <- readRDS(toRead)
rds <- GRangesList(rds)
rds <- unlist(rds)
#rds <- reduce(rds) # added to join adjacent windows
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

common <- findOverlaps(BrCa_0_delta20.gr,OvCa_0_delta20.gr)
brca.common.index <- queryHits(common)
ovca.common.index <- subjectHits(common)

cancerTypes <- c("BrCa","OvCa")
patientTypes <- c("cases","controls")
regionTypes <- c("common","unique")

annotations <- list.files("..")[grep("20_homer",list.files(".."))]
for(i in 1:length(annotations)){
	toRead <- paste0("../",annotations[i])
	toName <- gsub(".txt","",annotations[i])
	reading <- read.table(toRead,head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
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


regionsOfInterest <- c("promoter","exon","intron","tts","intergenic","LINE","SINE","LTR","satellite","cpg-island","low_complexity","simple_repeat")

for(i in 1:length(regionTypes)){
	region <- regionTypes[i]
	if(region=="common"){
		simpleAnnotation <- BrCa_0_delta20_homer[brca.common.index,"Annotation"]
		detailedAnnotation <- BrCa_0_delta20_homer[brca.common.index,"Detailed.Annotation"]
	} else {
		simpleAnnotation <- c(BrCa_0_delta20_homer[which(!rownames(BrCa_0_delta20_homer) %in% brca.common.index),"Annotation"],OvCa_0_delta20_homer[which(!rownames(OvCa_0_delta20_homer) %in% ovca.common.index),"Annotation"])
		detailedAnnotation <- c(BrCa_0_delta20_homer[which(!rownames(BrCa_0_delta20_homer) %in% brca.common.index),"Detailed.Annotation"],OvCa_0_delta20_homer[which(!rownames(OvCa_0_delta20_homer) %in% ovca.common.index),"Detailed.Annotation"])
	}
	simpleAnnotation <- lapply(strsplit(simpleAnnotation," "), `[`, 1)
	simpleAnnotation <- unlist(simpleAnnotation)
	simpleAnnotation[grep("intron",simpleAnnotation)] <- "Intron"
	simpleAnnotation[grep("non-coding",simpleAnnotation)] <- "Intergenic"
	simpleAnnotation[grep("promoter-TSS",simpleAnnotation)] <- "Promoter"
	simpleAnnotation[grep("exon",simpleAnnotation)] <- "Exon"

	for(j in 1:length(regionsOfInterest)){
		annotation <- regionsOfInterest[j]
		detailedAnnotation[grep(annotation,detailedAnnotation,ignore.case=T)] <- annotation
	}
	detailedAnnotation[grep("CpG",detailedAnnotation)] <- "CpG-Island"
	detailedAnnotation[grep("\\|",detailedAnnotation)] <- simpleAnnotation[grep("\\|",detailedAnnotation)]
	
	simpleName <- paste0(region,".simple")
	detailedName <- paste0(region,".detailed")
	
	assign(simpleName,simpleAnnotation)
	assign(detailedName,detailedAnnotation)

	annotationTable <- cbind(Local=table(detailedAnnotation),Global=NA)

	for(j in 1:nrow(annotationTable)){
		toFind <- rownames(annotationTable)[j]
		annotationTable[j,2] <- as.numeric(annStats[grep(toFind,annStats[,1],ignore.case=T)[1],3])
	}
		
	toName <- paste0(region,".AnnotationTable")
	assign(toName,annotationTable)
}

annotationMatrix <- array(NA,dim=c(nrow(annotationTable),4))
colnames(annotationMatrix) <- c("OR","LCI","UCI","P")
rownames(annotationMatrix) <- rownames(annotationTable)

alpha <- 0.05

for(j in 1:nrow(annotationMatrix)){

	signifBPwAnnot <- as.numeric(common.AnnotationTable[j,1])
	signifBPwoAnnot <- sum(as.numeric(common.AnnotationTable[,1]))-signifBPwAnnot
	allBPwAnnot <- as.numeric(common.AnnotationTable[j,2])/200
	allBPwoAnnot <- sum(as.numeric(common.AnnotationTable[,2])/200)-allBPwAnnot
	
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

annotationToPlot <- c(which(annotationMatrix[,2]<1 & annotationMatrix[,3]<1),which(annotationMatrix[,2]>1 & annotationMatrix[,3]>1))
	
annotationMatrix <- annotationMatrix[annotationToPlot,]
annotationMatrix <- annotationMatrix[order(annotationMatrix[,"OR"],decreasing=T),]

annotationMatrix <- annotationMatrix[2:nrow(annotationMatrix),]

plotName <- paste0("Barplot_common_annotation_OR.png")
png(plotName,h=6,w=12,res=300,unit="in")
par(mar=c(7,5,2,2))
x <- barplot(annotationMatrix[,"OR"])
ylim <- c(0,max(annotationMatrix[,"UCI"])*1.1)
barplot(annotationMatrix[,"OR"],ylim=ylim,names=rownames(annotationMatrix),las=2,main="Common",ylab="Odds Ratio")
arrows(x, annotationMatrix[,"UCI"], x, annotationMatrix[,"LCI"], length=0.05, angle=90, code=3)
abline(h=1,lty=2,lwd=2,col="red")

dev.off()


}

