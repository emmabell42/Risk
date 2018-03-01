# This code performs Fisher's tests on the distribution of 
# delta20 windows in repetitive elements

library(GenomicRanges)
library(stringr)
library(vioplot)
library(epitools)

hg19_repeat <- read.table(gzfile("Resources/hg19_rmsk.tsv.gz"),header=T,comment='$',stringsAsFactors=F)

rdsFiles <- list.files()[grep("allWindows.RDS",list.files())]
rdsFiles <- c(rdsFiles,list.files()[grep("_20.RDS",list.files())])

for(i in 1:length(rdsFiles)){
rds <- readRDS(rdsFiles[i])
rds <- GRangesList(rds)
rds <- unlist(rds)
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

cancerType <- c("BrCa","OvCa")

goodChrs <- paste0("chr",c(1:22,"X"))
goodRanges <- which(hg19_repeat$genoName %in% goodChrs)

chr <- hg19_repeat$genoName[goodRanges]
starts <- as.numeric(hg19_repeat$genoStart[goodRanges])
ends <- as.numeric(hg19_repeat$genoEnd[goodRanges])
repFamily <- hg19_repeat$repFamily[goodRanges]
repeats.gr <- GRanges(seqnames=chr,ranges=IRanges(start=starts,end=ends))
mcols(repeats.gr) <- data.frame(repFamily)
repTypes <- names(table(repFamily))

progress <- seq(5,100,5)

for(i in 1:length(cancerType)){
	
	#cat("\n","Began working on",cancerType[i],"at",as.character(Sys.time()),
	#"\n","Percentage complete",
	#"\n","0","               ","100","\n")

	cancer <- cancerType[i]
	#allWindows <- paste0(cancer,"_allWindows.gr")
	#delta20 <- paste0(cancer,"_0_delta20.gr")
	#allWindows.gr <- get(allWindows)
	#delta20.gr <- get(delta20)
	#exp1 <- c()
	#obs1 <- c()
	
	#previous.progress <- 0
	#cat("Calculating the observed and expected distributions","\n")
	
	#for(j in 1:length(repTypes)){

	#	repType <- repTypes[j]
	#	index <- which(repFamily==repType)
	#	repType.gr <- repeats.gr[index]
	#
	#	x <- findOverlaps(repType.gr, allWindows.gr)
	#	y <- findOverlaps(repType.gr, delta20.gr)
	#
	#	exp1 <- c(exp1,length(x))
	#	obs1<-c(obs1, length(y))
	#	
	#	pc <- (j/length(repTypes))*100
	#	current.progress <- which.min(abs(pc-progress))
	#		
	#	if(current.progress>previous.progress){
	#		
	#		cat("█")
	#		
	#	previous.progress <- current.progress
	#	
	#	}
	#
	#}
	
	#exp2 <- c(exp1, length(allWindows.gr)-sum(exp1))
	#obs2 <- c(obs1, length(delta20.gr)-sum(obs1))

	#oe <- (obs2/length(delta20.gr)*100)/ (exp2/length(allWindows.gr)*100)

	#chi.table <- array(NA,c(length(repTypes),3))
	#rownames(chi.table) <- repTypes
	#colnames(chi.table) <- c("est","lower","upper")
	
	#previous.progress <- 0
	#cat("\n","Performing Chi Square tests","\n")

	#for(j in 1:length(repTypes)){
	#
	#	chi <- matrix(c(obs2[j],exp2[j],length(delta20.gr),length(allWindows.gr)),ncol=2,nrow=2)
	#	chi.table[j,1] <- fisher.test(chi)[[3]]
	#	chi.table[j,2:3] <- fisher.test(chi)[[2]]
		
	#	pc <- (j/length(repTypes))*100
	#	current.progress <- which.min(abs(pc-progress))
			
	#	if(current.progress>previous.progress){
			
	#		cat("█")
			
	#	previous.progress <- current.progress
		
	#	}

	#}
	
	chiName <- paste0(cancer,".chiTable")
	#assign(chiName,chi.table)
	
	chi.table <- get(chiName)
	
	chi.table2 <- chi.table
	chi.toPlot <- c(which(chi.table[,2]<1 & chi.table[,3]<1),which(chi.table[,2]>1 & chi.table[,3]>1))
	chi.toPlot <- sort(chi.toPlot)
	
	chi.table2 <- chi.table[chi.toPlot,]
	chi.table2 <- chi.table2[order(chi.table2[,"est"],decreasing=T),]
	
	plotName <- paste0("Barplot_",cancer,"_repeats_OR.png")
	png(plotName,h=6,w=12,res=300,unit="in")
	par(mar=c(7,5,2,2))
	x <- barplot(chi.table2[,"est"])
	ylim <- c(0,max(chi.table2[,"upper"])*1.1)
	barplot(chi.table2[,"est"],ylim=ylim,names=rownames(chi.table2),las=2,main=cancer,ylab="Odds Ratio")
	arrows(x, chi.table2[,"lower"], x, chi.table2[,"upper"], length=0.05, angle=90, code=3)
	abline(h=1,lty=2,lwd=2,col="red")
	
	dev.off()

}
