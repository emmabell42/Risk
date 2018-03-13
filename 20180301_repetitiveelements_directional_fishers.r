# This code performs Fisher's tests on the distribution of 
# delta20 windows in repetitive elements

library(GenomicRanges)
library(stringr)
library(vioplot)
library(epitools)

hg19_repeat <- read.table(gzfile("Resources/hg19_rmsk.tsv.gz"),header=T,comment='$',stringsAsFactors=F)

rdsFiles <- list.files()[grep("allWindows.RDS",list.files())]
rdsFiles <- c(rdsFiles,list.files()[grep("_delta20.RDS",list.files())])

for(i in 1:length(rdsFiles)){
rds <- readRDS(rdsFiles[i])
rds <- GRangesList(rds)
rds <- unlist(rds)
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

cancerType <- c("BrCa","OvCa")
directions <- c("hypo","hyper")

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
	
	cat("\n","Began working on",cancerType[i],"at",as.character(Sys.time()),
	"\n","Percentage complete",
	"\n","0","               ","100","\n")

	cancer <- cancerType[i]
	allWindows <- paste0(cancer,"_allWindows.gr")
	delta20 <- paste0(cancer,"_0_delta20.gr")
	allWindows.gr <- get(allWindows)
	delta20.gr <- get(delta20)
	exp1 <- c()
	obs1 <- c()
	
	for(j in 1:length(directions)){
	
	direction <- directions[j]
	means <- mcols(delta20.gr)[,1:2]
	
	if(direction=="hypo"){
	
	index <- which(means[,2]<means[,1])
	direction.gr <- delta20.gr[index]
	
	}
	else{
	
	index <- which(means[,2]>means[,1])
	direction.gr <- delta20.gr[index]
	
	}
	cat("Direction of methylation change:",direction,"\n")
	previous.progress <- 0
	cat("Calculating the observed and expected distributions...","\n")
	
	for(k in 1:length(repTypes)){

		repType <- repTypes[k]
		index <- which(repFamily==repType)
		repType.gr <- repeats.gr[index]
	
		x <- findOverlaps(repType.gr, allWindows.gr)
		y <- findOverlaps(repType.gr, direction.gr)
	
		exp1 <- c(exp1,length(x))
		obs1<-c(obs1, length(y))
		
		pc <- (k/length(repTypes))*100
		current.progress <- which.min(abs(pc-progress))
			
		if(current.progress>previous.progress){
			
			cat("█")
			
		previous.progress <- current.progress
		
		}
	
	}
	
	exp2 <- c(exp1, length(allWindows.gr)-sum(exp1))
	obs2 <- c(obs1, length(direction.gr)-sum(obs1))

	oe <- (obs2/length(direction.gr)*100)/ (exp2/length(allWindows.gr)*100)

	chi.table <- array(NA,c(length(repTypes),4))
	rownames(chi.table) <- repTypes
	colnames(chi.table) <- c("est","lower","upper","p-val")
	
	previous.progress <- 0
	cat("\n","Performing Chi Square tests...","\n")

	for(l in 1:length(repTypes)){
	
		chi <- matrix(c(obs2[l],exp2[l],length(direction.gr),length(allWindows.gr)),ncol=2,nrow=2)
		chi.table[l,1] <- fisher.test(chi)[[3]]
		chi.table[l,2:3] <- fisher.test(chi)[[2]]
		chi.table[l,4] <- fisher.test(chi)[[1]]
		
		pc <- (l/length(repTypes))*100
		current.progress <- which.min(abs(pc-progress))
			
		if(current.progress>previous.progress){
			
			cat("█")
			
		previous.progress <- current.progress
		
		}

	}
	
	chi.table[,4] <- p.adjust(chi.table[,4],method="BH")
	
	chiName <- paste0(cancer,".",direction,".chiTable")
	assign(chiName,chi.table)
	
	chi.table <- get(chiName)
	
	chi.table2 <- chi.table
	chi.toPlot <- c(which(chi.table[,2]<1 & chi.table[,3]<1),which(chi.table[,2]>1 & chi.table[,3]>1 & chi.table[,4]<0.05))
	chi.toPlot <- sort(chi.toPlot)
		
	chi.table2 <- chi.table[chi.toPlot,]
	chi.table2 <- chi.table2[order(chi.table2[,"est"],decreasing=T),]
	
	plotName <- paste0("Barplot_",cancer,"_",direction,"_repeats_OR.png")
	png(plotName,h=6,w=12,res=300,unit="in")
	par(mar=c(7,5,2,2))
	x <- barplot(chi.table2[,"est"])
	ylim <- c(0,max(chi.table2[,"upper"])*1.1)
	barplot(chi.table2[,"est"],ylim=ylim,names=rownames(chi.table2),las=2,main=cancer,ylab="Odds Ratio")
	arrows(x, chi.table2[,"lower"], x, chi.table2[,"upper"], length=0.05, angle=90, code=3)
	abline(h=1,lty=2,lwd=2,col="red")
	legendText <- paste0("N = ",length(direction.gr)," of ",length(allWindows.gr))
	legend("topright",legend=legendText,bty="n")
	
	dev.off()
	cat("\n")

}
}
