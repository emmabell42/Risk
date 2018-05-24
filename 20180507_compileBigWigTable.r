# This script is part of the Hyperactive Enhancer project. It:
# 	1. Reads in the methylation score and read count BigWig files;
#	2. Compiles a table of methylation scores and read counts.

setwd("/data/emmabell42/Hyperactive enhancers/MethBase")


# Read in PU and DM
library(GenomicRanges)
pu <- read.table("../../resources/protectedregions.txt")
dm <- read.table("../../resources/DMRregions.txt")
enhancers <- rbind(pu,dm)
pu <- GRanges(seqnames=pu[,1],ranges=IRanges(pu[,2],pu[,3]))
dm <- GRanges(seqnames=dm[,1],ranges=IRanges(dm[,2],dm[,3]))
enhancers <- GRanges(seqnames=enhancers[,1],ranges=IRanges(enhancers[,2],enhancers[,3]))

# Read in BigWig files 
# Subset for CpGs with with 10+ reads
# Subset for CpGs that overlaps with my enhancers
# Get a list of all CpGs within those enhancers
# Create a GRangesList from the methylation GRanges objects
library(rtracklayer)
bws <- list.files(pattern=".bw",recursive=T)
toReplace <- c(".read.bw",".meth.bw")
paths <- gsub(paste(toReplace,collapse="|"),"",bws)
paths <- paths[!duplicated(paths)]
rNames <- gsub("-","_",paths)
rNames <- gsub("\\/",".",rNames)

all.cpgs <- c()
good.cpgs <- c()

progress <- seq(5,100,5)
previous.progress <- 0
for(i in 1:length(paths)){
	studyToRead <- paths[i]
	rName <- rNames[i]
	
	if(i==1){
	cat("\n","Began working at",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","                ","100","\n","")
	}
    	
	methToRead <- paste0(studyToRead,".meth.bw")
	coverageToRead <- paste0(studyToRead,".read.bw")
	
	meth <- import(con=methToRead,format="bw")
	read <- import(con=coverageToRead,format="bw")
	
	coverage <- mcols(read)[[1]]
	goodCoverage <- which(coverage>=10)
	
	goodMeth <- meth[goodCoverage]
	
	#chr <- seqnames(goodMeth)
	#starts <- start(goodMeth)
	#all.cpgs_id <- paste(chr,starts,sep="_")
	#all.cpgs <- c(all.cpgs,all.cpgs_id)
	
	goodMethEnhancers <- subsetByOverlaps(goodMeth,enhancers)
	
	chr <- seqnames(goodMethEnhancers)
	starts <- start(goodMethEnhancers)
	good.cpgs_id <- paste(chr,starts,sep="_")
	good.cpgs <- c(good.cpgs,good.cpgs_id)
		
	assign(rName,goodMethEnhancers)
	
	if(i==1){
		#allMethList <- list(meth)
		#goodMethList <- list(goodMeth)
		enhancerMethList <- list(goodMethEnhancers)
	} else {
		#allMethList <- c(allMethList,meth)
		#goodMethList <- c(goodMethList,goodMeth)
		enhancerMethList <- c(enhancerMethList,goodMethEnhancers)
	}
	
	pc <- (i/length(paths))*100
	current.progress <- which.min(abs(pc-progress))
    if(current.progress>previous.progress){
		cat("█")
	previous.progress <- current.progress
	}
	
}

#cat("\n","Began working on allMeth",as.character(Sys.time()))
#png("Boxplot_allMeth.png",h=6,w=12,unit="in",res=300)
#boxplot(lapply(allMethList,function(x) mcols(x)[[1]]),col="lightgrey",xlab="Data set",ylab="CpG methylation",pch=20)
#dev.off()
#cat("\n","Began working on goodMeth",as.character(Sys.time()))
#png("Boxplot_goodMeth.png",h=6,w=12,unit="in",res=300)
#boxplot(goodMethList,col="lightgrey",xlab="Data set",ylab="CpG methylation",pch=20)
#dev.off()

metadata <- read.csv("MethBase_index.csv")
col <- rainbow(length(table(metadata$Source)))

png("Boxplot_enhancerMeth.png",h=6,w=12,unit="in",res=300)
boxplot(lapply(enhancerMethList,function(x) mcols(x)[[1]]),col=col[as.numeric(metadata$Source)],xlab="Data set",ylab="CpG methylation",pch=19)
dev.off()

for(i in 1:length(paths)){
	dataset <- enhancerMethList[[i]]
	hypo <- subsetByOverlaps(dataset,pu)
	hyper <- subsetByOverlaps(dataset,dm)
	rName <- rNames[i]
	toName <- strsplit(rName,"\\.")[[1]][2]
	fileName <- paste0("Boxplot_",rName,"_hypoVsHyper.png")
	png(fileName,h=6,w=3,res=300,unit="in")
	boxplot(mcols(hypo)[[1]],mcols(hyper)[[1]],col=c("green","magenta"),main=toName,ylim=c(0,1),pch=19,ylab="mCpG/CpG",names=c("PU","DM"))
	dev.off()
}



