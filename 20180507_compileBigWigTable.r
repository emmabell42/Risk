# This script is part of the Hyperactive Enhancer project. It:
# 	1. Reads in the methylation score and read count BigWig files;
#	2. Compiles a table of methylation scores and read counts.

setwd("/data/emmabell42/H3K4me3/BigWigs")

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

cpgs <- c()

progress <- seq(5,100,5)
for(i in 1:length(paths)){
	studyToRead <- paths[i]
	rName <- rNames[i]
	
	if(i==1){
	cat("\n","Began working at",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","                ","100","\n","")
	}
    previous.progress <- 0
	
	methToRead <- paste0(studyToRead,".meth.bw")
	coverageToRead <- paste0(studyToRead,".read.bw")
	
	meth <- import(con=methToRead,format="bw")
	read <- import(con=coverageToRead,format="bw")
	
	coverage <- mcols(read)[[1]]
	goodCoverage <- which(coverage>=10)
	
	goodMeth <- meth[goodCoverage]
	
	goodMethEnhancers <- subsetByOverlaps(goodMeth,enhancers)
	
	chr <- seqnames(goodMethEnhancers)
	starts <- start(goodMethEnhancers)
	cpgs_id <- paste(chr,starts,sep="_")
	cpgs <- c(cpgs,cpgs_id)
		
	assign(rName,goodMethEnhancers)
	
	if(i==1){
		allMethList <- list(meth)
		readList <- list(read)
		
		methList <- list(goodMethEnhancers)
	} else {
		allMethList <- c(allMethList,meth)
		readList <- c(readList,read)
		
		methList <- c(methList,goodMethEnhancers)
	}
	
	pc <- (i/length(paths))*100
	current.progress <- which.min(abs(pc-progress))
    if(current.progress>previous.progress){
		cat("█")
	previous.progress <- current.progress
	}
	
}
cpgs.unique <- cpgs[!duplicated(cpgs)]

allMethList <- GRangesList(allMethList)
readList <- GRangesList(readList)
methList <- GRangesList(methList)
names(allMethList) <- rNames
names(readList) <- rNames
names(methList) <- rNames

# Save the methylation GRangesList as an RDS
saveRDS(allMethList,"MethBase_mouse_methylation.RDS")
saveRDS(readList,"MethBase_mouse_coverage.RDS")
saveRDS(methList,"MethBase_mouse_enhancerMethylation.RDS")
