# This code compares the GC% and normalised CpG fraction in
# unique and overlapping delta20 windows

setwd("/data/emmabell42/risk/BrCa_x_OvCa")

# Load packages

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
library(stringr)

# Read in files
rdsFiles <- list.files("..")[grep("20.RDS",list.files(".."))]
for(i in 1:length(rdsFiles)){
toRead <- paste0("../",rdsFiles[i])
rds <- readRDS(toRead)
rds <- GRangesList(rds)
rds <- unlist(rds)
#rds <- reduce(rds) # added to join adjacent windows
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

# Retrieve the DNA sequences of your regions 

cancerTypes <- c("BrCa","OvCa")

for(i in 1:length(cancerTypes)){

	cancer <- cancerTypes[i]
	gr <- paste0(cancer,"_0_delta20.gr")
	gr <- get(gr)
	
	seqs <- getSeq(Hsapiens, seqnames(gr), start(gr), end(gr), as.character = T)
	seqs <- DNAStringSet(x=seqs,use.names=T)
	
	toName <- paste0(cancer,".seqs")
	assign(toName,seqs)
	
}

#calculate CpG, C and G frequency
progress <- seq(5,100,5)
for(i in 1:length(cancerTypes)){

	cancer <- cancerTypes[i]
	seqs <- paste0(cancer,".seqs")
	seqs <- get(seqs)
	lengths <- width(seqs)
	
	cat("\n","Began working on",cancer,"at",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","        ","100","\n","")
    previous.progress <- 0

	cg.count <- rep(NA,length(seqs))
	c.count <- rep(NA,length(seqs))
	g.count <- rep(NA,length(seqs))

	for(j in 1:length(seqs)){
		cg.count[j] <- str_count(seqs[j], fixed("CG"))
		c.count[j] <- str_count(seqs[j], fixed("C"))
		g.count[j] <- str_count(seqs[j], fixed("G"))
		
		pc <- (j/length(seqs))*100
        current.progress <- which.min(abs(pc-progress))
        if(current.progress>previous.progress){
            cat("█")
        previous.progress <- current.progress
		}
	
	}

	#calculate GC%
	GC <- (c.count+g.count)/lengths
	
	toName <- paste0(cancer,".gc")
	assign(toName,GC)

	#calculate the normalised fraction of CpGs
	#Saxonov method
	oe <- (cg.count/(lengths))/((GC/2)^2)
	
	toName <- paste0(cancer,".oe")
	assign(toName,oe)

}
cat("\n","Done!","\n")



