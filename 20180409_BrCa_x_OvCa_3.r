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

mcols(BrCa_0_delta20.gr) <- cbind(mcols(BrCa_0_delta20.gr),CpG_norm=BrCa.oe,GC=BrCa.gc)
mcols(OvCa_0_delta20.gr) <- cbind(mcols(OvCa_0_delta20.gr),CpG_norm=OvCa.oe,GC=OvCa.gc)

# Generate "common" and "unique" objects
cancerTypes <- c("BrCa","OvCa")
patientTypes <- c("cases","controls")
regionTypes <- c("common","unique")
#deltaMethDirections <- c("hypo","hyper")

common.BrCa <- subsetByOverlaps(BrCa_0_delta20.gr,OvCa_0_delta20.gr)
common.BrCa <- sort(common.BrCa)
ids.BrCa <- paste(seqnames(common.BrCa),start(common.BrCa),sep="_")
common.OvCa <- subsetByOverlaps(OvCa_0_delta20.gr,BrCa_0_delta20.gr)
common.OvCa <- sort(common.OvCa)
ids.OvCa <- paste(seqnames(common.OvCa),start(common.OvCa),sep="_")
identical(ids.BrCa,ids.OvCa)

BrCa.direction <- mcols(common.BrCa)[,2]-mcols(common.BrCa)[,1]
OvCa.direction <- mcols(common.OvCa)[,2]-mcols(common.OvCa)[,1]

hypo <- which(BrCa.direction<0 & OvCa.direction<0)
hyper <- which(BrCa.direction>0 & OvCa.direction>0)
index <- sort(c(hypo,hyper))

#common.BrCa <- common.BrCa[index]
#common.OvCa <- common.OvCa[index]

# Create function - custom vioplot

for(i in 1:length(cancerTypes)){

	cancer <- cancerTypes[i]
	cancer.gr <- paste0(cancer,"_0_delta20.gr")
	cancer.gr <- get(cancer.gr)
	
	common <- paste0("common.",cancer)
	common <- get(common)
	
	unique. <- cancer.gr[!(cancer.gr %over% common)]
	
	toName <- paste0("unique.",cancer)
	assign(toName,unique.)
}

common <- mcols(common.BrCa)[[7]]	
unique.both <- c(mcols(unique.BrCa)[[7]],mcols(unique.OvCa)[[7]])

pngName <- paste0("Vioplot_CpGnorm_commonVsUnique.png")
png(pngName,w=6,h=6,unit="in",res=300)
vioplot(common,unique.both,col="grey",names=c("Common","Unique"))
title(main="CpG content",xlab="Region type",ylab="Norm. CpG fraction")
dev.off()

unique.b <- mcols(unique.BrCa)[[7]]
unique.o <- mcols(unique.OvCa)[[7]]

pngName <- paste0("Vioplot_CpGnorm_commonVsUnique_byCancer.png")
png(pngName,w=6,h=6,unit="in",res=300)
vioplot(common,unique.b,unique.o,col="grey",names=c("Common","Unique (BrCa)","Unique (OvCa)"))
title(main="CpG content",xlab="Region type",ylab="Norm. CpG fraction")
dev.off()

pngName <- paste0("Vioplot_CpGnorm_commonVsUnique_byCancer_log.png")
png(pngName,w=6,h=6,unit="in",res=300)
vioplot(log10(common),log10(unique.b),log10(unique.o),col="grey",names=c("Common","Unique (BrCa)","Unique (OvCa)"))
title(main="CpG content",xlab="Region type",ylab="Log10(Norm. CpG fraction)")
dev.off()

common <- mcols(common.BrCa)[[8]]	
unique.both <- c(mcols(unique.BrCa)[[8]],mcols(unique.OvCa)[[8]])

pngName <- paste0("Vioplot_GC_commonVsUnique.png")
png(pngName,w=6,h=6,unit="in",res=300)
vioplot(common,unique.both,col="grey",names=c("Common","Unique"))
title(main="GC content",xlab="Region type",ylab="Prop. GC")
dev.off()

unique.b <- mcols(unique.BrCa)[[8]]
unique.o <- mcols(unique.OvCa)[[8]]

pngName <- paste0("Vioplot_GC_commonVsUnique_byCancer.png")
png(pngName,w=6,h=6,unit="in",res=300)
vioplot(common,unique.b,unique.o,col="grey",names=c("Common","Unique (BrCa)","Unique (OvCa)"))
title(main="GC content",xlab="Region type",ylab="Prop GC")
dev.off()

pngName <- paste0("Vioplot_GC_commonVsUnique_byCancer_log.png")
png(pngName,w=6,h=6,unit="in",res=300)
vioplot(log10(common),log10(unique.b),log10(unique.o),col="grey",names=c("Common","Unique (BrCa)","Unique (OvCa)"))
title(main="GC content",xlab="Region type",ylab="Log10(Prop. GC)")
dev.off()

