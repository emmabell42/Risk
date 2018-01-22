# This code compiles the pyrosequencing candidate table
# Candidates must have a high number of CpGs, high difference in methylation, not occur in repetitive genomic regions, and have coverage of 10+.

cancerType <- c("OvCa")
min.diff <- 0

gr <- list.files()[grep("RDS",list.files())]
for(i in 1:length(gr)){
	reading <- readRDS(gr[i])
	toName <- paste0(cancerType[i],".diff.gr")
	assign(toName,reading)
}

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
library(stringr)

for(i in 1:length(cancerType)){

	type <- cancerType[i]
	gr <- get(paste0(cancerType[i],".diff.gr"))
		
	cat("\n","Analysing cancer type: ",type)
	
	for(j in 1:23){
	gr.chr <- gr[[j]]
	mcols(gr.chr) <- data.frame(mcols(gr.chr))[,1:6]
	seqs <- getSeq(Hsapiens, seqnames(gr.chr), start(gr.chr), end(gr.chr), as.character = T)

	# Calculate CpG, C and G frequency

	cg.count <- rep(NA,length(gr.chr))
	
	for(k in 1:length(gr.chr)){
		cg.count[k] <- str_count(seqs[k], fixed("CG"))
	}
	
	mcols(gr.chr) <- cbind(mcols(gr.chr),data.frame(CpGs=cg.count))
	
	gr[[j]] <- gr.chr
	}
	
	# Write to R object
	toName <- paste0(cancerType[i],".diff.gr")
	assign(toName,gr)
	
}
cat("\n")

# Read in HOMER annotations

annotations <- list.files()[grep("homer",list.files())]
for(i in 1:length(annotations)){
	reading <- read.table(annotations[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
	reading$Annotation[grep("intron",reading$Annotation)] <- "intron"
	reading$Annotation[grep("promoter-TSS",reading$Annotation)] <- "promoter-TSS"
	reading$Annotation[grep("exon",reading$Annotation)] <- "exon"
	reading$Annotation[grep("non-coding",reading$Annotation)] <- "non-coding"
	reading$Annotation[grep("TTS",reading$Annotation)] <- "TTS"
	cat(annotations[i]," ",table(reading$Annotation),"\n")
	assign(annotations[i],reading)
}

for(i in 1:length(cancerType)){

	type <- cancerType[i]
	gr <- get(paste0(cancerType[i],".diff.gr"))
	df.diff <- c() 
	
	for(j in 1:23){
		df.diff <- rbind(df.diff,data.frame(seqnames=seqnames(gr[[j]]),starts=start(gr[[j]])-1,ends=end(gr[[j]]),mcols(gr[[j]])))
		
	}
	
	id <- paste(df.diff[,1],df.diff[,2],sep="_")
	df.diff <- cbind(df.diff,id)
	
	toName <- paste0(cancerType[i],".diff.df")
	assign(toName,df.diff)
		
}

colnames(homer_OvCa_0_delta20_CI.txt)[1] <- "id"
OvCa.diff.df <- merge(OvCa.diff.df,homer_OvCa_0_delta20_CI.txt,by="id")

# Exclude windows residing within repetitive elements (i.e. SINE, LINE, LTRs)

OvCa.diff.df <- cbind(OvCa.diff.df,Repetitive=0)
reps <- c("\\|LINE\\|","\\|SINE\\|","\\|LTR\\|")
OvCa.diff.df[grep(paste0(reps,collapse="|"),OvCa.diff.df[,"Detailed.Annotation"]),"Repetitive"] <- 1

OvCa.pyro.candidates <- OvCa.diff.df[,c(1:12,34)]
OvCa.pyro.candidates <- OvCa.pyro.candidates[which(OvCa.pyro.candidates[,"Casemean"]<OvCa.pyro.candidates[,"Controlmean"]),]
OvCa.pyro.candidates <- OvCa.pyro.candidates[which(OvCa.pyro.candidates[,"Repetitive"]==0),]
OvCa.pyro.candidates <- OvCa.pyro.candidates[which(OvCa.pyro.candidates[,"CpGs"]>4),]

# Write out table candidates at text file

write.table(OvCa.pyro.candidates,"20180122_OvCaPyroCandidates.txt",sep="\t",row.names=F,quote=F)

# Assess the read coverage of each CpG within the pyro candidate windows

pyro <- OvCa.pyro.candidates

# Read in the table of CpGs filtered for 10+ read coverage

OvCa <- readRDS("/data/SHARE/GINA/pALL4_ovca_150217_filtered10coverage.rds")

OvCa.gr <- GRanges(OvCa[,1],IRanges(OvCa[,2],OvCa[,3]))
mcols(OvCa.gr) <- as.data.frame(cbind(cover55=OvCa[,6],cover56=OvCa[,13]))
pyro.gr <- GRanges(pyro[,2],IRanges(pyro[,3],pyro[,4]))

# Find the pyro candidate regions that contain CpGs with 10+ read coverage

pyro.subset.gr <- subsetByOverlaps(pyro.gr,OvCa.gr)

# Find the CpGs within pyro candidate regions that have 10+ read coverage

pyro.cpgs.gr <- subsetByOverlaps(OvCa.gr,pyro.gr)

pyro.cpgs.df <- data.frame(seqnames=seqnames(pyro.cpgs.gr),starts=start(pyro.cpgs.gr)-1,ends=end(pyro.cpgs.gr),mcols(pyro.cpgs.gr))
pyro.cpgs.df <- cbind(pyro.cpgs.df,id=NA)

# Add the ID to the candidate CpGs

for(i in 1:nrow(pyro.cpgs.df)){
cpg_chr <- pyro.cpgs.df[i,1]
cpg_start <- pyro.cpgs.df[i,2]
cpg_end <- pyro.cpgs.df[i,3]

candidate_match <- which(pyro[,2]==cpg_chr & cpg_start>=pyro[,3] & pyro[,4]>=cpg_end)

id <- pyro[candidate_match,1]

pyro.cpgs.df[i,"id"] <- id
}

write.table(pyro.cpgs.df,"20180122_OvCaPyroCandidates_CpGcoverage.txt",sep="\t",row.names=F,quote=F)
