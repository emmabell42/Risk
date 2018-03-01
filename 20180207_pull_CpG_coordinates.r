# This code gets the genomic co-ordinates of CpGs within your region of interest

# Load the required R packages
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
library(stringr)

# Replace "XXXXXXX" with the IDs of the window you're interested in
candidates <- c("chr9_6948600","chr4_92392600")

# Get the co-ordinates of the 200 bp window from their IDs
chr <- unlist(lapply(str_split(candidates,"_"),`[[`, 1))
starts <- as.numeric(unlist(lapply(str_split(candidates,"_"),`[[`, 2)))
ends <- starts+200

# Get the DNA sequences of each 200 bp window
candidates.gr <- GRanges(chr,IRanges(starts,ends))
seqs <- getSeq(Hsapiens, seqnames(candidates.gr), start(candidates.gr), end(candidates.gr), as.character = T)

# Record the position of each CpG in each 200 bp window
cpg.coordinates <- c()
for(i in 1:length(seqs)){
	positions <- gregexpr("CG",seqs[i])[[1]]
	cpgs.inWindow <- array(NA,dim=c(length(positions),4))
	for(j in 1:length(positions)){
		pos <- positions[j]
		cpgs.inWindow[,1] <- candidates[i]
		cpgs.inWindow[j,2] <- chr[i]
		cpgs.inWindow[j,3] <- starts[i]+pos-1
		cpgs.inWindow[j,4] <- starts[i]+pos
	}
cpg.coordinates <- rbind(cpg.coordinates,cpgs.inWindow)
}
colnames(cpg.coordinates) <- c("ID","CpG_chr","CpG_start","CpG_end")

# Write out the table of CpG co-ordinates

write.table(cpg.coordinates,"CpG_coordinates.txt",sep="\t",quote=F)