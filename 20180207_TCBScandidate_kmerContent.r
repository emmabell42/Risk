# Assessing the AAAAA and TTTTT of candidate TCBS probes

candidates <- read.table("OvCa_0_delta20.candidates_methylationSummaryScore.txt",head=T,sep="\t",stringsAsFactors=F)

library(stringr)
library(GenomicRanges)

chr <- unlist(lapply(str_split(rownames(candidates),"_"),`[[`,1))
starts <- as.numeric(unlist(lapply(str_split(rownames(candidates),"_"),`[[`,2)))

candidates.gr <- GRanges(chr,IRanges(starts,starts+200))
mcols(candidates.gr) <- data.frame(ids=rownames(candidates))

library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
library(Biostrings)

seqs <- getSeq(Hsapiens, seqnames(candidates.gr), start(candidates.gr), end(candidates.gr), as.character = T)
seqs <- DNAStringSet(x=seqs,use.names=T)

oligos <- c()

width <- 5
for(i in 1:length(seqs)){
oligos <- rbind(oligos,oligonucleotideFrequency(seqs[[i]],width=width))
}

baduns <- oligos[,1]+oligos[,1024]
candidates <- cbind(candidates,"AAAAA/TTTTT"=baduns)
write.table(candidates,"OvCa_0_delta20.candidates_methylationSummaryScore.txt",sep="\t",quote=F)