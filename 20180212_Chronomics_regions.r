# This code creates a dummy list of candidate genomic regions
# to send to Chronomics, so I can get an idea of how well they
# cover my regions of interest

library(GenomicRanges)
library(stringr)

candidates <- read.table("OvCa_0_delta20.candidates_methylationSummaryScore.txt",sep="\t",head=T,stringsAsFactors=F)

candidates <- rownames(candidates)
chr <- unlist(lapply(str_split(candidates,"_"),`[[`, 1))
starts <- as.numeric(unlist(lapply(str_split(candidates,"_"),`[[`, 2)))
ends <- starts+200

candidates <- cbind(chr,starts,ends)
write.table(candidates,"candidates.bed",sep="\t",col.names=F,row.names=F,quote=F)

candidates.gr <- GRanges(chr,IRanges(starts,ends))
