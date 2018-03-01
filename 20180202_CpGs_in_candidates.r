

# Sequence composition analysis

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
library(stringr)

candidates <- read.table("OvCa_0_delta20.candidates_methylationSummaryScore.txt",sep="\t",head=T,stringsAsFactors=F)
chr <- unlist(lapply(str_split(rownames(candidates),"_"),`[[`, 1))
starts <- as.numeric(unlist(lapply(str_split(rownames(candidates),"_"),`[[`, 2)))
ends <- starts+200

coverage <- readRDS("/data/SHARE/GINA/pALL4_ovca_270317.rds")

coverage.gr <- GRanges(coverage[,1],IRanges(coverage[,2],coverage[,3]))
mcols(coverage.gr) <- as.data.frame(cbind(methCase=coverage[,4],coverageCase=coverage[,7],methControl=coverage[,11],coverageControl=coverage[,14]))

ids <- paste(chr,starts,sep="_")

OvCa_0_delta20.candidates.gr <- GRanges(chr,IRanges(starts,ends))
mcols(OvCa_0_delta20.candidates.gr) <- data.frame(ids)

coverage.subset.gr <- subsetByOverlaps(coverage.gr,OvCa_0_delta20.candidates.gr)

coverage.subset.gr <- coverage.subset.gr[which(mcols(coverage.subset.gr)$coverageCase>=10 & mcols(coverage.subset.gr)$coverageControl>=10)]

coverage.subset.df <- cbind(seqnames(coverage.subset.gr),ranges(coverage.subset.gr),mcols(coverage.subset.gr))


