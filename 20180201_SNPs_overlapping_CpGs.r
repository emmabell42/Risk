# Check for SNPs at the pyro candidate CpGs

library(biomaRt)
snp.db <- useMart("ENSEMBL_MART_SNP",host="grch37.ensembl.org", dataset="hsapiens_snp")

pyro.snps <- as.list(rep(NA,nrow(OvCa_0_delta20.pyro)))

for(i in 1:nrow(OvCa_0_delta20.pyro)){
cat(i,"\n")
chr <- as.character(OvCa_0_delta20.pyro[i,2])
chr <- gsub("chr","",chr)
	if(chr %in% 1:22) {
	chr <- as.numeric(chr)
	}
starts <- OvCa_0_delta20.pyro[i,3]
ends <- OvCa_0_delta20.pyro[i,4]
snps <- getBM(c("chr_name","chrom_start","chrom_end",
"refsnp_id","allele","minor_allele_freq","validated"),
filters=c("chr_name","start","end"),
values=list(chr_name=chr,chrom_start=starts,chrom_end=ends),
mart=snp.db)
pyro.snps[[i]] <- snps
}

lengths <- OvCa_0_delta20.pyro[,4]-OvCa_0_delta20.pyro[,3]
snp.rate <- c()

for(i in 1:nrow(OvCa_0_delta20.pyro)){
snps <- pyro.snps[[i]]
snp.rate <- c(snp.rate,nrow(snps)/lengths[i])
}

rm(lengths)

# Find CpGs that overlap with SNPs

cpg.snps <- as.list(rep(NA,nrow(OvCa_0_delta20.pyro)))
snpsPerWindow <- c()

for(i in 1:nrow(OvCa_0_delta20.pyro)){

snps.candidate <- pyro.snps[[i]]
snps.candidate.gr <- GRanges(paste0("chr",snps.candidate[,1]),IRanges(snps.candidate[,2],snps.candidate[,3]))
mcols(snps.candidate.gr) <- snps.candidate[,4:ncol(snps.candidate)]

pyro.candidate.gr <- GRanges(OvCa_0_delta20.pyro[i,2],IRanges(OvCa_0_delta20.pyro[i,3],OvCa_0_delta20.pyro[i,4]))

cpgs.candidate <- subsetByOverlaps(coverage.gr,pyro.candidate.gr)
cpg.snps[[i]] <- subsetByOverlaps(snps.candidate.gr,cpgs.candidate)
cpg.snps[[i]] <- cpg.snps[[i]][which(!cpg.snps[[i]]$validated=="")]

snpsPerWindow <- c(snpsPerWindow,length(cpgs.candidate)-length(cpg.snps[[i]]))

}
