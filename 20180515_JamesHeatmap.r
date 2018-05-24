##
## James' marvellous gene-centric heatmap
##

# Read in the table of CpGs filtered for 10+ read coverage

library(GenomicRanges)

cpgs <- readRDS("/data/SHARE/GINA/pALL4_ovca_150217_filtered10coverage.rds")

cpgs.gr <- GRanges(cpgs[,1],IRanges(cpgs[,2],cpgs[,3]))
mcols(cpgs.gr) <- as.data.frame(cbind(methCase=cpgs[,4],cpgsCase=cpgs[,7],methControl=cpgs[,11],cpgsControl=cpgs[,14]))

rm(cpgs)

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
att <- listAttributes(mart)
genes <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position", "end_position","5_utr_start","5_utr_end","3_utr_start","3_utr_end"), mart=mart)
genes[,"chromosome_name"] <- paste0("chr",genes[,"chromosome_name"])

chromosomes <- levels(seqnames(cpgs.gr))

genes <- genes[which(genes[,"chromosome_name"] %in% chromosomes),]
genes <- na.omit(genes)

intergenic.extension <- 5000

genes.gr <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"],genes[,"end_position"]))
mcols(genes.gr) <- genes[,"ensembl_gene_id"]
genes.extended <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"]-intergenic.extension,genes[,"end_position"]+intergenic.extension))

cpgs.genes <- subsetByOverlaps(cpgs.gr,genes.extended)

exons <- getBM(attributes=c("ensembl_gene_id","ensembl_exon_id","chromosome_name","start_position", "end_position","exon_chrom_start","exon_chrom_end","rank"),mart=mart)
exons[,"chromosome_name"] <- paste0("chr",exons[,"chromosome_name"])

exons <- exons[which(exons[,"chromosome_name"] %in% chromosomes),]

exon.1 <- exons[which(exons[,"rank"]==1),]
exons.not1 <- exons[which(exons[,"rank"]>1),]

upstream <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"]-integenic.extension))
utr5 <- GRanges(genes[,"chromosome_name"],IRanges(
