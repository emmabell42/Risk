##
## James' marvellous gene-centric heatmap
##

# Read in the table of CpGs filtered for 10+ read coverage

library(GenomicRanges)

cpgs <- readRDS("/data/SHARE/GINA/pALL4_ovca_150217_filtered10coverage.rds")

cpgs.gr <- GRanges(cpgs[,1],IRanges(cpgs[,2],cpgs[,3]))
mcols(cpgs.gr) <- as.data.frame(cbind(methCase=cpgs[,"meth55"],methControl=cpgs[,"meth56"],deltaMeth=cpgs[,"meth55"]-cpgs[,"meth56"]))

library(biomaRt)

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
att <- listAttributes(mart)
genes <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position", "end_position","5_utr_start","5_utr_end","3_utr_start","3_utr_end"), mart=mart)
genes[,"chromosome_name"] <- paste0("chr",genes[,"chromosome_name"])

chromosomes <- levels(seqnames(cpgs.gr))

genes <- genes[which(genes[,"chromosome_name"] %in% chromosomes),]
genes <- genes[which(!duplicated(genes[,"ensembl_gene_id"])),]

intergenic.extension <- 5000

genes.gr <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"],genes[,"end_position"]))
mcols(genes.gr) <- genes[,"ensembl_gene_id"]

genes.extended.gr <- genes.gr+intergenic.extension

cpgs.genes <- subsetByOverlaps(cpgs.gr,genes.extended)

genebody.cpgs <- subsetByOverlaps(cpgs.gr,genes.gr)
genebody.cpgs <- sort(genebody.cpgs)
cpgs.in.genebody <- countOverlaps(genes.gr,genebody.cpgs)

genes.100.CpGs <- genes.gr[which(cpgs.in.genebody>=100)]

genes.100.CpGs.upstream <- GRanges(seqnames(genes.100.CpGs),IRanges(start(genes.100.CpGs)-intergenic.extension,start(genes.100.CpGs)-1))
genes.100.CpGs.downstream <- GRanges(seqnames(genes.100.CpGs),IRanges(end(genes.100.CpGs)+1,end(genes.100.CpGs)+intergenic.extension))

precision <- 100

meth.matrix <- array(NA,dim=c(length(genes.100.CpGs),precision))

progress <- seq(5,100,5)
previous.progress <- 0

for(n in 1){
cat("\n","Began working at",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","                ","100","\n","")
for(i in 1:length(genes.100.CpGs)){
	
	pc <- (i/length(genes.100.CpGs))*100
	current.progress <- which.min(abs(pc-progress))
	if(current.progress>previous.progress){
		cat("█")
	previous.progress <- current.progress
	}

	gene <- genes.100.CpGs[i]
	cpgs <- subsetByOverlaps(genebody.cpgs,gene)
	
	meth.values <- mcols(cpgs)[[3]]
	
	if(length(meth.values)==precision){
		meth.matrix[i,] <- meth.values
	} else {
		gene.length <- width(gene)
		new.meth.values <- sapply(split(meth.values,cut(1:length(meth.values),precision)),mean)
		meth.matrix[i,] <- new.meth.values
	}
}
}

saveRDS(meth.matrix,"methMatrix_100bins.RDS")

library(gplots)

rownames(meth.matrix) <- mcols(genes.100.CpGs)[[1]]
breaks = c(seq(min(meth.matrix),-10,length=100),seq(-9,9,length=100),seq(10,max(meth.matrix),length=100))
cols <- colorRampPalette(c("red", "white", "blue"))(n = 299)

rowsToPlot <- nrow(meth.matrix)

#cor.m <- cor(t(meth.matrix[1:rowsToPlot,]),method="spearman")
#dendro <- hclust(as.dist(1-cor.m))

png("James_marvellous_heatmap.png",h=12,w=12,res=300,unit="in")
heatmap.2(meth.matrix[1:rowsToPlot,],col=cols,breaks=breaks,symm=F,symkey=F,symbreaks=T, scale="none",Colv=F,trace="n",
#Rowv=as.dendrogram(dendro),
labRow="",labCol="")
dev.off()
