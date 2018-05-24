##
## James' marvellous gene-centric heatmap
##

# Read in the table of CpGs filtered for 10+ read coverage

# Load the required packages
library(GenomicRanges)
library(biomaRt)
library(gplots)

# Read in the CpG-level methylation data
cpgs <- readRDS("/data/SHARE/GINA/pALL4_ovca_150217_filtered10coverage.rds")

# Create a GRanges object from the CpG-level methylation data
cpgs.gr <- GRanges(cpgs[,1],IRanges(cpgs[,2],cpgs[,3]))
mcols(cpgs.gr) <- as.data.frame(cbind(methCase=cpgs[,"meth55"],methControl=cpgs[,"meth56"],deltaMeth=cpgs[,"meth55"]-cpgs[,"meth56"]))

# Get a list of all Ensembl genes in hg19 with their genomic co-ordinates
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
att <- listAttributes(mart)
genes <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position", "end_position","5_utr_start","5_utr_end","3_utr_start","3_utr_end"), mart=mart)
genes[,"chromosome_name"] <- paste0("chr",genes[,"chromosome_name"])

# Get a list of the relevant chromosomes in hg19
chromosomes <- levels(seqnames(cpgs.gr))

# Subset the list of Ensembl genes to just thouse on the relevant chromosomes
genes <- genes[which(genes[,"chromosome_name"] %in% chromosomes),]
# Remove duplicate Ensembl IDs
genes <- genes[which(!duplicated(genes[,"ensembl_gene_id"])),]

# Define the length of the upstream and downstream regions
intergenic.extension <- 5000

# Define the length of the promoter region
promoter <- c(1000,500)

# Create a GRanges object from the list of Ensembl genes with the intergenic extension window
genes.gr <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"],genes[,"end_position"]))
mcols(genes.gr) <- genes[,"ensembl_gene_id"]
genes.extended.gr <- genes.gr+intergenic.extension
mcols(genes.extended.gr) <- genes[,"ensembl_gene_id"]

# Subset the CpG GRanges objects to just those CpGs within the Ensembl genes
cpgs.genes <- subsetByOverlaps(cpgs.gr,genes.extended.gr)
cpgs.genes <- sort(cpgs.genes)
cpgs.genes.counts <- countOverlaps(genes.extended.gr,cpgs.genes)

# Find those genes with at least 100 CpGs
min.100.CpGs <- genes.gr[which(cpgs.genes.counts>=100)]
genes.min.100.CpGs <- genes.extended.gr[cpgs.genes.counts]

# Get the genomic locations of each exon within Ensembl annotated genes
exons <- getBM(attributes=c("ensembl_gene_id","ensembl_exon_id","chromosome_name","start_position", "end_position","exon_chrom_start","exon_chrom_end","rank"),mart=mart)
exons[,"chromosome_name"] <- paste0("chr",exons[,"chromosome_name"])
# Remove weird chromosomes
exons <- exons[which(exons[,"chromosome_name"] %in% chromosomes),]
# Convert to GRanges
exons.gr <- GRanges(exons[,"chromosome_name"],IRanges(exons[,"exon_chrom_start"],exons[,"exon_chrom_end"]))
mcols(exons.gr) <- cbind(exons[,"ensembl_gene_id"],exons[,"rank"])

# Get just exons 1 and 2
exon1.gr <- exons.gr[which(mcols(exons.gr)[[2]]==1)]
exon2.gr <- exons.gr[which(mcols(exons.gr)[[2]]==2)]

# Find which genes contain at least 2 exons
min.2.exons <- mcols(exon1.gr)[[1]][which(mcols(exon1.gr)[[1]] %in% mcols(exon2.gr)[[1]])]

# Infer the genomic coordinates of intron 1 based on exon 1 and 2
ids <- c()
chr <- c()
intron.start <- c()
intron.end <- c()
for(n in 1:1){
	progress <- seq(5,100,5)
	cat("\n","Began working at",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","        ","100","\n","")
    previous.progress <- 0

	for(i in 1:length(min.2.exons)){
		ensembl.id <- min.2.exons[i]
		index <- which(mcols(exon1.gr)[[1]]==ensembl.id)
		ids <- c(ids,as.character(mcols(exon1.gr)[[1]][index]))
		chr <- c(chr,as.character(seqnames(exon1.gr)[index]))
		intron.start <- c(intron.start,as.numeric(end(exon1.gr)[index])+1)
		intron.end <- c(intron.end,as.numeric(start(exon2.gr)[index])-1)
		
		pc <- (i/length(min.2.exons))*100
		current.progress <- which.min(abs(pc-progress))
			if(current.progress>previous.progress){
				cat("█")
			previous.progress <- current.progress
			}
	}
}

intron1.gr <- GRanges(chr,IRanges(intron.start,intron.end))
mcols(intron1.gr) <- ids
# Remove any duplicate Ensembl IDs
intron1.gr <- intron1.gr[which(!duplicated(ids))]

# Get the coordinates of other exons
exonsOther.gr <- exons.gr[which(mcols(exons.gr)[[2]]>1)]

# Get the coordinates of the 5' UTRs
utr5.index <- which(!is.na(genes[,"5_utr_start"])
utr5.gr <- GRanges(genes[utr5.index,"chromosome_name"],IRanges(genes[utr5.index,"5_utr_start"],genes[utr5.index,"5_utr_end"]))
mcols(utr5.gr) <- genes[utr5.index,"ensembl_gene_id"]
utr5.gr <- c(utr5.gr,utr5.inferred)
# Infer the coordinates of the 5'UTRs for those where this information is missing
utr5.inferred.index <- which(exon$rank==1 & !is.na(exon[,"5_utr_start"]) & exon[,"start_position"]!=exon[,"exon_chrom_start"])
utr5.inferred <- GRanges(exon[utr5.inferred.index,"chromosome_name"],IRanges(exon[utr5.inferred.index,"start_position"],exon[utr5.inferred.index,"exon_chrom_start"]))
mcols(utr5.inferred) <- exon[utr5.inferred.index,"ensembl_gene_id"]


# Get the coordinates of the promoter
promoter.gr <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"]-promoter[1],genes[,"start_position"]+promoter[2]))
mcols(promoter.gr) <- genes[,"ensembl_gene_id"]

# Get the coordinates of the 3' UTR
utr3.index <- which(!is.na(exon[,"3_utr_start"]))
utr3.gr <- GRanges(utr3[utr3.index,"chromosome_name"],IRanges(utr3[utr3.index,"3_utr_start"],utr3[utr3.index,"3_utr_end"]))
mcols(utr3.gr) <- utr3[utr3.index,"ensembl_gene_id"]

# Infer the coordinates of other 3' UTRs
# Find the last exon of every gene
# Test whether the end coordinate of the gene is the same as that of the final exon
# If not, find the difference and infer that as the 3' UTR
exons.ranked <- exons[order(exons$rank,descending=T),]
exons.notDup <- exons.ranked[which(!duplicated(exons.ranked$ensembl_gene_id))]
utr3.index <- exons.notDup[which(exons.notDup[,"end_position"]>exons.notDup[,"exon_chrom_end"]),"ensembl_gene_id"]

for(n in 1:1){
	progress <- seq(5,100,5)
	cat("\n","Began working at",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","        ","100","\n","")
    previous.progress <- 0

	for(i in 1:length(utr3.index)){
		ensembl.id <- utr3.index[i]
		index <- which(exon.notDup$ensembl_gene_id==ensembl.id)
		ids <- c(ids,exon.notDup$ensembl_gene_id[index])
		chr <- c(chr,as.character(exon.notDup$chromosome_name[index]))
		exon.end <- c(exon.end,as.numeric(exon.notDup$exon_chrom_end[index])+1)
		gene.end <- c(gene.end,as.numeric(exon.notDup$end_position[index])-1)
		
		pc <- (i/length(min.2.exons))*100
		current.progress <- which.min(abs(pc-progress))
			if(current.progress>previous.progress){
				cat("█")
			previous.progress <- current.progress
			}
	}
}

utr3.inferred <- GRanges(chr,IRanges(exon.end,gene.end))
mcols(utr3.inferred) <- ids
utr3.gr <- c(utr3.gr,utr3.inferred)

# Get the coordinates of other introns
notIntrons <- c(exons.gr,intron1.gr,utr5.gr,utr3.gr)
intronOther <- setdiff(gene.gr,notIntrons.gr)

# Get the coordinates of the upstream and downstream windows
upstream.gr <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"start_position"]-intergenic.extension,genes[,"start_position"]-promoter[1]))
mcols(upstream.gr) <- genes[,"ensembl_gene_id"]
downstream.gr <- GRanges(genes[,"chromosome_name"],IRanges(genes[,"end_position"]+1,genes[,"end_position"]+intergenic.extension))
mcols(downstream.gr) <- genes[,"ensembl_gene_id"]

#genicRegions <- GRangesList(upstream.gr,promoter.gr,utr5.gr,exon1.gr,intron1.gr,exonOther.gr,intronOther.gr,utr3.gr,downstream.gr)

regions <- c("upstream","promoter","utr5","exon1","intron1","exonOther","intronOther","utr3","downstream")

regional.mean.meth <- array(NA,dim=c(length(genes.100.CpGs),length(regions)))
colnames(meth.matrix) <- regions


for(i in 1:length(genes.100.CpGs)){
	gene <- genes.100.CpGs[i]
	ensembl.id <- mcols(gene)[[1]]
		
	for(j in 1:length(regions)){
		
		region <- regions[j]
		region.gr <- get(paste0(region,".gr"))
		
		if(ensembl.id %in% mcols(region.gr)[[1]]){
			index <- which(mcols(region.gr)[[1]]==ensembl.id)[1]
			region.cpgs <- subsetByOverlaps(gene.cpgs,region.gr[index])
			meth.values <- mcols(region.cpgs)[[3]]
			mean.meth <- mean(meth.values)
			regional.mean.meth[i,j] <- mean.meth
		}
	}
}

precision <- 100

progress <- seq(5,100,5)

for(i in 1:length(regions)){
	
	region <- regions[i]
	region.gr <- get(paste0(region,".gr"))
	
	meth.matrix <- array(NA,dim=c(length(region.gr),precision))
	
	
	cat("\n","Began working on ",region," at ",as.character(Sys.time()),
    "\n","Percentage complete",
    "\n","0","        ","100","\n","")
    previous.progress <- 0
	
	for(j in 1:length(region.gr)){
		
		local.region.gr <- region.gr[j]
		local.region.cpgs <- subsetByOverlaps(gene.cpgs,local.region.gr)
		
		meth.values <- mcols(local.region.cpgs)[[3]]
		
		if(length(meth.values)==precision){
			meth.matrix[j,] <- meth.values
		} else {
			local.region.length <- width(local.region.gr)
			new.meth.values <- sapply(split(meth.values,cut(1:length(meth.values),precision)),mean)
			meth.matrix[j,] <- new.meth.values
		}
		
		pc <- (i/length(region.gr))*100
		current.progress <- which.min(abs(pc-progress))
			if(current.progress>previous.progress){
				cat("█")
			previous.progress <- current.progress
			}
		
	}
	
	r.name <- paste0(region,".methMatrix")
	assign(matrix.name,meth.matrix)
	rds.name <- paste0(r.name,".RDS")
	saveRDS(meth.matrix,rds.name)

}

for(i in 1:length(regions)){
	
	region <- regions[i]
	meth.matrix <- paste0(region,".methMatrix")
	meth.matrix <- get(meth.matrix)
	
	pngName <- paste0("20180521_heatmap_deltaMeth_",region,".png")
	png(pngName,h=12,w=12,res=300,unit="in")
	heatmap.2(meth.matrix,trace="n",col=bluered(100),cexRow="",cexCol="",Colv=NA)
	dev.off()
	
}