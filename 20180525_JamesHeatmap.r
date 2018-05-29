# ---
# title: "James' Marvellous Gene-Centered Heatmap"
# output:
  # html_document:
    # df_print: paged
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# ---
# # Aim

# The project aims to create a gene-centered heatmap illustrating the difference in CpG methylation between ovarian cancer cases and controls along the length of the gene. I will initially consider the following regions:
# -Promoter
# -5' UTR
# -Exon 1
# -Intron 1

# ## Requirements

# R packages:

# ```{r echo=TRUE}
 library(GenomicRanges)
 library(biomaRt)
 library(gplots)
# ```

# Data:
# - CpG-level methylation from ovarian cancer cases and controls;
# - Coordinates of all Ensembl annotated genes in hg19.

# ## Output

# Working directory:

# ```{r echo=TRUE}
 setwd("/data/emmabell42/risk")
# ```

# # Methodology

# 1. Prepare the CpG methylation data.
# 	1.1. Read in the CpG-level methylation data pertaining to ovarian cancer cases and controls.
# 	1.2. Convert the table from 1.1 to a GRanges object.
# 2. Get a list of all Ensembl genes in hg19 with the genomic co-ordinates of the TSS and TES.
# 	2.1. Use biomaRt to retrieve the Ensembl ID, TSS, and TES of all genes in hg19.
# 	2.2. Remove entries on non-canonical chromosomes, ncRNAs, and duplicate Ensembl IDs.
# 	2.3. Define an upstream and downstream window for each gene.
# 	2.4. Subset the methylation data to just those CpGs within that overlap with the upstream/downstream windows or overlap a gene body.
# 3. Create a matrix of methylation values centred around the TSS of each gene.
# 	3.1. For a given Ensembl ID, get the delta methylation scores for all CpGs within all regions.
# 	3.2. For a given Ensembl ID, average the methylation scores from 3.1.
# 	3.3. Plot the resulting average methylation scores for each region set.

# # Results

# 1. Prepare the CpG-level methylation data.
# 	1.1. Read in the CpG methylation data pertaining to ovarian cancer cases and controls.

# The file *pALL4_ovca_150217_filtered10coverage.rds* is a table where the rows are CpGs and the columns are genomic coordinates, read counts, and methylation percentages. 
# James created this table by subsetting the original table to only those CpGs with 10X coverage.
# The columns containing *55* and *56* respectively refer to the cases and controls.

# ```{r echo=TRUE}
cpgs <- readRDS("/data/SHARE/GINA/pALL4_ovca_150217_filtered10coverage.rds")
# ```

#	1.2. Create a GRanges object from the CpG-level methylation data

# As well as the percentage methylation scores, I've included a metadata column called 'deltaMeth' containing the difference in methylation between cases and controls.

# ```{r echo=TRUE}
cpgs.gr <- GRanges(cpgs[,"chr.x"],IRanges(cpgs[,"start.x"],cpgs[,"end.x"]))
mcols(cpgs.gr) <- as.data.frame(cbind(methCase=cpgs[,"meth55"],methControl=cpgs[,"meth56"],deltaMeth=cpgs[,"meth55"]-cpgs[,"meth56"]))
# ```

# 2. Get a list of all Ensembl genes in hg19 with the genomic co-ordinates of the TSS and TES.
# 	2.1. Use biomaRt to retrieve the Ensembl ID, TSS, TES, exon start/end coordinates, and exon rank of all genes in hg19.

# 'host="grch37.ensembl.org"' specifies the hg19 version of Ensembl.

# ```{r echo=TRUE}
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
# ```

# 	2.2. Remove entries on non-canonical chromosomes and duplicate Ensembl IDs.

# I only want canonical chromosomes so I'll create a vector containing their names to filter the *getBM()* output. 
# The 'chromosome_name' attribute does not include the string 'chr', so I'll remove this from my vector of canonical chromosomes.
# I only want transcribed genes, so I'll add "with_hgnc","with_refseq_mrna","with_refseq_ncrna" to my list of filters.
# 

# ```{r echo=TRUE}
chromosomes <- levels(seqnames(cpgs.gr))
chromosomes <- gsub("chr","",chromosomes)
filters <- c("chromosome_name","with_hgnc","with_refseq_mrna","with_refseq_ncrna")
values <- list(chromosomes,TRUE,TRUE,FALSE)
genes <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","strand"),filter=filters,value=values,mart=mart,uniqueRows=T)
# ```

# These filters return 18910 Ensembl IDs.

# ```{r echo=TRUE}
nrow(genes)
# ```

# The CpG methylation data uses "chr" in the chromosome names. I'll paste "chr" at the start of the chromosome names from biomaRt.

# ```{r echo=TRUE}
genes[,"chromosome_name"] <- paste0("chr",genes[,"chromosome_name"])
# ```

# 3. Create a matrix of methylation values centred around the TSS of each gene.

# To take into account direction of transcription, I'll use the "end_position" as the TSS for those where strand=-1

# ```{r echo=TRUE}
start_position <- genes[,"start_position"]
antisense.index <- which(genes[,"strand"]<0)
start_position[antisense.index] <- genes[antisense.index,"end_position"]
# ```

# I'll create a GRanges object centered on the TSS of each Ensembl gene.

# ```{r echo=TRUE}
tss.gr <- GRanges(genes[,"chromosome_name"],IRanges(start_position,start_position))
mcols(tss.gr) <- genes[,"ensembl_gene_id"]
# ```

# Set the parameters of the extension around the TSS and the size of the bins.

# ```{r echo=TRUE}
extensions <- c(500,1000,2000,3000,5000)
precisions <- c(50,100,200,500)
# ```

# 	3.1. Subset the methylation data to just those CpGs that overlap with the TSS.
# 	3.2. For a given Ensembl ID, get methylation information for all CpGs within its TSS.
# 	3.3. Plot the resulting average methylation scores for each region set.
# There are lots of bins containing NaN. These bins presumably have no CpGs or no measured variation in methylation.
# *heatmap.2()* can't handle NaNs when constructing dendrograms. I will replace NaNs with 0.
# I will use Euclidean distances clustered using the Ward method to create a dendrogram.
# The "ward.D2" method is for use with non-squared distances.
# I will use a blue-white-red colour scheme. Values between -10 and 10 will be coloured white.

# ```{r echo=TRUE}
for(e in 1:length(extensions)){
	
	extend <- extensions[e]
	cat("\n",as.character(Sys.time())," Analysing methylation at TSS +/-",extend," BP...","\n")
	
	tss.extend.gr <- tss.gr+extend
	
	cpgs.tss.gr <- subsetByOverlaps(cpgs.gr,tss.extend.gr)
	
	delta.meths <- lapply(tss.extend.gr,function(x) subsetByOverlaps(cpgs.tss.gr,x))
	
	for(p in 1:length(precisions)){
		
		precision <- precisions[p]
		
		cat("\n",as.character(Sys.time())," Calculating means over ",precision," BP windows...","\n")
	
		delta.meth.means <- array(NA,dim=c(length(delta.meths),2*extend/precision))
	
		progress <- seq(5,100,5)
		previous.progress <- 0
	
		for(n in 1){
	
			cat("\n","Extension =",extend," Precision =",precision,"\n","Task started at",as.character(Sys.time()),
			"\n","Percentage complete",
			"\n","0","                ","100","\n","")

			for(i in 1:length(delta.meths)){
				
				if(i %in% antisense.index){
				
				window.start <- end(tss.extend.gr)[i]
				cpg.position <- window.start-start(delta.meths[[i]])
				
				} else {
								
				window.start <- start(tss.extend.gr)[i]
				cpg.position <- start(delta.meths[[i]])-window.start
				
				}
				
				delta.meth <- mcols(delta.meths[[i]])[[3]]
							
				breaks <- cut(cpg.position,breaks=seq(0,2*extend,precision))
				bins <- split(delta.meth,breaks)
				
				mean.meth <- sapply(bins,mean,na.rm=T)
				
				delta.meth.means[i,] <- mean.meth
				
				pc <- (i/length(delta.meths))*100
				current.progress <- which.min(abs(pc-progress))
				if(current.progress>previous.progress){
					cat("█")
				previous.progress <- current.progress
				}
			
			}
			
		}
		
		delta.meth.means.noNAs <- delta.meth.means
		delta.meth.means.noNAs[is.na(delta.meth.means.noNAs)] <- 0

		variation <- rowSums(delta.meth.means.noNAs)
		length(which(variation==0))
		
		cat("\n",as.character(Sys.time())," Calculating distance matrix and clustering...","\n")
		
		distance.matrix <- dist(delta.meth.means.noNAs, method = "euclidean") 
		clustering <- hclust(1-distance.matrix, method="ward.D2") 

		breaks <- c(seq(range(delta.meth.means,na.rm=T,finite=T)[1],-10,length=4),seq(-9,9,length=2),seq(10,range(delta.meth.means,na.rm=T,finite=T)[2],length=4))
		cols <- colorRampPalette(c("blue", "white", "red"))(n = 9)

		labCol <- rep("",2*(extend/precision))
		labCol[1] <- paste0("-",(extend/1000),"KB")
		labCol[length(labCol)] <- paste0("+",(extend/1000),"KB")

		RowSideColors <- as.character(genes[,"strand"])
		RowSideColors <- gsub("-1","darkgrey",RowSideColors)
		RowSideColors <- gsub("1","lightgrey",RowSideColors)

		cat("\n",as.character(Sys.time())," Plotting heatmap...","\n")
		
		pngName <- paste0("Heatmap_OvCa_TSS_plus",extend,"_",precision,"BP.png")
		png(pngName,h=12,w=12,res=300,unit="in")
		heatmap.2(delta.meth.means.noNAs,
		trace="n",
		lhei=c(1,4), # layout
		col=cols,breaks=breaks,scale="none", # colour
		symm=F,symkey=F,symbreaks=T,density.info="none",keysize=1,key.title="",key.xlab=expression(paste(Delta," methylation")),key.par=list(cex=1), # key
		dendrogram="row",Colv=F,Rowv=as.dendrogram(clustering), # dendrograms
		labRow="",labCol=labCol,cexCol=1.25,srtCol=0,adjCol = c(0.5,1),xlab="TSS",mar=c(2,2) # row and column labels
		)
		dev.off()
	}
}
# ```