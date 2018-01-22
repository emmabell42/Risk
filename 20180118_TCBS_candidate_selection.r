# This code compiles the target capture bisulphite seq candidate table
# Candidates must have a high number of CpGs, high difference in methylation and adjacent windows

cancerType <- c("OvCa","BrCa")
min.diff <- 0

gr <- list.files()[grep("RDS",list.files())]
for(i in 1:length(gr)){
	reading <- readRDS(gr[i])
	toName <- paste0(cancerType[i],".diff.gr")
	assign(toName,reading)
}

# Sequence composition analysis

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
	c.count <- rep(NA,length(gr.chr))
	g.count <- rep(NA,length(gr.chr))

	for(k in 1:length(gr.chr)){
		cg.count[k] <- str_count(seqs[k], fixed("CG"))
		c.count[k] <- str_count(seqs[k], fixed("C"))
		g.count[k] <- str_count(seqs[k], fixed("G"))
	}
	
	# Calculate GC%
	GC <- (c.count+g.count)/200
	toName <- paste0(type,".GC")
	assign(toName,GC)
	
	# Calculate normalised CpG fraction
	
	# Saxonov method
	#oe <- (cg.count/(200))/((GC/2)^2)
	
	# Gardiner-Garden and Frommer method
	oe <- (cg.count/(c.count*g.count))*200
	toName <- paste0(type,".CpG.norm")
	assign(toName,oe)
	
	# Which qualify as a CGI?
	
	CGI <- GC>=0.5 & oe >=0.6
	
	mcols(gr.chr) <- cbind(mcols(gr.chr),data.frame(CpGs=cg.count,GC=GC,NormCpG=oe,CGI=CGI))
	
	gr[[j]] <- gr.chr
	}
	
	# Write to R object
	toName <- paste0(cancerType[i],".diff.gr")
	assign(toName,gr)
	
	# Write to RDS
	toName <- paste0(cancerType[i],"_",min.diff,"_delta20_gr.RDS")
	saveRDS(gr,toName)
}
cat("\n")

Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
chr.names <- paste("Chr",c(1:22,"X"),sep=" ")

for(i in 1:length(cancerType)){
type <- cancerType[i]
toName <- paste0(type,"_",min.diff,"_GCvsCpG.png")
GC <- get(paste0(type,".GC"))
CpG <- get(paste0(type,".CpG"))
CGI <- length(which(GC>=0.5 & CpG >=0.6))
png(toName,h=6,w=6,unit="in",res=300)
smoothScatter(GC,CpG,colramp = Lab.palette,pch=NA,cex.axis=0.6)
abline(h=0.6,v=0.5,col="white",lty=2)
legend("topright",legend=CGI,text.col="white",bg=NA,bty="n")
dev.off()
}

# Find adjacent windows

for(i in 1:length(cancerType)){

	type <- cancerType[i]
	gr <- get(paste0(cancerType[i],".diff.gr"))
		
	cat("\n","Analysing cancer type: ",type)
	
	for(j in 1:23){
		gr.chr <- gr[[j]]
		gr.reduce <- reduce(gr.chr)
		gr.reduce <- gr.reduce[which((end(gr.reduce)-start(gr.reduce))>200)]
		
		overlaps <- countOverlaps(gr.chr,gr.reduce)
		
		mcols(gr.chr) <- cbind(mcols(gr.chr),data.frame(adjacent=overlaps))
		gr[[j]] <- gr.chr
	}

	# Write to R object
	toName <- paste0(cancerType[i],".diff.gr")
	assign(toName,gr)
	
	# Write to RDS
	toName <- paste0(cancerType[i],"_",min.diff,"_delta20_gr.RDS")
	saveRDS(gr,toName)	
	
}

chr.names <- c(1:22,"X")

for(i in 1:length(cancerType)){

	type <- cancerType[i]
	gr <- get(paste0(cancerType[i],".diff.gr"))
	adj <- c()
	cgi <- c()
	both <- c()
	gr.pyro.candidates <- gr
	
	for(j in 1:23){
		gr.chr <- gr[[j]]
		adj <- c(adj,table(mcols(gr.chr)[[11]])[[2]])
		cgi <- c(cgi,table(mcols(gr.chr)[[10]])[[2]])
		both <- c(both,length(which(mcols(gr.chr)[[10]]==1 & mcols(gr.chr)[[11]]==1)))
		gr.pyro.candidates[[j]] <- gr.chr[which(mcols(gr.chr)[[10]]==1 & mcols(gr.chr)[[11]]==1)]
	}
	
	toName <- paste0(type,"_",min.diff,"_numAdjacentWindows.png")
	png(toName,h=3,w=6,unit="in",res=300)
	barplot(adj,names=chr.names,ylab="# delta 20 windows",cex.names=0.5)
	dev.off()

	toName <- paste0(type,"_",min.diff,"_pcAdjacentWindows.png")
	png(toName,h=3,w=6,unit="in",res=300)
	barplot((adj/sum(adj))*100,names=chr.names,ylab="PC of delta 20 windows",cex.names=0.5)
	dev.off()
	
	toName <- paste0(type,"_",min.diff,"_numCGI.png")
	png(toName,h=3,w=6,unit="in",res=300)
	barplot(cgi,names=chr.names,ylab="# delta 20 windows",cex.names=0.5)
	dev.off()
	
	toName <- paste0(type,"_",min.diff,"_numAdjCGI.png")
	png(toName,h=3,w=6,unit="in",res=300)
	barplot(both,names=chr.names,ylab="# delta 20 windows",cex.names=0.5)
	legend("top",legend=sum(both),bg=NA,bty="n")
	dev.off()
}

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
	
	#df.diff.towrite <- cbind(df.diff[,1:3],id,"","0")
	#toName <- paste0(cancerType[i],"_",min.diff,"_delta20_CI.txt")
	#write.table(df.diff.towrite,toName,sep="\t",quote=F,row.names=F,col.names=F)

}

colnames(homer_BrCa_0_delta20_CI.txt)[1] <- "id"
colnames(homer_OvCa_0_delta20_CI.txt)[1] <- "id"

BrCa.diff.df <- merge(BrCa.diff.df,homer_BrCa_0_delta20_CI.txt,by="id")
OvCa.diff.df <- merge(OvCa.diff.df,homer_OvCa_0_delta20_CI.txt,by="id")
BrCa.diff.df <- BrCa.diff.df[which(BrCa.diff.df$CGI==TRUE & BrCa.diff.df$adjacent==1),]
OvCa.diff.df <- OvCa.diff.df[which(OvCa.diff.df$CGI==TRUE & OvCa.diff.df$adjacent==1),]
write.table(BrCa.diff.df,"BrCa_0_pyroCandidates.txt",sep="\t",row.names=F,quote=F)
write.table(OvCa.diff.df,"OvCa_0_pyroCandidates.txt",sep="\t",row.names=F,quote=F)

gwas <- read.table("gwas-downloaded_2018-01-12-ovarian_cancer.tsv",sep="\t",head=T,fill=T,stringsAsFactors=F)

gwas.genes <- c(gwas$MAPPED_GENE,gwas$REPORTED.GENE.S.)
gwas.genes <- unlist(str_split(gwas.genes,","))
gwas.genes <- unlist(str_split(gwas.genes,"-"))
gwas.genes <- unlist(str_split(gwas.genes," "))
gwas.genes <- unique(gwas.genes)
gwas.genes[which(gwas.genes=="")] <- NA
gwas.genes <- as.character(na.omit(gwas.genes))

which(OvCa.diff.df$Gene.name %in% gwas.genes)