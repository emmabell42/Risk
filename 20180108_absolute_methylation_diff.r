# This code plots the absolute differences in methylation values for cases and controls
# It requires a vectors called 'cancerType' containing either "OvCa", "BrCa", or c("OvCa","BrCa") depending on which data are being used
# It also requires a vector called 'min.diff' that specifies the minimum overlap in CI that was used in reading in the RDS

hg19 <- read.table("/data/emmabell42/resources/hg19.chrom.sizes.txt",stringsAsFactors=F)
hg19 <- hg19[c(1:7,9:19,22,20,24,23,8),]

for(i in 1:length(cancerType)){

	bin_0_2 <- c()
	bin_2_5 <- c()
	bin_5_10 <- c()
	bin_10_20 <- c()
	bin_20_n <- c()
	
	gr.all <- get(paste0("gr.",cancerType[i]))

	
	for(j in 1:length(gr.all)){
		chr <- gr.all[[j]]
		diffs <- chr$Casemean-chr$Controlmean
		abs.diffs <- abs(diffs)
		bin_0_2 <- length(which(abs.diffs<2))
		bin_2_5 <- c(bin_2_5,length(which(abs.diffs>=2 & abs.diffs<5)))
		bin_5_10 <- c(bin_5_10,length(which(abs.diffs>=5 & abs.diffs<10)))
		bin_10_20 <- c(bin_10_20,length(which(abs.diffs>=10 & abs.diffs<20)))
		bin_20_n <- c(bin_20_n,length(which(abs.diffs>20)))
	}

	bins <- rbind(bin_0_2,bin_2_5,bin_5_10,bin_10_20,bin_20_n)
	colnames(bins) <- c(1:22,"X")

	toName <- paste0(cancerType[i],"_",min.diff,"_barplot_abs_methDiff_byChr.txt")
	write.table(bins,toName,sep="\t",quote=F)

	n_bins <- colSums(bins)
	bins.prop <- bins

	for(j in 1:23){
		bins.prop[,j] <- bins.prop[,j]/n_bins[j]
	}

	plotName <- paste0(cancerType[i],"_",min.diff,"_barplot_absMethDiff_byChr.png")
		
	png(plotName,h=3,w=6,unit="in",res=300)
	par(mar=c(3,5,2,2))
	barplot(bins.prop,xlab="Chr",ylab="Proportion of windows",cex.axis=0.6,cex.names=0.5)
	dev.off()

	bins.total <- rowSums(bins)
	bins.total.prop <- bins.total/sum(bins.total)

	plotName <- paste0(cancerType[i],"_",min.diff,"_barplot_absMethDiff.png")
	
	png(plotName,h=6,w=3,unit="in",res=300)
	par(mar=c(3,5,2,2))
	barplot(cbind(bins.total.prop),xlab="Chr",ylab="Proportion of windows",cex.axis=0.6,cex.names=0.5)
	dev.off()

	# Write out genomic regions with >20% delta meth

	bins <- as.list(rep(NA,23))
	
	for(j in 1:length(gr.all)){
		chr <- gr.all[[j]]
		diffs <- chr$Casemean-chr$Controlmean
		abs.diffs <- abs(diffs)
		bins[[j]] <- which(abs.diffs>=20)
	}

	gr.diff <- gr.all
	lengths <- c()
	for(j in 1:23){
		gr.diff[[j]] <- gr.diff[[j]][bins[[j]]]
		lengths <- c(lengths,length(bins[[j]]))
	}
	
	# Write out GRanges object as RDS
	toName <- paste0(cancerType[i],"_",min.diff,"_delta20_gr.RDS")
	saveRDS(gr.diff,toName)

	# Calculate the proportion of the genome represented by delta20 windows
	lengths.200 <- lengths*200
	chr.props <- lengths.200/hg19[,3]
	chr.pc <- chr.props*100
	plotName <- paste0(cancerType[i],"_",min.diff,"_barplot_pc_delta20_byChr.png")

	png(plotName,h=3,w=6,unit="in",res=300)
	par(mar=c(3,5,1,1))
	barplot(chr.pc,ylab="Percentage of chr",names=c(1:22,"X"),cex.names=0.5,cex.axis=0.6)
	dev.off()

	# Create a data frame from the GRanges object and write out to run through HOMER for annotation
	
	df.diff <- c() 
	for(j in 1:23){
		df.diff <- rbind(df.diff,data.frame(seqnames=seqnames(gr.diff[[j]]),starts=start(gr.diff[[j]])-1,ends=end(gr.diff[[j]]),mcols(gr.diff[[j]])))
		
	}
	id <- paste(df.diff[,1],df.diff[,2],sep="_")
	df.diff <- cbind(df.diff,id)
	
	toName <- paste0(cancerType[i],".diff.gr")
	assign(toName,df.diff)
	
	df.diff.towrite <- cbind(df.diff[,1:3],id,"","0")
	toName <- paste0(cancerType[i],"_",min.diff,"_delta20_CI.txt")
	write.table(df.diff.towrite,toName,sep="\t",quote=F,row.names=F,col.names=F)
}

