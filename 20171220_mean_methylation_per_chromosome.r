# This code plots the mean methylation values per chromosome for cases and controls
# It requires a vectors called 'cancerType' containing either "OvCa", "BrCa", or c("OvCa","BrCa") depending on which data are being used
# It also requires a vector called 'min.diff' that specifies the minimum overlap in CI that was used in reading in the RDS

Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
chr.names <- paste("Chr",c(1:22,"X"),sep=" ")

ats.lab <- seq(1,by=3,length.out=23)
ats <- rep(ats.lab,each=2)
for(i in seq(1,length(ats),2)){
ats[i] <- ats[i]-0.5
ats[i+1] <- ats[i+1]+0.5
}

for(i in 1:length(cancerType))
	{
	plotName <- paste0(cancerType[i],"_",min.diff,"_densityScatterplot_meanMeth_byChr.png")
	gr.all <- get(paste0("gr.",cancerType[i]))

	# Plot density scatterplots
	
	png(plotName,h=6,w=6,unit="in",res=300)
	par(mfrow=c(5,5),mar=c(2,2,1,1))
		for(j in 1:23){
		smoothScatter(gr.all[[j]]$Controlmean,gr.all[[j]]$Casemean,colramp = Lab.palette,pch=NA,cex.axis=0.6,main=chr.names[j])
		}
	dev.off()

	gr <- data.frame(c())
		for(j in 1:23){
		gr.chr <- gr.all[[j]]
		Controlmean <- data.frame(chr=j,cohort="Control",means=as.numeric(gr.chr$Controlmean))
		Casemean <- data.frame(chr=j,cohort="Cases",means=as.numeric(gr.chr$Casemean))
		gr <- rbind(gr,Controlmean,Casemean)
		}
	
	plotName <- paste0(cancerType[i],"_",min.diff,"_boxplot_meanMeth_byChr.png")
	
	# Plot boxplots by chromosome
	
	png(plotName,h=6,w=12,unit="in",res=300)
	boxplot(as.numeric(means)~cohort*chr,data=gr,pch=20,at=ats,xaxt="n",ylab="Methylation percentage (mean/window)",col=c("lightgrey","darkgrey"))
	axis(1,ats.lab,chr.names,las=2)
	dev.off()
	
	plotName <- paste0(cancerType[i],"_",min.diff,"_boxplot_meanMeth.png")
	
	# Plot boxplots overall
	
	png(plotName,h=6,w=3,unit="in",res=300)
	boxplot(as.numeric(means)~cohort,data=gr,pch=20,ylab="Methylation percentage (mean/window)",col=c("lightgrey","darkgrey"),names=c("Controls","Cases"))
	dev.off()
	
}

