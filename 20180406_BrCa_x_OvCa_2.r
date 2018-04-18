# This code compares the distribution of delta20 windows in BrCa
# vs OvCa

setwd("/data/emmabell42/risk/BrCa_x_OvCa")

# Load packages
library(GenomicRanges)
library(vioplot)

# Read in files
rdsFiles <- list.files("..")[grep("20.RDS",list.files(".."))]
for(i in 1:length(rdsFiles)){
	toRead <- paste0("../",rdsFiles[i])
	rds <- readRDS(toRead)
	rds <- GRangesList(rds)
	rds <- unlist(rds)
	toName <- gsub("RDS","gr",rdsFiles[i])
	assign(toName,rds)
}

cancerTypes <- c("BrCa","OvCa")
patientTypes <- c("cases","controls")
regionTypes <- c("common","unique")
deltaMethDirections <- c("hypo","hyper")

common.BrCa <- subsetByOverlaps(BrCa_0_delta20.gr,OvCa_0_delta20.gr)
common.BrCa <- sort(common.BrCa)
ids.BrCa <- paste(seqnames(common.BrCa),start(common.BrCa),sep="_")
common.OvCa <- subsetByOverlaps(OvCa_0_delta20.gr,BrCa_0_delta20.gr)
common.OvCa <- sort(common.OvCa)
ids.OvCa <- paste(seqnames(common.OvCa),start(common.OvCa),sep="_")
identical(ids.BrCa,ids.OvCa)

BrCa.direction <- mcols(common.BrCa)[,2]-mcols(common.BrCa)[,1]
OvCa.direction <- mcols(common.OvCa)[,2]-mcols(common.OvCa)[,1]

hypo <- which(BrCa.direction<0 & OvCa.direction<0)
hyper <- which(BrCa.direction>0 & OvCa.direction>0)
index <- sort(c(hypo,hyper))

common.BrCa <- common.BrCa[index]
common.OvCa <- common.OvCa[index]

for(i in 1:length(cancerTypes)){

	cancer <- cancerTypes[i]
	cancer.gr <- paste0(cancer,"_0_delta20.gr")
	cancer.gr <- get(cancer.gr)
	
	common <- paste0("common.",cancer)
	common <- get(common)
	
	unique. <- cancer.gr[!(cancer.gr %over% common)]
	
	toName <- paste0("unique.",cancer)
	assign(toName,unique.)
	
	common.deltaMeth <- mcols(common)[,2]-mcols(common)[,1]
	unique.deltaMeth <- mcols(unique.)[,2]-mcols(unique.)[,1]
	
	toName <- paste0("unique.",cancer,".deltaMeth")
	assign(toName,unique.deltaMeth)
	
	toName <- paste0("common.",cancer,".deltaMeth")
	assign(toName,common.deltaMeth)
	
	pngName <- paste0("Vioplot_deltaMeth_",cancer,"_commonVsUnique.png")
	png(pngName,w=6,h=6,unit="in",res=300)
	vioplot(common.deltaMeth,unique.deltaMeth,col="grey",names=c("Common","Unique"))
	title(main=cancer,xlab="Region type",ylab=expression(paste(Delta, "methylation (cases-controls)")))
	dev.off()
		
	}

png("Vioplot_deltaMeth_commonVsUnique.png",w=6,h=6,unit="in",res=300)
vioplot(c(common.BrCa.deltaMeth,common.OvCa.deltaMeth),c(unique.BrCa.deltaMeth,unique.OvCa.deltaMeth),col="grey",names=c("Common","Unique"))
title(main="BrCa + OvCa",xlab="Region type",ylab=expression(paste(Delta, "methylation (cases-controls)")))
dev.off()
t.test(c(common.BrCa.deltaMeth,common.OvCa.deltaMeth),c(unique.BrCa.deltaMeth,unique.OvCa.deltaMeth))
#        Welch Two Sample t-test
#
#data:  c(common.BrCa.deltaMeth, common.OvCa.deltaMeth) and #c(unique.BrCa.deltaMeth, unique.OvCa.deltaMeth)
#t = -18.706, df = 3099.4, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -9.240920 -7.487513
#sample estimates:
#mean of x mean of y
#-19.94807 -11.58385

cor.test(common.BrCa.deltaMeth,common.OvCa.deltaMeth)

#        Pearson's product-moment correlation
#
#data:  common.BrCa.deltaMeth and common.OvCa.deltaMeth
#t = 77.793, df = 1506, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8843036 0.9044621
#sample estimates:
#      cor
#0.8948381


common.BrCa.cases <- mcols(common.BrCa)[,2]
common.BrCa.controls <- mcols(common.BrCa)[,1]
common.OvCa.cases <- mcols(common.OvCa)[,2]
common.OvCa.controls <- mcols(common.OvCa)[,1]
unique.BrCa.cases <- mcols(unique.BrCa)[,2]
unique.BrCa.controls <- mcols(unique.BrCa)[,1]
unique.OvCa.cases <- mcols(unique.OvCa)[,2]
unique.OvCa.controls <- mcols(unique.OvCa)[,1]

for(i in 1:length(cancerTypes))
{
	cancer <- cancerTypes[i]
	for(j in 1:length(regionTypes))
	{
		region <- regionTypes[j]
		regions <- paste(region,cancer,sep=".")
		regions <- get(regions)
		for(k in 1:length(patientTypes))
		{
			patients <- patientTypes[k]
			if(patients=="cases"){
				meth <- mcols(regions)[,2]
				ci <- mcols(regions)[,5]-mcols(regions)[,6]
			}
			else {
				meth <- mcols(regions)[,1]
				ci <- mcols(regions)[,3]-mcols(regions)[,4]
			}
			methName <- paste(region,cancer,patients,sep=".")
			assign(rName,meth)
			ciName <- paste(region,cancer,patients,"ci",sep=".")
			assign(ciName,ci)
		}
	}
}

png("Vioplot_absMeth_common.png",w=6,h=6,unit="in",res=300)
vioplot(common.BrCa.cases,common.BrCa.controls,common.OvCa.cases,common.OvCa.controls,col="grey",names=rep(c("Cases","Controls"),2))
title("Common",xlab="Cancer and patient type",ylab="Absolute methylation")
axis(1, at = c(1.5, 3.5), labels = c("BrCa","OvCa"), padj = 1.5,tick=F) 
dev.off()

png("Vioplot_absMeth_unique.png",w=6,h=6,unit="in",res=300)
vioplot(unique.BrCa.cases,unique.BrCa.controls,unique.OvCa.cases,unique.OvCa.controls,col="grey",names=rep(c("Cases","Controls"),2))
title("Unique",xlab="Cancer and patient type",ylab="Absolute methylation")
axis(1, at = c(1.5, 3.5), labels = c("BrCa","OvCa"), padj = 1.5,tick=F) 
dev.off()

png("Vioplot_absMeth_BrCa.png",w=6,h=6,unit="in",res=300)
vioplot(common.BrCa.cases,common.BrCa.controls,unique.BrCa.cases,unique.BrCa.controls,col="grey",names=rep(c("Cases","Controls"),2))
title("BrCa",xlab="Region and patient type",ylab="Absolute methylation")
axis(1, at = c(1.5, 3.5), labels = c("Common","Unique"), padj = 1.5,tick=F) 
dev.off()

png("Vioplot_absMeth_OvCa.png",w=6,h=6,unit="in",res=300)
vioplot(common.OvCa.cases,common.OvCa.controls,unique.OvCa.cases,unique.OvCa.controls,col="grey",names=rep(c("Cases","Controls"),2))
title("OvCa",xlab="Region and patient type",ylab="Absolute methylation")
axis(1, at = c(1.5, 3.5), labels = c("Common","Unique"), padj = 1.5,tick=F) 
dev.off()

png("Vioplot_absMeth_OvCa_2.png",w=6,h=6,unit="in",res=300)
vioplot(common.OvCa.cases,unique.OvCa.cases,common.OvCa.controls,unique.OvCa.controls,col="grey",names=rep(c("Common","Unique"),2))
title("OvCa",xlab="Region and patient type",ylab="Absolute methylation")
axis(1, at = c(1.5, 3.5), labels = c("Cases","Controls"), padj = 1.5,tick=F) 
dev.off()

for(i in 1:length(cancerTypes))
{
	cancer <- cancerTypes[i]
	pngName <- paste0("DensityPlot_absMeth_",cancer,".png")
	png(pngName,w=6,h=6,unit="in",res=300)
	par(mfrow=c(2,1))
	
	for(j in 1:length(patientTypes))
	{
		patients <- patientTypes[j]
		
		common <- paste("common",cancer,patients,sep=".")
		unique. <- paste("unique",cancer,patients,sep=".")
		common <- get(common)
		unique. <- get(unique.)
		
		p <- t.test(common,unique.)[[3]]
		p <- format(p,digits=3)
		medians <- c(median(common),median(unique.))
		
		common.density <- density(common,bw=2)
		unique.density <- density(unique.,bw=2)
		
		densities <- c(common.density$y,unique.density$y)
		ylim <- c(0,max(densities)*1.1)
		main <- paste(cancer,patients,sep=" ")
		
		plot(unique.density,type="l",main=main,xlab=paste0("P = ",p),ylim=ylim)
		points(common.density,type="l",col="red")
		abline(v=medians,col=c("red","black"),lty=2)
		
		if(j==1)
		{
			legend("topleft",legend=c("Unique","Common"),col=c("black","red"),lwd="1",bty="n")
		}		
		
		}
	
	dev.off()
	
	}
	
# Variation

png("Vioplot_methVariation_BrCa.png",w=6,h=6,unit="in",res=300)
vioplot(common.BrCa.cases.ci,unique.BrCa.cases.ci,common.BrCa.controls.ci,unique.BrCa.controls.ci,col="grey",names=rep(c("Common","Unique"),2))
title("BrCa",xlab="Region and patient type",ylab="Upper CI-Lower CI")
axis(1, at = c(1.5, 3.5), labels = c("Cases","Controls"), padj = 1.5,tick=F) 
dev.off()
png("Vioplot_methVariation_OvCa.png",w=6,h=6,unit="in",res=300)
vioplot(common.OvCa.cases.ci,unique.OvCa.cases.ci,common.OvCa.controls.ci,unique.OvCa.controls.ci,col="grey",names=rep(c("Common","Unique"),2))
title("OvCa",xlab="Region and patient type",ylab="Upper CI-Lower CI")
axis(1, at = c(1.5, 3.5), labels = c("Cases","Controls"), padj = 1.5,tick=F) 
dev.off()
