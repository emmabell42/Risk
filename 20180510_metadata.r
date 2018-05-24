# This script is part of the Hyperactive Enhancer project. It:
# 	1. Reads in the metadata csv;
#	2. Summarises the metadata;
#	3. Plots the metadata.

setwd("/data/emmabell42/Hyperactive enhancers/MethBase")


# Read in metadata
metadata <- read.csv("MethBase_index.csv")

png("Metadata_source.png",h=6,w=12,unit="in",res=300)
par(mar=c(15,5,5,2))
barplot(sort(table(metadata$Source),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Source of data")
dev.off()

png("Metadata_model.png",h=6,w=6,unit="in",res=300)
par(mar=c(7,5,5,2))
barplot(sort(table(metadata$Source.1),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Model")
dev.off()

png("Metadata_developmentalStage.png",h=6,w=6,unit="in",res=300)
par(mar=c(7,5,5,2))
barplot(sort(table(metadata$Developmental.stage),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Developmental Stage")
dev.off()

png("Metadata_Tissue.png",h=6,w=12,unit="in",res=300)
par(mar=c(12,5,5,2))
barplot(sort(table(metadata$Tissue),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Tissue",cex.names=1)
dev.off()

png("Metadata_platform.png",h=6,w=6,unit="in",res=300)
par(mar=c(7,5,5,2))
barplot(sort(table(metadata$Platform),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Platform")
dev.off()

png("Metadata_age.png",h=6,w=6,unit="in",res=300)
par(mar=c(7,5,5,2))
barplot(sort(table(metadata$Age),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Age")
dev.off()

png("Metadata_parity.png",h=6,w=6,unit="in",res=300)
par(mar=c(7,5,5,2))
barplot(sort(table(metadata$Parity),de=T),las=2,col="lightgrey",ylab="Number of experiments",main="Parity")
dev.off()
