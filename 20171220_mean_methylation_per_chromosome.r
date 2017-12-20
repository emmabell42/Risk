# This code plots the mean methylation values per chromosome for cases and controls

png("Boxplot_mean_meth_chromosome.png",h=3,w=6,unit="in",res=300)
boxplot(gr.all[[1]]$Casemean,gr.all[[1]]$Controlmean,
gr.all[[2]]$Casemean,gr.all[[2]]$Controlmean,
gr.all[[3]]$Casemean,gr.all[[3]]$Controlmean,
gr.all[[4]]$Casemean,gr.all[[4]]$Controlmean,
gr.all[[5]]$Casemean,gr.all[[5]]$Controlmean,
gr.all[[6]]$Casemean,gr.all[[6]]$Controlmean,
gr.all[[7]]$Casemean,gr.all[[7]]$Controlmean,
gr.all[[8]]$Casemean,gr.all[[8]]$Controlmean,
gr.all[[9]]$Casemean,gr.all[[9]]$Controlmean,
gr.all[[10]]$Casemean,gr.all[[10]]$Controlmean,
gr.all[[11]]$Casemean,gr.all[[11]]$Controlmean,
gr.all[[12]]$Casemean,gr.all[[12]]$Controlmean,
gr.all[[13]]$Casemean,gr.all[[13]]$Controlmean,
gr.all[[14]]$Casemean,gr.all[[14]]$Controlmean,
gr.all[[15]]$Casemean,gr.all[[15]]$Controlmean,
gr.all[[16]]$Casemean,gr.all[[16]]$Controlmean,
gr.all[[17]]$Casemean,gr.all[[17]]$Controlmean,
gr.all[[18]]$Casemean,gr.all[[18]]$Controlmean,
gr.all[[19]]$Casemean,gr.all[[19]]$Controlmean,
gr.all[[20]]$Casemean,gr.all[[20]]$Controlmean,
gr.all[[21]]$Casemean,gr.all[[21]]$Controlmean,
gr.all[[22]]$Casemean,gr.all[[22]]$Controlmean,
gr.all[[23]]$Casemean,gr.all[[23]]$Controlmean,
pch=20,col=c("grey","darkgrey"),xant="",at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,))
dev.off()

Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
chr.names <- paste("Chr",c(1:22,"X"),sep=" ")


png("scatterplot_plot_mean_meth_chromosome.png",h=6,w=6,unit="in",res=300)
par(mfrow=c(5,5),mar=c(2,2,1,1))
for(i in 1:23){
smoothScatter(gr.all[[i]]$Controlmean,gr.all[[i]]$Casemean,colramp = Lab.palette,pch=NA,cex.axis=0.6,main=chr.names[i])
}
dev.off()

gr <- data.frame(c())
for(i in 1:23){
gr.chr <- gr.all[[i]]
Controlmean <- data.frame(chr=i,cohort="Control",means=as.numeric(gr.chr$Controlmean))
Casemean <- data.frame(chr=i,cohort="Cases",means=as.numeric(gr.chr$Casemean))
gr <- rbind(gr,Controlmean,Casemean)
}

ats.lab <- seq(1,by=3,length.out=23)
ats <- rep(ats.lab,each=2)
for(i in seq(1,length(ats),2)){
ats[i] <- ats[i]-0.5
ats[i+1] <- ats[i+1]+0.5
}

png("Boxplot_mean_meth_chromosome.png",h=6,w=12,unit="in",res=300)
boxplot(as.numeric(means)~cohort*chr,data=gr,pch=20,at=ats,xaxt="n",ylab="Methylation percentage (mean/window)",col=c("lightgrey","darkgrey"))
axis(1,ats.lab,chr.names,las=2)
dev.off()

