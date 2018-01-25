# This code identifies and writes out the windows with >=10% and >=20% differences in methylation between cases and controls
# It requires a vectors called 'cancerType' containing either "OvCa", "BrCa", or c("OvCa","BrCa") depending on which data are being used
# It also requires a vector called 'min.diff' that specifies the minimum overlap in CI that was used in reading in the RDS

# Define the minimum differences in methylation that you want to search for
delta <- c(10,20)

for(i in 1:length(cancerType)){
	
	grName <- paste0("gr.",cancerType[i])
	gr <- get(grName)
	
	for(j in 1:length(delta)){
		
		windowsToKeep <- as.list(rep(NA,23))
		deltaToSearch <- delta[j]
		gr.diff <- gr
		
		for(k in 1:length(gr)){
			chr <- gr[[k]]
			diffs <- chr$Casemean-chr$Controlmean
			abs.diffs <- abs(diffs)
			windowsToKeep[[k]] <- which(abs.diffs>=deltaToSearch)
			gr.diff[[k]] <- gr[[k]][windowsToKeep[[k]]]
		}

		# Write out GRanges object as RDS
		toName <- paste0(cancerType[i],"_",min.diff,"_delta",deltaToSearch,".RDS")
		saveRDS(gr.diff,toName)

		# Create a data frame from the GRanges object and write out to run through HOMER for annotation
		
		df.diff <- c() 
		for(k in 1:23){
			df.diff <- rbind(df.diff,data.frame(seqnames=seqnames(gr.diff[[k]]),starts=start(gr.diff[[k]])-1,ends=end(gr.diff[[k]]),mcols(gr.diff[[k]])))
			
		}
		id <- paste(df.diff[,1],df.diff[,2],sep="_")
		df.diff <- cbind(df.diff,id)
		
		toName <- paste0(cancerType[i],".",deltaToSearch,".df")
		assign(toName,df.diff)
		
		df.diff.towrite <- cbind(df.diff[,1:3],id,"","0")
		toName <- paste0(cancerType[i],"_",min.diff,"_delta",deltaToSearch,".txt")
		write.table(df.diff.towrite,toName,sep="\t",quote=F,row.names=F,col.names=F)
	}
}
