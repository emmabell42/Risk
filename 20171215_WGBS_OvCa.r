# This code reads into R the WGBS CI data
# It then looks for windows where the confidence intervals between cases and controls don't overlap by X amount
# This code is adapted from 20171215_WGBS_BrCa_JF.r
# The RDS files for each chromosome are stored here: /data/SHARE/GINA/methpredict/

# Define the cancer type you want to read in as "OvCa", "BrCa", or c("OvCa","BrCa") if you want both data sets
cancerType <- "BrCa"

# Set a minimum difference for the control and case CIs
min.diff <- 0
gr.all <- c()

# Read in the rds for each chromosome and concatenate into a single list of GRanges objects
for(i in 1:length(cancerType)){

	type <- cancerType[i]

	for (j in c(1:22,"X")){
		filez<-paste("//data/SHARE/GINA/methpredict/WGBS_",type,"_200bpWindowCIs_chr",j,".rds", sep="")
		cat("\n","Reading ",filez)
		
		system.time({
			wg.wgbs.windows<-readRDS(filez)
		})
		
		wh1<-which(wg.wgbs.windows$ControlLCI-wg.wgbs.windows$CaseUCI>min.diff)
		wh2<-which(wg.wgbs.windows$CaseLCI-wg.wgbs.windows$ControlUCI>min.diff)
		gr.all.diff<-wg.wgbs.windows[c(wh1,wh2),]
		gr.all<-c(gr.all,gr.all.diff)

}

toName <- paste0("gr.",type)
names(gr.all) <- paste("chr",c(1:22,"X"),sep="")
assign(toName,gr.all)

}

rm(i)
rm(filez)
rm(wg.wgbs.windows)
rm(wh1)
rm(wh2)
rm(gr.all)
rm(type)
rm(gr.all.diff)
rm(toName)
rm(j)
