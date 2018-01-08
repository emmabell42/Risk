# This code reads the ovarian cancer WGBS into R
# It then looks for windows where the confidence intervals between cases and controls don't overlap by X amount
# This code is adapted from 20171215_WGBS_BrCa_JF.r
# The RDS files for each chromosome are stored here: /data/SHARE/GINA/methpredict/

# Set a minimum difference for the control and case CIs
min.diff <- 0
gr.all <- c()

# Read in the rds for each chromosome and concatenate into a single list of GRanges objects
for (i in c(1:22,"X")){
	filez<-paste("//data/SHARE/GINA/methpredict/","WGBS_OvCa_200bpWindowCIs_chr",i,".rds", sep="")
	cat("\n","Reading ",filez)
	system.time({
	wg.wgbs.windows<-readRDS(filez)
	})
	wh1<-which(wg.wgbs.windows$ControlLCI-wg.wgbs.windows$CaseUCI>min.diff)
	wh2<-which(wg.wgbs.windows$CaseLCI-wg.wgbs.windows$ControlUCI>min.diff)
	gr.all.diff<-wg.wgbs.windows[c(wh1,wh2),]
	gr.all<-c(gr.all,gr.all.diff)
}
names(gr.all) <- paste("chr",c(1:22,"X"),sep="")
rm(min.diff)
rm(i)
rm(filez)
rm(wg.wgbs.windows)
rm(wh1)
rm(wh2)
rm(gr.all.diff)
