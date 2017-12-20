# This is the code James used to read the breast cancer WGBS into R and look for windows where the confidence intervals between cases and controls don't overlap
# The RDS files for each chromosome are stored here: /data/SHARE/GINA/methpredict/

wg.wgbs.windows<-readRDS("/data/SHARE/GINA/methpredict/WGBS_BrCa_200bpWindowCIs_chr1.rds")
sum(wg.wgbs.windows$ControlLCI-wg.wgbs.windows$CaseUCI>0,na.rm=T)
sum(wg.wgbs.windows$CaseLCI-wg.wgbs.windows$ControlUCI>0,na.rm=T)

path
data/SHARE/GINA/methpredict/

for (i in c(1:22,"X"){
	filez<-paste("WGBS_BrCa_200bpWindowCIs_chr",i,".rds", sep="")
	wg.wgbs.windows<-readRDS(filez)
	wh1<-which(wg.wgbs.windows$ControlLCI-wg.wgbs.windows$CaseUCI>5)
	wh2<-which(wg.wgbs.windows$CaseLCI-wg.wgbs.windows$ControlUCI>5)
	gr.all.diff<-wg.wgbs.windows[c(wh1,wh2),]
	gr.all.5<-c(gr.all.5,gr.all.diff)
	gr.all<-c(gr.all,wg.wgbs.windows)
}

sum(wg.wgbs.windows$ControlLCI-wg.wgbs.windows$CaseUCI>5,na.rm=T)
sum(wg.wgbs.windows$CaseLCI-wg.wgbs.windows$ControlUCI>5,na.rm=T)

gr.all<-wg.wgbs.windows

###select which dataset to use
gr.all<-wg.wgbs.windows
wh1<-which(wg.wgbs.windows$ControlLCI-wg.wgbs.windows$CaseUCI>0)
gr.all.20<-wg.wgbs.windows[wh1,]

###select which dataset to use
gr.all<-wg.wgbs.windows
wh2<-which(wg.wgbs.windows$CaseLCI-wg.wgbs.windows$ControlUCI>0)
gr.all.20<-wg.wgbs.windows[wh2,]


###select which dataset to use
gr.all<-wg.wgbs.windows
gr.all.20<-wg.wgbs.windows[c(wh1,wh2),]










####################
x1<-findOverlaps(gr.cgi, gr.all) 

x2<-findOverlaps(gr.tss, gr.all) 

x3<-findOverlaps(gr.450k, gr.all) 

x4<-findOverlaps(gr.epic, gr.all) 

exp1<-c(length(x1), length(x2),length(x3),length(x4))

exp2<-c(length(x1), length(x2),length(x3),length(x4), length(gr.all)-sum(exp1))


y1<-findOverlaps(gr.cgi, gr.all.5) 

y2<-findOverlaps(gr.tss, gr.all.5) 

y3<-findOverlaps(gr.450k, gr.all.5) 

y4<-findOverlaps(gr.epic, gr.all.5) 

obs1<-c(length(y1), length(y2),length(y3),length(y4))
obs2<-c(length(y1), length(y2),length(y3),length(y4), length(gr.all.5)-sum(obs1))

oe<-(obs2/length(gr.all.5)*100)/ (exp2/length(gr.all)*100)

chi.table<-array(NA,c(5,3))
rownames(chi.table)<-c("CGI","Genebody","450k","EPIC","Intergenic")
colnames(chi.table)<-c("est","lower","upper")

library(epitools)

for(i in 1:5)
{
chi<-c(length(gr.all),length(gr.all.5),exp2[i],obs2[i])
chi.table[i,]<-oddsratio(chi/10)$measure[2,]
}


chi.table.down<-as.data.frame(chi.table)

chi.table.up<-as.data.frame(chi.table)




library(ggplot2)

ggplot(cbind(chi.table.up[,1]), geom_bar(position=position_dodge(), geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))


