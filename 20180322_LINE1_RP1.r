# LINE-1 investigation
# Research plan 1: Which class of LINE-1 elements appear differentially methylated in OvCa patients compared to controls?
# Subsection 1.1: How many differentially methylated windows overlap each class?

# Load libraries
library(GenomicRanges)

# Read in RDS files
rds <- c("../OvCa_0_delta20.RDS","../OvCa_allWindows.RDS")
r.name <- c("delta20","all")
for(i in 1:length(rds)){
reading <- readRDS(rds[i])
reading <- GRangesList(reading)
reading <- unlist(reading)
assign(r.name[i],reading)
}

rm(rds)

ids <- paste(seqnames(delta20),start(delta20),sep="_")

# Read in LINE-1 class bed files
beds <- list.files()[grep(".bed",list.files())]
gr.name <- gsub("bed","gr",beds)
for(i in 1:length(beds)){
reading <- read.table(beds[i],sep="\t",stringsAsFactors=F)
colnames(reading)[1:3] <- c("chr","start","end")
gr <- GRanges(seqnames=reading$chr,IRanges(reading$start,reading$end))
assign(gr.name[i],gr)
}

rm(beds)
rm(r.name)
rm(gr.name)
rm(reading)
rm(gr)

directions <- c("all","hypo","hyper")

# Find overlaps between the delta20 windows and the three classes of LINE-1 element
l1s <- ls()[grep(".gr",ls())]
means <- mcols(delta20)[,2]-mcols(delta20)[,1]
means.all <- mcols(all)[,2]-mcols(all)[,1]

for(d in 1:length(directions)){

direction <- directions[d]
cat("\n","Analysing ",direction,"\n")

if(direction=="hypo"){
	index <- which(means<0)
}
if(direction=="hyper"){
	index <- which(means>0)
} else {
	index <- 1:length(means)
}

toOverlap <- delta20[index]

overlaps <- array(NA,dim=c(length(toOverlap),3))
colnames(overlaps) <- l1s
rownames(overlaps) <- ids[index]

cat("Counting overlaps...","\n")
for(i in 1:ncol(overlaps)){
	l1.class <- get(l1s[i])
	overlap <- countOverlaps(toOverlap,l1.class)
	overlap[overlap>1] <- 1
	overlaps[,i] <- overlap
}

num.overlaps <- colSums(overlaps)

fileName <- paste0("Barplot_class_num_overlaps_",direction,".png")
png(fileName,h=6,w=3,res=300,unit="in")
barplot(c(num.overlaps[1],num.overlaps[3],num.overlaps[2]),xlab="Class",ylab="# overlapping delta20 windows",names=c("FLI-L1","ORF2-L1","FLnI-L1"),cex.names=0.7,main=direction)
dev.off()

# What proportion of each class of LINE-1 elements contains a delta20 window?
l1.overlaps <- as.list(rep(NA,3))
names(l1.overlaps) <- l1s

cat("Calculating proportions...","\n")
for(i in 1:ncol(overlaps)){
	l1.class <- get(l1s[i])
	overlap <- countOverlaps(l1.class,toOverlap)
	overlap[overlap>1] <- 1
	l1.overlaps[[i]] <- overlap
}

prop.l1.overlaps <- sapply(l1.overlaps, function(x) sum(x)/length(x))

fileName <- paste0("Barplot_class_prop_overlaps_",direction,".png")
png(fileName,h=6,w=3,res=300,unit="in")
barplot(c(prop.l1.overlaps[1],prop.l1.overlaps[3],prop.l1.overlaps[2]),xlab="Class",ylab="Prop. L1 differentially methylated",names=c("FLI-L1","ORF2-L1","FLnI-L1"),cex.names=0.7,ylim=c(0,1),main=direction)
dev.off()

# Subsection 1.3: Are any of the classes of LINE-1 element overrepresented in the delta20 windows?

cat("Calculting odds ratios...","\n")
chi.table <- array(NA,c(length(l1s),4))
rownames(chi.table) <- l1s
colnames(chi.table) <- c("est","lower","upper","p-val")

if(direction=="hypo"){
	index <- which(means.all<0)
}
if(direction=="hyper"){
	index <- which(means.all>0)
} else {
	index <- 1:length(means.all)
}

all.toOverlap <- all[index]


for(i in 1:length(l1s)){
	l1.class <- get(l1s[i])
	
	all.overlap <- countOverlaps(all.toOverlap,l1.class)
	all.overlap[all.overlap>1] <- 1
	
	signifWAnnot <- num.overlaps[i]
	signifWOAnnot <- length(toOverlap)-signifWAnnot
	allWindowsWAnnot <- sum(all.overlap)
	allWindowsWOAnnot <- length(all.toOverlap)-allWindowsWAnnot

	chi <- matrix(c(signifWAnnot,signifWOAnnot,allWindowsWAnnot,allWindowsWOAnnot),ncol=2,nrow=2)
	chi.table[i,1] <- fisher.test(chi)[[3]]
	chi.table[i,2:3] <- fisher.test(chi)[[2]]
	chi.table[i,4] <- fisher.test(chi)[[1]]

}

chi.table[,4] <- p.adjust(chi.table[,4],method="BH")
chi.table <- chi.table[c(1,3,2),]

fileName <- paste0("Barplot_OR_class_overlap_",direction,".png")
png(fileName,h=6,w=6,res=300,unit="in")
x <- barplot(chi.table[,"est"])
ylim <- c(0,max(chi.table[,"upper"])*1.1)
barplot(chi.table[,"est"],ylim=ylim,cex.names=0.7,ylab="Odds Ratio",xlab="Class",names=c("FLI-L1","ORF2-L1","FLnI-L1"),main=direction)
arrows(x, chi.table[,"upper"], x, chi.table[,"lower"], length=0.05, angle=90, code=3)
abline(h=1,lty=2,lwd=2,col="red")

dev.off()
}
