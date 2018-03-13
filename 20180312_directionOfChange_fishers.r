# This code performs Fisher's tests on the direction of 
# delta20 windows

library(GenomicRanges)
library(stringr)

rdsFiles <- list.files()[grep("allWindows.RDS",list.files())]
rdsFiles <- c(rdsFiles,list.files()[grep("_delta20.RDS",list.files())])

for(i in 1:length(rdsFiles)){
rds <- readRDS(rdsFiles[i])
rds <- GRangesList(rds)
rds <- unlist(rds)
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

cancerType <- c("BrCa","OvCa")

# Equal proportion testing

props <- array(NA,dim=c(4,2))
colnames(props) <- cancerType
rownames(props) <- c("est","lci","uci","p")

for(i in 1:length(cancerType)){

	cancer <- cancerType[i]
	delta20 <- paste0(cancer,"_0_delta20.gr")
	delta20.gr <- get(delta20)
	
	means <- mcols(delta20.gr)[,1:2]
	hypo <- length(which(means[,2]<means[,1]))
	
	prop <- prop.test(hypo,nrow(means),p=0.5)
	props[1,i] <- prop[[4]]
	props[2:3,i] <- prop[[6]]
	props[4,i] <- prop[[3]]
	
}

# 

directions <- c("hypo","hyper")

chi.table <- array(NA,dim=c(4,2))
colnames(props) <- cancerType
rownames(props) <- c("est","lci","uci","p")

for(i in 1:length(cancerType)){

	cat("\n","Began working on",cancerType[i],"at",as.character(Sys.time()),
	"\n","Percentage complete",
	"\n","0","               ","100","\n")

	cancer <- cancerType[i]
	delta20 <- paste0(cancer,"_0_delta20.gr")
	delta20.gr <- get(delta20)
	allWindows <- paste0(cancer,"_allWindows.gr")
	allWindows.gr <- get(allWindows)
	exp1 <- c()
	obs1 <- c()
	
	means.delta20 <- mcols(delta20.gr)[,1:2]
	means.allWindows <- mcols(allWindows.gr)[,1:2]

	for(j in 1:length(directions)){
	
		direction <- directions[j]
		cat("Direction of methylation change:",direction,"\n")
		previous.progress <- 0
		cat("Calculating the observed and expected distributions...","\n")
		
		if(direction=="hypo"){
	
		index <- which(means.allWindows[,2]<means.allWindows[,1])
				
		}
		
		else{
		
		index <- which(means.allWindows[,2]>means.allWindows[,1])
		
		}
		
		direction.gr <- allWindows.gr[index]
		
		x <- findOverlaps(direction.gr, allWindows.gr)
		y <- findOverlaps(direction.gr, delta20.gr)
				
		exp1 <- c(exp1,length(x))
		obs1<-c(obs1, length(y))
		
		exp2 <- c(exp1, length(allWindows.gr)-sum(exp1))
		obs2 <- c(obs1, length(direction.gr)-sum(obs1))

		oe <- (obs2/length(direction.gr)*100)/ (exp2/length(allWindows.gr)*100)

		chi.table <- array(NA,c(length(directions),4))
		rownames(chi.table) <- directions
		colnames(chi.table) <- c("est","lower","upper","p-val")
	
		for(k in 1:length(directions)){
		
		cat("\n","Performing Chi Square tests...","\n")
	
		chi <- matrix(c(obs2[k],exp2[k],length(direction.gr),length(allWindows.gr)),ncol=2,nrow=2)
		chi.table[k,1] <- fisher.test(chi)[[3]]
		chi.table[k,2:3] <- fisher.test(chi)[[2]]
		chi.table[k,4] <- fisher.test(chi)[[1]]
		
		}
	
	chi.table[,4] <- p.adjust(chi.table[,4],method="BH")
	
	chiName <- paste0(cancer,".chiTable")
	assign(chiName,chi.table)
	
}
