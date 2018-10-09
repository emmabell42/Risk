#for(1){
  #previous.progress <- 0
  #progress <- seq(5,100,5)
  #cat("\n","Began running at ",as.character(Sys.time()),
   # "\n","Percentage complete",
    #"\n","0","               ","100","\n",sep="")
  for(i in 1:length(bedFiles)){
    cat("\n","\t",i,"."," Performing CREAM on bed file ",i," of ",length(bedFiles)," at ",as.character(Sys.time()),"\n",sep="")
    read <- bedFiles[i]
    in_path <- paste0("data/H3K27ac_bed/short_bed/",read)
    cores <- CREAM(in_path=in_path)
    cores.gr <- GRanges(seqnames=cores[,1],ranges=IRanges(as.numeric(cores[,2]),as.numeric(cores[,3])))
    if(i==1){
      cores.gl <- cores.gr
      } else {
      cores.gl <- list(cores.gl,cores.gr)
      }
    #pc <- (i/length(bedFiles))*100
    #current.progress <- which.min(abs(pc-progress))
    #if(current.progress>previous.progress){
    #  cat("█")
    #  previous.progress <- current.progress
  }
  names(cores.gl) <- gsub(".bed",".cores",bedFiles)
#}
