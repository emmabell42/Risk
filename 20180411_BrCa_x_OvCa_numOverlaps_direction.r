# This code compares the distribution of delta20 windows in BrCa
# vs OvCa

setwd("/data/emmabell42/risk/BrCa_x_OvCa")

# Load packages
library(GenomicRanges)
library(venneuler)

# Read in files
rdsFiles <- list.files("..")[grep("20.RDS",list.files(".."))]
for(i in 1:length(rdsFiles)){
toRead <- paste0("../",rdsFiles[i])
rds <- readRDS(toRead)
rds <- GRangesList(rds)
rds <- unlist(rds)
#rds <- reduce(rds) # added to join adjacent windows
toName <- gsub("RDS","gr",rdsFiles[i])
assign(toName,rds)
}

common.brca <- subsetByOverlaps(BrCa_0_delta20.gr,OvCa_0_delta20.gr)
common.brca <- sort(common.brca)
common.ovca <- subsetByOverlaps(OvCa_0_delta20.gr,BrCa_0_delta20.gr)
common.ovca <- sort(common.ovca)

brca.direction <- mcols(common.brca)[,2]-mcols(common.brca)[,1]
ovca.direction <- mcols(common.ovca)[,2]-mcols(common.ovca)[,1]

hypo <- which(brca.direction<0 & ovca.direction<0)
hyper <- which(brca.direction>0 & ovca.direction>0)

index <- sort(c(hypo,hyper))

common.brca <- common.brca[index]
common.ovca <- common.ovca[index]

brca <- length(BrCa_0_delta20.gr)
ovca <- length(OvCa_0_delta20.gr)
overlap <- length(index)

combinations <- c("BrCa"=brca,"OvCa"=ovca,"BrCa&OvCa"=overlap)
v <- venneuler(combinations)
png("Venn_sameDirection.png",w=6,h=6,unit="in",res=300)
plot(v)
dev.off()

# P-value according to a hypergeometric distribution
1-phyper(q=overlap,m=brca,n=7729274-brca,k=ovca)
#[1] 1

# Observed/expected
observed <- length(index)
expected <- (20020633152)/7729274 # integer overflow - did brca*ovca with calculator

oe <- observed/expected
#[1] 0.5821866
