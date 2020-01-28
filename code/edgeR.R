library(tximport)
library(readr)
library(edgeR)
#import the transcripts list
gene <- read.csv("gene.csv")$gene
#Define the gene name if your transcript name is different from gene name
tx2gene <- data.frame(GENEID=gene,TXNAME=gene)
#import run ID (SRR)
run <- read.csv("comm.csv",header=F)
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,".tsv",sep="")
names(file) <- run$V1
#import the output from kallisto. It include raw counts, efficient length, and tpm 
txi.kallisto.tsv <- tximport(file, type = "kallisto", tx2gene = tx2gene)
df <- txi.kallisto.tsv$counts
#import raw counts into edgeR
y <- DGEList(df)
#may have runs without any expression, just discard them
keep <- y$samples$lib.size != 0
y <- y[ ,keep]
no <-  y$samples$lib.size == 0
colnames(y[,no])
#calculate TMM scaling factor, after that, the matrix is still raw count, it just add a column of scaling factor in the y$sample
y <- calcNormFactors(y)
#calculate cpm, the matrix is after TMM and cpm
cpm <- cpm(y,normalized.lib.sizes=T)
logcpm <- log(cpm+1)

#function for sort cpm from high to low for each gene/ORF
sortbyrow <- function(matrix,nrow,ncol){
  newM <- matrix(nrow = nrow,ncol=ncol)
  for (i in 1:nrow){
    newM[i,] <- sort(as.numeric(matrix[i,]),decreasing = T)
  }
  return(newM)
}

#cpm_new is a cpm matrix sorted by gene/ORF mean expression from high to low 
cpm_order_new <- sortbyrow(cpm_new,36046,3457)

#Data is too large and hard to make a plot with whole data, split to several part.
#braker is the quantile value
png(filename="sgd.png",width=3000,height=1500,type="cairo")
par(mar=c(0,0,0,0))
image(t(cpm_order_new[6692:1,]),col = c("white",rev(heat.colors(4))), breaks=c(0,0.000000001,0.3089,2.712,18.030,1.13e7),axes=F)
box()
dev.off()
par(mar=c(0,0,0,0))
image(t(cpm_order_new[7831:6693,]),col = c("white",rev(heat.colors(4))), breaks=c(0,0.000000001,0.3089,2.712,18.030,1.13e7),axes=F)
box()
par(mar=c(0,0,0,0))
image(t(cpm_order_new[20759:7832,]),col = c("white",rev(heat.colors(4))), breaks=c(0,0.000000001,0.3089,2.712,18.030,1.13e7),axes=F)
box()
png(filename="un.1.png",width=3000,height=1500,type="cairo")
par(mar=c(0,0,0,0))
image(t(cpm_order_new[36046:20760,]),col = c("white",rev(heat.colors(4))), breaks=c(0,0.000000001,0.3089,2.712,18.030,1.13e7),axes=F)
box()
dev.off()
