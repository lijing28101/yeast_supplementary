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

