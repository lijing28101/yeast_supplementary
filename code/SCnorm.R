#!/usr/bin/env Rscript

library(SCnorm)
counts <- read.table("rawcount.txt",row.names = 1, stringsAsFactors = FALSE)
counts <- as.matrix(counts)
ExampleSimSCData1 <- SingleCellExperiment::SingleCellExperiment(assays = list('counts' = counts))
Conditions <- read.csv("condition.csv",header=F)$V1
pdf("MyNormalizedData_k_evaluation.pdf")
par(mfrow=c(2,2))
DataNorm <- SCnorm(Data=ExampleSimSCData1,Conditions=Conditions,PrintProgressPlots = TRUE,FilterCellNum = 10,NCores=3, useZerosToScale=TRUE)
dev.off()
GenesNotNormalized <- results(DataNorm, type="GenesFilteredOut")
ScaleFactors <- results(DataNorm, type="ScaleFactors")
pdf("check_exampleData_count-depth_evaluation3.pdf", height=5, width=7)
countDeptEst.SCNORM <- plotCountDepth(Data = ExampleSimSCData1, NormalizedData = NormalizedData,
                                      Conditions = Conditions,
                                      FilterCellProportion = .1, NCores=3)
dev.off()
save.image("scnorm.RData")
