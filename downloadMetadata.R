#!/usr/bin/env Rscript

library(SRAdb)
sqlfile <-'SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),'SRAmetadb.sqlite')
dbListTables(sra_con)
#newSRR.txt is the list of SRR run ID
ID <- as.character(read.table("newSRR.txt",header=F)$V1)
c <- sraConvert(ID,sra_con=sra_con)
IDList <- as.character(c$experiment)
rr<-c()
#3457 is the number of runs
for(i in 1:3457)
{
  currentID <- IDList[i]
  rw <- dbGetQuery(sra_con, paste("select * from experiment where ", "experiment_accession like '", currentID, "'",sep=""))
  rr[[i]]<-rw
}
df<-data.frame(matrix(unlist(rr), nrow=length(rr), byrow=TRUE))
write.csv(df,"newInfo.csv")
