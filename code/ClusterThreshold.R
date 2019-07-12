library(igraph)
sgdcpm_matrix <- cpm[1:6692,2:3458]
sgdcpm_matrix <- as.matrix(sgdcpm_matrix)
sgdcpm_matrix <- log(sgdcpm_matrix+1)
cormat <- cor(t(sgdcpm_matrix))
diag(cormat) <- 0

cutoff <- seq(0,100)/100
max.comp <- c()
edge.den <- c()
no.nodes <- c()
no.edges <- c()
clusters <- c()
for (i in 1:101) {
  cormat1 <- cormat 
  cormat1[ cormat1 < cutoff[i] ] <- 0
  graph <- graph.adjacency(cormat1, mode="lower")
  a <- clusters(graph)
  max.comp[i] <- max(a$csize)
  edge.den[i] <- edge_density(graph)
  no.nodes[i] <- sum(a$csize[a$csize!=1])
  no.edges[i] <- ecount(graph)
  clusters[i] <- a$no
}

cormat1 <- cormat 
cormat1[ cormat1 <= 0 ] <- 0
graph <- graph.adjacency(cormat1, weighted=TRUE, mode="lower")
modules <- decompose.graph(graph)
largest <- which.max(sapply(modules, vcount))
max.comp[i] <- vcount(modules[[largest]])

par(mar=c(5,5,1,5))
plot(summary$cutoff,summary$max.comp,ylab="largest connected component", xaxt = "n", xlab="PCC cutoff",col="red",cex=0.5,pch=16)
par(new=T)
plot(summary$cutoff,summary$edge.den,xaxt = "n", yaxt = "n",ylab = "", xlab = "",col="blue",cex=0.5,pch=16)
axis(side=4)
axis(side=1,at=seq(0,1,by=0.1))
mtext("Network Density", side = 4, line = 2.5)
legend("topright", c("LargestComponent", "NetworkDensity"),col = c("red", "blue"),pch=16)

summary <- read.csv("cutoff_summary.csv")
summary <- summary[,-1]

plot(summary$cutoff,summary$no.nodes,cex=0.5,pch=16,xlab="PCC cutoff", ylab="Number of Nodes", xaxt="n")
axis(side=1,at=seq(0,1,by=0.1))

plot(summary$cutoff,summary$no.edges,cex=0.5,pch=16,xlab="PCC cutoff", ylab="Number of Edges", xaxt="n",yaxt="n")
axis(side=1,at=seq(0,1,by=0.1))
axis(side=2,at=seq(0,20000000,by=2000000))

plot(summary$cutoff,summary$max.comp,cex=0.5,pch=16,xlab="PCC cutoff", ylab="Largest Connected Component", xaxt="n")
axis(side=1,at=seq(0,1,by=0.1))

library(coexnet)
pdf("sgdcutoff.nolog.pdf")
cutoffsgdnolog <- findThreshold(sgdcpm_matrix,method = "correlation",plotting = T)
cutoffsgdnolog
dev.off()

