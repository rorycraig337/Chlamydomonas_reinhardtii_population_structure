setwd("/Users/rory/Documents/research/projects/chlamy_biogeography/github/population_structure/isolation_by_distance/")
library(vegan)

#plot A and B
IBDaP <- read.csv("NA1.4D.R.csv", header=FALSE)
IBDbP <- read.csv("NA2.4D.R.csv", header=FALSE)
plot(IBDaP$V4,IBDaP$V3, col="black",bg="#2E509E",pch=21,cex=1.25, xlab = "distance (km)", ylab = "genetic distance", ylim = c(0.01,0.032), xlim = c(0,2200))
abline(lm(IBDbP$V3~IBDbP$V4), col="#CA1D17", lwd=1.25)
points(IBDbP$V4,IBDbP$V3, col="black",bg="#CA1D17", pch=22,cex=1.25)
clip(-100,600, -100, 100)
abline(lm(IBDaP$V3~IBDaP$V4), col="#2E509E", lwd=1.25)

IBDa <- read.csv("A_4D.all_sites.nuclear.matrix.R2.csv", header=FALSE)
IBDb <- read.csv("B_4D.all_sites.nuclear.matrix.R2.csv", header=FALSE)

#test for only A
A <- with(IBDa, sort(unique(c(as.character(V1), as.character(V2)))))
A
Ma <- array(0, c(length(A), length(A)), list(A, A))
ia <- match(IBDa$V1, A)
ja <- match(IBDa$V2, A)
Ma[cbind(ia,ja)] <- Ma[cbind(ja,ia)] <- IBDa$V3
Ma
Qa <- array(0, c(length(A), length(A)), list(A, A))
Qa[cbind(ia,ja)] <- Qa[cbind(ja,ia)] <- IBDa$V4
Qa
mantel(as.dist(Ma), as.dist(Qa), method="spearman", permutations=999)

#test for only B
B <- with(IBDb, sort(unique(c(as.character(V1), as.character(V2)))))
B
Mb <- array(0, c(length(B), length(B)), list(B, B))
ib <- match(IBDb$V1, B)
jb <- match(IBDb$V2, B)
Mb[cbind(ib,jb)] <- Mb[cbind(jb,ib)] <- IBDb$V3
Mb
Qb <- array(0, c(length(B), length(B)), list(B, B))
Qb[cbind(ib,jb)] <- Qb[cbind(jb,ib)] <- IBDb$V4
Qb
mantel(as.dist(Mb), as.dist(Qb), method="spearman", permutations=999)

#Mantel tests
IBD <- read.csv("4D.all_sites.nuclear.matrix.R2.csv", header=FALSE)
key <- read.table("clade_key.txt")
key
an <- with(IBD, sort(unique(c(as.character(V1), as.character(V2)))))
an
M <- array(0, c(length(an), length(an)), list(an, an))
i <- match(IBD$V1, an)
j <- match(IBD$V2, an)
M[cbind(i,j)] <- M[cbind(j,i)] <- IBD$V3
M
Q <- array(0, c(length(an), length(an)), list(an, an))
Q[cbind(i,j)] <- Q[cbind(j,i)] <- IBD$V4
Q
R <- array(0, c(length(an), length(an)), list(an, an))
R[cbind(i,j)] <- R[cbind(j,i)] <- IBD$V5
R

mantel(as.dist(M), as.dist(Q), method="pearson", permutations=999)
mantel(as.dist(M), as.dist(Q), method="pearson", permutations=999, strata = key$V2)
mantel.partial(as.dist(M), as.dist(R), as.dist(Q), method="pearson", permutations=999)
