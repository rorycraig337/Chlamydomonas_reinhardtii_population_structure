setwd("/Users/rory/Documents/research/projects/chlamydomonas_reinhardtii_biogeography/github/population_structure/")
library("gdsfmt")
library("SNPRelate")
library("MASS")
library("calibrate")

vcf <- c("4D.snps_20kb.1lab.vcf")
snpgdsVCF2GDS(vcf, "4D_20kb.1lab.gds", method="biallelic.only")
snpgdsSummary("4D_20kb.1lab.gds")
genofile <- snpgdsOpen("4D_20kb.1lab.gds")
pca <- snpgdsPCA(genofile, autosome.only=FALSE)
plot(pca)

pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id, # make data frame
                  EV1 = pca$eigenvect[,1],   # 1st eigenvector
                  EV2 = pca$eigenvect[,2],   # 2nd eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.table(tab, "PCA_4D.txt", sep="\t")

tab <- read.table("PCA_4D.txt", sep="\t")

plot(tab$EV1,tab$EV2,xlab="eigenvector 1 (14.0%)",ylab="eigenvector 2 (7.6%)",
     xlim=c(-0.2,0.35),ylim=c(-0.7,0.2),
     col="white")
points(tab$EV1[1:2],tab$EV2[1:2],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV1[3:8],tab$EV2[3:8],col="black",bg="#CA1D17", pch=22,cex=1.5)
points(tab$EV1[9:23],tab$EV2[9:23],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV1[24],tab$EV2[24],col="black",bg="#CA1D17", pch=22,cex=1.5)
points(tab$EV1[25:26],tab$EV2[25:26],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV1[27],tab$EV2[27],col="black",bg="#CA1D17", pch=22,cex=1.5)
points(tab$EV1[28:34],tab$EV2[28:34],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV1[35:36],tab$EV2[35:36],col="black",bg="#FEC44F", pch=24,cex=1.5)

plot(tab$EV3,tab$EV4,xlab="eigenvector 3 (5.9%)",ylab="eigenvector 4 (4.1%)",
     xlim=c(-0.35,0.45),ylim=c(-0.25,0.8),
     col="white")
points(tab$EV3[1:2],tab$EV4[1:2],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV3[3:8],tab$EV4[3:8],col="black",bg="#CA1D17", pch=22,cex=1.5)
points(tab$EV3[9:23],tab$EV4[9:23],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV3[24],tab$EV4[24],col="black",bg="#CA1D17", pch=22,cex=1.5)
points(tab$EV3[25:26],tab$EV4[25:26],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV3[27],tab$EV4[27],col="black",bg="#CA1D17", pch=22,cex=1.5)
points(tab$EV3[28:34],tab$EV4[28:34],col="black",bg="#2E509E",pch=21,cex=1.5)
points(tab$EV3[35:36],tab$EV4[35:36],col="black",bg="#FEC44F", pch=24,cex=1.5)
