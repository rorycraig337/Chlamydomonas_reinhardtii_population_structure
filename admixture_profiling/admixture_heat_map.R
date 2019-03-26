setwd("/Users/rory/Documents/research/projects/chlamy_biogeography/github/admixture_profiling/")
library(ape)
library(phytools)
library(grid)
library(ggplot2)
library(ggdendro)
library(reshape2)

totals <- read.table("admixed_totals.txt", header = FALSE)
totals
Mb = totals[,3] / 1000000
Mb
barplot(Mb, ylim = c(0,25))

tree <- read.tree("tree.nwk")
chronogram <- chronos(tree)
tree_dendogram <- as.dendrogram(as.hclust.phylo(chronogram))
dendro.plot <- ggdendrogram(data = tree_dendogram, rotate = TRUE)
print(dendro.plot)

data <- read.csv("chromosome_17.20kb.R.csv", comment.char="#", na.strings=c("NA"))
data$isolates

data.matrix <- as.matrix(data[,2:ncol(data)])
rnames <- data[,1]
rnames

data.long <- melt(data, id = c("isolates"))
data.long$isolates
order <- c(4,6,8,28,9,5,7,25,18,1,2,3,11,32,13,26,27,29,33,24,31,15,14,16,21,30,34,35,22,10,12,19,20,17,23)
data.long$isolates = factor(data.long$isolates,levels(data$isolates)[order])
heatmap.plot <- ggplot(data = data.long, aes(x = variable, y = isolates)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(high = "#e31a1c", low = "#1f78b4", na.value="grey", limit = c(0,1), midpoint = 0.5) +
  theme(axis.text.y = element_text(size = 6),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "top")

grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 0.98))
dev.off()
plot(x_values,admix,type="l", lwd=1.5, cex.axis=1,frame=F,xlab="Mb",ylab = "admix(%)",col="#00441b",xlim=c(0,9.730733),ylim=c(0,45),xaxt="n")
axis(1,at=c(0,3.500558), labels=c("",""), lwd.ticks = 0)
axis(1,at=seq(0,3.5, by=0.5), lwd=0, lwd.ticks = 1, labels=FALSE, tck=-0.02)
