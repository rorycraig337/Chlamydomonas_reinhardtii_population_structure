setwd("/Users/rory/Documents/research/projects/chlamy_biogeography/github/identity_by_descent/NA1/")
library(ggplot2)
library(reshape2)

WG <- read.table("NA1.IBD_tracts.within_genome.tsv", header = TRUE)
WG_long <- melt(WG)
WG_long

pl <- ggplot(WG_long,aes(x=value, fill=variable)) + geom_density(alpha=0.35) + theme_bw()
pl + scale_fill_manual(values=c("#1f78b4", "#6a3d9a"))

require(gridExtra)
combined <- read.table("NA1_IBD_summary.txt", header = FALSE)
cor.test(combined$V3,combined$V4)
cor.test(combined$V3,combined$V5)
cor.test(combined$V4,combined$V5)
combined_long <- melt(combined)
combined_long

m <- ggplot(combined_long,aes(x=value, fill=variable)) + geom_density(alpha=0.35) + theme_bw()
m + scale_fill_manual(values=c("#33a02c", "#1f78b4", "#6a3d9a"))

plot1 <- m + scale_fill_manual(values=c("#33a02c", "#1f78b4", "#6a3d9a")) +  theme(legend.position="none")
plot2 <- pl + scale_fill_manual(values=c("#1f78b4", "#6a3d9a")) + theme(legend.position="none")
grid.arrange(plot1, plot2, ncol=2)

wfield <- read.table("within_field_cumulative.txt", header = FALSE)
bfield <- read.table("Farnham_MacDonald_cumulative.txt", header = FALSE) 
wilcox.test(wfield$V3,bfield$V3)

wtime <- read.table("within_time_cumulative.txt", header = FALSE)
btime <- read.table("Farnham_Farnham_cumulative.txt", header = FALSE) 
wilcox.test(wtime$V3,btime$V3)

wtime <- read.table("within_time_cumulative.txt", header = FALSE)
btime <- read.table("Farnham_Farnham_cumulative.txt", header = FALSE) 
wilcox.test(wtime$V3,btime$V3)

within_MA_QC <- read.table("within_MA_QC_cumulative.txt", header = FALSE)
between_MA_QC <- read.table("MA_QC_cumulative.txt", header = FALSE)
wilcox.test(within_MA_QC$V3,between_MA_QC$V3)

far93 <- read.table("Farnham_93_94_cumulative.txt", header = FALSE)
mac <- read.table("MacDonald_cumulative.txt", header = FALSE) 
far16 <- read.table("Farnham_2016_cumulative.txt", header = FALSE)
wilcox.test(far93$V3,mac$V3)
wilcox.test(far93$V3,far16$V3)

cai <- read.table("introgression_cohort_average.txt", header = FALSE)
cor.test(cai$V4,cai$V5)
