library(ggplot2)
library(ggpubr)

#read in data
ground_truth_expr<-read.table(gzfile("../data/SupplementaryFigure15.gz"))

#plot data and save
graph<-ggplot(data=ground_truth_expr, aes(x=estimates)) + geom_histogram(binwidth = 10)  +coord_cartesian(ylim=c(0, 100000), xlim=c(0,1500)) + xlab("100 x log2(Ground Truth Expression in TPM + 1)") + ylab("Frequency of Isoforms")
graph
ggsave("../pdfs/SupplementaryFigure15.pdf", plot=last_plot())
ggsave("../pngs/SupplementaryFigure15.png", plot=last_plot())
