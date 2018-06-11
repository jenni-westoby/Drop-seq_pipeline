library(scater)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(MASS)
library(ggplot2)
library(data.table)
######################################################################################
# RAW MT QC

dropseq_counts<-read.table("data/Dropseq_results/clean_Kallisto_real_Counts.txt")

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
ids<-names(dropseq_counts)
batch<-rep("batch",1000)

anno<-new("AnnotatedDataFrame", as.data.frame(cbind(batch,ids)))
rownames(anno)<-anno$ids
dropseq_scater <- scater::newSCESet(
  countData = dropseq_counts,
  phenoData = anno
)

dropseq_scater_QC <- scater::calculateQCMetrics(
  dropseq_scater,
  feature_controls = list(MT = mt_isoforms)
)

mt_reads<-scater::plotPhenoData(
  dropseq_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_MT",
             colour = "batch")
)

dropseq_scater_QC<-dropseq_scater_QC[,dropseq_scater_QC$pct_counts_feature_controls_MT<10]
dropseq_counts<-dropseq_counts[,colnames(dropseq_counts) %in% phenoData(dropseq_scater_QC)$ids]
#####################################################################################
# RAW READS QC
QC_raw<-read.csv("data/Dropseq_results/QC_stats/raw/read_alignment_qc.csv", header=T)
QC_raw<-QC_raw[complete.cases(QC_raw),]

QC_raw_x<-read.csv("data/Dropseq_results/QC_stats/raw/Dropseq_read_alignment_2.csv")
QC_raw<-rbind(QC_raw,QC_raw_x)

QC_raw<-as.data.frame(QC_raw)

#Alignment based QC plots
Unique<-ggplot(data=QC_raw, aes(x=reorder(Filename,Unique), y=Unique/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() +xlab(" ") + ggtitle("\nNumber of Uniquely Mapping Reads")
NonUnique<-ggplot(data=QC_raw, aes(x=reorder(Filename,NonUnique), y=NonUnique/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + xlab(" ") + ggtitle("Number of Non-Uniquely Mapping Reads")
Unmapped<-ggplot(data=QC_raw, aes(x=reorder(Filename,Unmapped), y=Unmapped)) + geom_point(stat="identity") + ylab("Number of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + ylim(0,100) + xlab(" ") + ggtitle("Number of Unmapped Reads")
NumAlign<-ggplot(data=QC_raw, aes(x=reorder(Filename,NumAlignments), y=NumAlignments/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + xlab(" ") + ggtitle("\nNumber of Alignments")

#Read number based QC plot
number_reads<-read.csv("Dropbox/Dropseq_num_reads.txt")
number_reads<-as.data.frame(number_reads)
colnames(number_reads)<-c("Filename", "NumReads")
NumReads<-ggplot(data=number_reads, aes(x=reorder(Filename,NumReads), y=NumReads/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_text(size=14)) + scale_x_discrete() + xlab(" ") + ggtitle("\nNumber of Reads")

ggarrange(mt_reads, NumReads, NumAlign, Unique, NonUnique, Unmapped, ncol = 2,nrow=3, labels = c("A","B","C","D","E","F"))
ggsave("Pictures/Benchmarking_paper_Figures/Supplementary_Figures/Supplementary_3_Dropseq_raw_QC.pdf", plot=last_plot(),height = (6.04*1.5), width=(8.57*1.5))

#remove based on graphs
QC_raw<-QC_raw[QC_raw$NumAlignments<25000 & QC_raw$Unique<15000 & QC_raw$NonUnique<3000 & QC_raw$NumReads<20000,]
QC_raw<-QC_raw[QC_raw$Filename %in% colnames(dropseq_counts),]
rm(list=setdiff(ls(), c("QC_raw")))
##########################################################################################
# SIM MT QC
dropseq_counts<-read.table("data/Dropseq_results/clean_ground_truth_Counts.txt")

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
ids<-names(dropseq_counts)
batch<-rep("batch",1000)

anno<-new("AnnotatedDataFrame", as.data.frame(cbind(batch,ids)))
rownames(anno)<-anno$ids
dropseq_scater <- scater::newSCESet(
  countData = dropseq_counts,
  phenoData = anno
)

dropseq_scater_QC <- scater::calculateQCMetrics(
  dropseq_scater,
  feature_controls = list(MT = mt_isoforms)
)

mt_reads<-scater::plotPhenoData(
  dropseq_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_MT",
             colour = "batch")
)

dropseq_scater_QC<-dropseq_scater_QC[,dropseq_scater_QC$pct_counts_feature_controls_MT<10]
dropseq_counts<-dropseq_counts[,colnames(dropseq_counts) %in% phenoData(dropseq_scater_QC)$ids]

##################################################################################################
# SIM READS QC

QC_sim<-read.csv("data/Dropseq_results/QC_stats/simulated/read_alignment_qc.csv", header=T)
QC_sim<-as.data.frame(QC_sim)
QC_sim$Filename<-sub(".fastq", "", QC_sim$Filename)
QC_sim<-QC_sim[QC_sim$Filename %in% QC_raw$Filename,]

#Alignment based QC plots
Unique<-ggplot(data=QC_sim, aes(x=reorder(Filename,Unique), y=Unique/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() +xlab(" ") + ggtitle("\nNumber of Uniquely Mapping Reads")
NonUnique<-ggplot(data=QC_sim, aes(x=reorder(Filename,NonUnique), y=NonUnique/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + xlab(" ") + ggtitle("Number of Non-Uniquely Mapping Reads")
Unmapped<-ggplot(data=QC_sim, aes(x=reorder(Filename,Unmapped), y=Unmapped)) + geom_point(stat="identity") + ylab("Number of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + ylim(0,100) + xlab(" ") + ggtitle("Number of Unmapped Reads")
NumAlign<-ggplot(data=QC_sim, aes(x=reorder(Filename,NumAlignments), y=NumAlignments/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + xlab(" ") + ggtitle("\nNumber of Alignments")

#Read number based QC plot
number_reads<-read.csv("Dropbox/Dropseq_num_reads.txt")
number_reads<-as.data.frame(number_reads)
colnames(number_reads)<-c("Filename", "NumReads")
NumReads<-ggplot(data=number_reads, aes(x=reorder(Filename,NumReads), y=NumReads/1000)) + geom_point(stat="identity") + ylab("Thousands of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_text(size=14)) + scale_x_discrete() + xlab(" ") + ggtitle("\nNumber of Reads")

ggarrange(mt_reads, NumReads, NumAlign, Unique, NonUnique, Unmapped, ncol = 2,nrow=3, labels = c("A","B","C","D","E","F"))
ggsave("Pictures/Benchmarking_paper_Figures/Supplementary_Figures/Supplementary_3_Dropseq_sim_QC.pdf", plot=last_plot(),height = (6.04*1.5), width=(8.57*1.5))

QC_sim<-QC_sim[QC_sim$NumAlignments>1000 & QC_sim$Unique>1000 & QC_sim$NonUnique<500,]

rm(list=setdiff(ls(), c("QC_sim")))

#########################################################################

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter){
  results<-read.table(path)
  colnames(results)<-sub(".fastq", "", colnames(results))
  colnames(results)<-sub(".fq", "", colnames(results))
  results<-results[,colnames(results) %in% filter$Filename]
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}


ground_truth<-data_processing("data/Dropseq_results/clean_ground_truth_TPM.txt", QC_sim)
setDT(ground_truth, keep.rownames = TRUE)[]
  
cells<-colnames(ground_truth)
  
ground_truth_expr<-ground_truth %>% gather(cells[2:ncol(ground_truth)], key="cell", value="estimates")
ground_truth_expr$estimates<-100*log2(ground_truth_expr$estimates +1)

write.table(ground_truth_expr, gzfile("Benchmarking_paper/processed_files/data/SupplementaryFigure15.gz"))