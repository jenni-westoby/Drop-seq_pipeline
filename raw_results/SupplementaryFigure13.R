library(scater)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(MASS)
library(ggplot2)
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

save(dropseq_scater_QC, file="Benchmarking_paper/processed_files/data/SupplementaryFigure13_scater_object.RData")

QC_raw<-read.csv("data/Dropseq_results/QC_stats/raw/read_alignment_qc.csv", header=T)
QC_raw<-QC_raw[complete.cases(QC_raw),]

QC_raw_x<-read.csv("data/Dropseq_results/QC_stats/raw/Dropseq_read_alignment_2.csv")
QC_raw<-rbind(QC_raw,QC_raw_x)

QC_raw<-as.data.frame(QC_raw)

#Read number based QC plot
number_reads<-read.csv("Dropbox/Dropseq_num_reads.txt")
number_reads<-as.data.frame(number_reads)
colnames(number_reads)<-c("Filename", "NumReads")
number_reads$Filename<-gsub(".fastq", "", number_reads$Filename)

QC_raw<-merge(QC_raw[,1:5],number_reads)

write.table(QC_raw, file = "Benchmarking_paper/processed_files/data/SupplementaryFigure13_reads_alignment_data.txt")