library(scater)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(MASS)
library(ggplot2)

######################################################################################
# RAW MT QC
dropseq_counts<-read.table("data/clean_Kallisto_real_Counts.txt")

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
QC_raw<-read.csv("../Simulation/QC_stats/raw/read_alignment_qc.csv", header=T)

#remove based on SupplementaryFigure13
QC_raw<-QC_raw[QC_raw$NumAlignments<25000 & QC_raw$Unique<15000 & QC_raw$NonUnique<3000 & QC_raw$NumReads<20000,]
QC_raw<-QC_raw[QC_raw$Filename %in% colnames(dropseq_counts),]
rm(list=setdiff(ls(), c("QC_raw")))

##########################################################################################
# SIM MT QC
dropseq_counts<-read.table("data/clean_ground_truth_Counts.txt")

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

save(dropseq_scater_QC, file="../figures/data/SupplementaryFigure14_scater_object.RData")

#####################################################################################
# SIM READS QC
QC_sim<-read.csv("../Simulation/QC_stats/simulated/read_alignment_qc.csv", header=T)
QC_sim<-as.data.frame(QC_sim)
QC_sim$Filename<-sub(".fastq", "", QC_sim$Filename)
QC_sim<-QC_sim[QC_sim$Filename %in% QC_raw$Filename,]

write.table(QC_sim, file = "../figures/data/SupplementaryFigure14_reads_alignment_data.txt")
