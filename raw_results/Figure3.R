library(scater)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(MASS)
library(ggplot2)
library(SingleCellExperiment)
######################################################################################
# RAW MT QC

dropseq_counts<-read.table("data/clean_Kallisto_real_Counts.txt")

ids<-names(dropseq_counts)
batch<-rep("batch",ncol(dropseq_counts))

anno<- as.data.frame(cbind(batch,ids))
rownames(anno)<-anno$ids

dropseq_scater <- SingleCellExperiment(
  assays = list(counts = as.matrix(dropseq_counts)),
  colData = anno
)

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
isSpike(dropseq_scater, "MT") <- rownames(dropseq_scater) %in% mt_isoforms

dropseq_scater_QC <- calculateQCMetrics(
  dropseq_scater,
  feature_controls = list(MT = isSpike(dropseq_scater, "MT"))
)

mt_reads<-plotPhenoData(
  dropseq_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_MT",
             colour = "batch")
)

dropseq_scater_QC<-dropseq_scater_QC[,dropseq_scater_QC$pct_counts_MT<10]
dropseq_counts<-dropseq_counts[,colnames(dropseq_counts) %in% dropseq_scater_QC$ids]

#####################################################################################
# RAW READS QC
QC_raw<-read.csv("data/raw/read_alignment_qc.csv", header=F)
names(QC_raw)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")
QC_raw<-QC_raw[complete.cases(QC_raw),]

#remove cells based on Figure S13
QC_raw<-QC_raw[QC_raw$NumAlignments<25000 & QC_raw$Unique<15000 & QC_raw$NonUnique<3000 & QC_raw$NumReads<20000,]
QC_raw<-QC_raw[QC_raw$Filename %in% colnames(dropseq_counts),]
rm(list=setdiff(ls(), c("QC_raw")))

##########################################################################################
# SIM MT QC
dropseq_counts<-read.table("data/clean_ground_truth_Counts.txt")

ids<-names(dropseq_counts)
batch<-rep("batch",ncol(dropseq_counts))

anno<-as.data.frame(cbind(batch,ids))
rownames(anno)<-anno$ids

dropseq_scater <- SingleCellExperiment(
  assays = list(counts = as.matrix(dropseq_counts)),
  colData = anno
)

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
isSpike(dropseq_scater, "MT") <- rownames(dropseq_scater) %in% mt_isoforms

dropseq_scater_QC <- calculateQCMetrics(
  dropseq_scater,
  feature_controls = list(MT = isSpike(dropseq_scater, "MT"))
)

mt_reads<-scater::plotPhenoData(
  dropseq_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_MT",
             colour = "batch")
)

dropseq_scater_QC<-dropseq_scater_QC[,dropseq_scater_QC$pct_counts_MT<10]
dropseq_counts<-dropseq_counts[,colnames(dropseq_counts) %in% dropseq_scater_QC$ids]

##################################################################################################
# SIM READS QC

QC_sim<-read.csv("data/simulated/read_alignment_qc.csv", header=F)
names(QC_sim)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")
QC_sim<-as.data.frame(QC_sim)
QC_sim$Filename<-sub(".fastq", "", QC_sim$Filename)
QC_sim<-QC_sim[QC_sim$Filename %in% QC_raw$Filename,]
QC_sim<-QC_sim[QC_sim$Filename %in% colnames(dropseq_counts),]

QC_sim<-QC_sim[QC_sim$NumAlignments>1000 & QC_sim$Unique>1000 & QC_sim$NonUnique<500,]

rm(list=setdiff(ls(), c("QC_sim")))

################################################################################
#LOAD DATA AND PROCESS IT

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


#load data
ground_truth<-data_processing("data/clean_ground_truth_TPM.txt", QC_sim)
RSEM<-data_processing("data/clean_RSEM_TPM.txt", QC_sim)
Salmon_align<-data_processing("data/clean_Salmon_align_TPM.txt", QC_sim)
Salmon_quasi<-data_processing("data/clean_Salmon_quasi_TPM.txt", QC_sim)
Salmon_SMEM<-data_processing("data/clean_Salmon_SMEM_TPM.txt", QC_sim)
Sailfish<-data_processing("data/clean_Sailfish_TPM.txt", QC_sim)
eXpress<-data_processing("data/clean_eXpress_TPM.txt", QC_sim)
Kallisto<-data_processing("data/clean_Kallisto_TPM.txt", QC_sim)

################################################################################
# SPEARMAN'S GRAPH

#Function to find correlation
correlation<-function(x,y) {
  (diag(cor(y,x,method="spearman")))
}

#Find Spearman's rho for each method
RSEM_cor<-correlation(RSEM,ground_truth)
Salmon_align_cor<-correlation(Salmon_align, ground_truth)
Salmon_quasi_cor<-correlation(Salmon_quasi, ground_truth)
Salmon_SMEM_cor<-correlation(Salmon_SMEM, ground_truth)
Sailfish_cor<-correlation(Sailfish, ground_truth)
eXpress_cor<-correlation(eXpress, ground_truth)
Kallisto_cor<-correlation(Kallisto, ground_truth)



#get data into right format
spearmans_data<-melt(rbind(RSEM_cor,Salmon_align_cor, Salmon_quasi_cor, Salmon_SMEM_cor, Sailfish_cor, eXpress_cor, Kallisto_cor))

##############################################################################################
# NRMSE GRAPH

#find NRMSE
library(hydroGOF)

RSEM_nmrse<-nrmse(log2(RSEM+1),log2(ground_truth+1))
Salmon_align_nmrse<-nrmse(log2(Salmon_align+1),log2(ground_truth+1))
Salmon_quasi_nmrse<-nrmse(log2(1+Salmon_quasi),log2(ground_truth+1))
Salmon_SMEM_nmrse<-nrmse(log2(Salmon_SMEM+1),log2(ground_truth+1))
Sailfish_nmrse<-nrmse(log2(1+Sailfish),log2(ground_truth+1))
eXpress_nmrse<-nrmse(log2(eXpress+1),log2(ground_truth+1))
Kallisto_nmrse<-nrmse(log2(Kallisto+1),log2(ground_truth+1))

#get data into right format
nrmse_data<-melt(rbind(RSEM_nmrse,Salmon_align_nmrse, Salmon_quasi_nmrse, Salmon_SMEM_nmrse, Sailfish_nmrse, eXpress_nmrse, Kallisto_nmrse))

###############################################################################################
# PRECISION GRAPH

#Function that returns precision
make_precision<-function(ground_truth, tool_estimates, threshold_unexpr){
  TP<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates>threshold_unexpr])
  #print(TP)
  FP<-length(ground_truth[ground_truth<=threshold_unexpr & tool_estimates>threshold_unexpr])
  return(TP/(TP+FP))
}

#Function that returns precision value per cell
return_precision_per_cell<-function(truth_input_data, estimate_input_data){

  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_precision(truth_input_data[,j], estimate_input_data[,j], 0)
  }

  return(do.call(rbind,results))

}

#Find precision for each method
RSEM_precision<-return_precision_per_cell(ground_truth,RSEM)[,1]
Salmon_align_precision<-return_precision_per_cell(ground_truth,Salmon_align)[,1]
Salmon_quasi_precision<-return_precision_per_cell(ground_truth,Salmon_quasi)[,1]
Salmon_SMEM_precision<-return_precision_per_cell(ground_truth,Salmon_SMEM)[,1]
Sailfish_precision<-return_precision_per_cell(ground_truth,Sailfish)[,1]
eXpress_precision<-return_precision_per_cell(ground_truth,eXpress)[,1]
Kallisto_precision<-return_precision_per_cell(ground_truth,Kallisto)[,1]

precision_data<-melt(rbind(RSEM_precision,Salmon_align_precision,  Salmon_quasi_precision, Salmon_SMEM_precision, Sailfish_precision, eXpress_precision, Kallisto_precision))

#################################################################################################
# RECALL PLOTS

#Function that returns precision
make_recall<-function(ground_truth, tool_estimates, threshold_unexpr){
  TP<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates>threshold_unexpr])
  #print(TP)
  FN<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates<=threshold_unexpr])
  return(TP/(TP+FN))
}

#Function that returns precision value per cell
return_recall_per_cell<-function(truth_input_data, estimate_input_data){

  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_recall(truth_input_data[,j], estimate_input_data[,j], 0)
  }

  return(do.call(rbind,results))

}

#Find precision for each method
RSEM_recall<-return_recall_per_cell(ground_truth,RSEM)[,1]
Salmon_align_recall<-return_recall_per_cell(ground_truth,Salmon_align)[,1]
Salmon_quasi_recall<-return_recall_per_cell(ground_truth,Salmon_quasi)[,1]
Salmon_SMEM_recall<-return_recall_per_cell(ground_truth,Salmon_SMEM)[,1]
Sailfish_recall<-return_recall_per_cell(ground_truth,Sailfish)[,1]
eXpress_recall<-return_recall_per_cell(ground_truth,eXpress)[,1]
Kallisto_recall<-return_recall_per_cell(ground_truth,Kallisto)[,1]

recall_data<-melt(rbind(RSEM_recall,Salmon_align_recall,  Salmon_quasi_recall, Salmon_SMEM_recall, Sailfish_recall, eXpress_recall, Kallisto_recall))

###############################################################
# F1 PLOTS

find_F1<-function(precision,recall){
  F1<-2*((precision*recall)/(precision + recall))
  return(F1)
}

RSEM_F1<-find_F1(RSEM_precision,RSEM_recall)
Salmon_align_F1<-find_F1(Salmon_align_precision, Salmon_align_recall)
Salmon_quasi_F1<-find_F1(Salmon_quasi_precision, Salmon_quasi_recall)
Salmon_SMEM_F1<-find_F1(Salmon_SMEM_precision, Salmon_SMEM_recall)
Sailfish_F1<-find_F1(Sailfish_precision, Sailfish_recall)
eXpress_F1<-find_F1(eXpress_precision, eXpress_recall)
Kallisto_F1<-find_F1(Kallisto_precision, Kallisto_recall)

F1_data<-melt(rbind(RSEM_F1,Salmon_align_F1, Salmon_quasi_F1, Salmon_SMEM_F1, Sailfish_F1, eXpress_F1, Kallisto_F1))

F1_data$Var2<-as.character(F1_data$Var2)
precision_data$Var2<-as.character(precision_data$Var2)
recall_data$Var2<-as.character(recall_data$Var2)

F1_data<-cbind(F1_data, statistic=rep("F1", nrow(F1_data)))
precision_data<-cbind(precision_data, statistic=rep("precision", nrow(precision_data)))
recall_data<-cbind(recall_data, statistic=rep("recall", nrow(recall_data)))
spearmans_data<-cbind(spearmans_data, statistic=rep("spearmans", nrow(spearmans_data)))
nrmse_data<-cbind(nrmse_data, statistic=rep("nrmse", nrow(nrmse_data)))

ggplot_results<-rbind(spearmans_data,nrmse_data,F1_data,precision_data,recall_data)
write.table(ggplot_results, "../figures/data/Figure3.txt")
