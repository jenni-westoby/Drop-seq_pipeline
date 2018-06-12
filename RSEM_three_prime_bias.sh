#!/bin/bash

cd Dropseq_Alignment_Cookbook

#Convert bamfile to fastq
$path_to_java -XX:MaxHeapSize=1000m -jar picard.jar SamToFastq \
I=temp_files/SRR3587500_BSE.bam \
FASTQ=final_SRR3587500.fastq \


cd ..

#RSEM
./Simulation/RSEM/rsem-calculate-expression --star\
      --star-path Simulation/STAR/bin/Linux_x86_64 \
      -p 8 \
      --fragment-length-mean 1600 --fragment-length-sd 50 \
      --append-names \
      --calc-pme \
      Dropseq_Alignment_Cookbook/final_SRR3587500.fastq \
      Simulation/ref/reference Dropseq_Alignment_Cookbook/final_SRR3587500
