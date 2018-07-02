#!/usr/bin/env bash

while read i;
do

  samtools view Dropseq_Alignment_Cookbook/temp_files/SRR3587500_BSE.bam | grep "XC:Z:"$i | cat header.sam - | samtools view -Sb - > Dropseq_Alignment_Cookbook/demultiplexed_bams/$i".bam"
  
done <$1
