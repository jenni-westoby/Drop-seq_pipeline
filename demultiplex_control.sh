#!/usr/bin/env bash

#loop through array of barcodes, pulling out relevant lines of bamfile for each barcode. Put output into a new directory every 1000 lines for indexing reasons
#Note to self - ran out of time doing 1000 barcodes per file on long job. Need to modify code so splits into 250 files per directory

counter=0

for f in Dropseq_Alignment_Cookbook/temp_files/split_barcodes/*
do
  counter=$((counter+1))

  ./demultiplex.sh $f 

done
