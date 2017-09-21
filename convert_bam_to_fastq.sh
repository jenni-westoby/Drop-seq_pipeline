#!/usr/bin/env bash

path_to_data=$1
path_to_java=$2
path_to_picard_jar=$3
outdir=$4

for file in $path_to_data/*
do
  filename=`echo $file | awk -F/ '{print $NF}' | awk -F.bam '{print $1}'`
  $path_to_java -XX:MaxHeapSize=1000m -jar $path_to_picard_jar SamToFastq \
        I=$file \
        FASTQ=$outdir/$filename"_1.fastq" \
done
