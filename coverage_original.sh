#!/bin/bash

source Simulation/venv/bin/activate

while read line; do
  ./Simulation/samtools-1.5/samtools sort $line -T $line'temp' -o $line".sorted.bam"
  ./Simulation/samtools-1.5/samtools  index $line".sorted.bam"
  filename=`echo $line | awk -F/ '{print $NF}'`
  geneBody_coverage.py -r $TEAM/genomes_release_89/mm10.HouseKeepingGenes.norchr.bed -i $line".sorted.bam"  -o $filename
done<original_bams.txt
