#!/bin/bash

filename=`echo $1 |awk -F_ '{print $1}'`

./quality_control.sh quality_control $filename Dropseq_Alignment_Cookbook/demultiplexed_fastqs '1.fastq' 'raw'
./simulate.sh simulate $filename Dropseq_Alignment_Cookbook/demultiplexed_fastqs
./quality_control.sh qualitycontrol $filename Simulation/data/simulated '1.fq' "simulated"
./quantify.sh Kallisto $filename
./quantify.sh eXpress $filename
./quantify.sh Salmon $filename
./quantify.sh RSEM $filename
./quantify.sh Sailfish $filename
./quantify_real_data.sh Kallisto $filename
