#!/bin/bash

filename=`echo $1 |awk -F. '{print $1}'`

./quality_control.sh qualitycontrol $filename Dropseq_Alignment_Cookbook/demultiplexed_fastqs '.fastq' 'raw'
./simulate.sh simulate $filename Dropseq_Alignment_Cookbook/demultiplexed_fastqs
./quality_control.sh qualitycontrol $filename Simulation/data/simulated '.fq' "simulated"
./quantify.sh Kallisto $filename
./quantify.sh eXpress $filename
./quantify.sh Salmon $filename
./quantify.sh RSEM $filename
./quantify.sh Sailfish $filename
./quantify_real_data.sh Kallisto $filename
