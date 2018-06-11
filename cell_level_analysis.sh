#!/bin/bash

./quality_control.sh quality_control $1 Dropseq_Alignment_Cookbook/demultiplexed_fastqs '_1.fq' 'raw'
./simulate.sh simulate $1 Dropseq_Alignment_Cookbook/demultiplexed_fastqs
./quality_control.sh qualitycontrol $1 Simulation/data/simulated '_1.fq' "simulated"
./benchmark.sh quantify Kallisto $1
./benchmark.sh quantify eXpress $1
./benchmark.sh quantify Salmon $1
./benchmark.sh quantify RSEM $1
./benchmark.sh quantify Sailfish $1
./quantify_real_data.sh Kallisto $1
