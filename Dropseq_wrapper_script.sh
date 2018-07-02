#!/bin/bash

#Rename command line arguments
path_to_java=$1
path_to_ref_fasta=$2
path_to_ref_gtf=$3

./setup.sh setup
./Dropseq_Alignment_Cookbook.sh
./RSEM_ref.sh make_ref path_to_ref_gtf path_to_ref_fasta
./make_indexes.sh
./RSEM_three_prime_bias.sh $path_to_java

for i in Dropseq_Alignment_Cookbook/demultiplexed_fastqs/*_1.fastq;
do
  ./cell_level_analysis.sh $i
done

#make clean results matrices
./make_matrix.sh make_matrix RSEM
./make_matrix.sh make_matrix eXpress
./make_matrix.sh make_matrix Kallisto
./make_matrix.sh make_matrix Sailfish
./make_matrix.sh make_matrix Salmon_align
./make_matrix.sh make_matrix Salmon_SMEM
./make_matrix.sh make_matrix Salmon_quasi
./make_matrix.sh make_matrix ground_truth
./make_matrix.sh make_matrix Kallisto_real
./clean_data.sh

#move data - some of this should move into setup.sh
cp Simulation/results_matrices/clean* raw_results/data/

#format data to make figures
cd raw_results
Rscript Figure3.R
Rscript SupplementaryFigure13.R
Rscript SupplementaryFigure14.R
Rscript SupplementaryFigure15.R

cd ../figures/scripts
Rscript Figure3.R
Rscript SupplementaryFigure13.R
Rscript SupplementaryFigure14.R
Rscript SupplementaryFigure15.R

#Still need to update the below
#S14
#S15

#make figure pdfs
