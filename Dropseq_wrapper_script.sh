# #!/bin/bash
#
# #Rename command line arguments
# path_to_java=$1
# path_to_ref_fasta=$2
# path_to_ref_gtf=$3
#
# ./setup.sh setup
# ./Dropseq_Alignment_Cookbook.sh
# ./RSEM_ref.sh make_ref path_to_ref_gtf path_to_ref_fasta
# ./make_indexes.sh $path_to_ref_gtf $path_to_ref_fasta
# ./RSEM_three_prime_bias.sh $path_to_java
#
# for i in Dropseq_Alignment_Cookbook/demultiplexed_fastqs/*.fastq;
# do
#   ./cell_level_analysis.sh $i
# done


for i in Dropseq_Alignment_Cookbook/small_test_data/*.fastq;
do
  num_jobs=`bjobs | wc -l`
  max_jobs=100

  # #This prevents the number of queued jobs greatly exceeding 30.
  # while [[ $num_jobs -gt $max_jobs ]];
  # do
  #   sleep 100
  #   num_jobs=`bjobs | wc -l`
  # done

  filename=`echo $i | awk -F/ '{print $3}'`
  bsub -n8 -R"span[hosts=1]" -c 99999 -G team_hemberg -q normal -o $TEAM/temp.logs/output.$filename -e $TEAM/temp.logs/error.$filename -R"select[mem>100000] rusage[mem=100000]" -M100000 ./cell_level_analysis.sh $filename
done

# #make clean results matrices
# ./make_matrix.sh make_matrix RSEM
# ./make_matrix.sh make_matrix eXpress
# ./make_matrix.sh make_matrix Kallisto
# ./make_matrix.sh make_matrix Sailfish
# ./make_matrix.sh make_matrix Salmon_align
# ./make_matrix.sh make_matrix Salmon_SMEM
# ./make_matrix.sh make_matrix Salmon_quasi
# ./make_matrix.sh make_matrix ground_truth
# ./make_matrix.sh make_matrix Kallisto_real
# ./clean_data.sh
#
# #move data - some of this should move into setup.sh
# cp Simulation/results_matrices/clean* raw_results/data/
# cp -r Simulation/QC_stats/raw raw_results/data/
# cp -r Simulation/QC_stats/simulated raw_results/data/
#
# #format data to make figures
# cd raw_results
# Rscript Figure3.R
# Rscript SupplementaryFigure13.R
# Rscript SupplementaryFigure14.R
# Rscript SupplementaryFigure15.R
#
# cd ../figures/scripts
# Rscript Figure3.R
# Rscript SupplementaryFigure13.R
# Rscript SupplementaryFigure14.R
# Rscript SupplementaryFigure15.R
#
# #Still need to update the below
# #S14
# #S15
#
# #make figure pdfs
