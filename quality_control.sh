#!/bin/bash

#Function which produces quality control statistics based on read and alignment
#quality. Takes 2 or 3 args. The first arg is the path to the reference gtf file,
#the second is the path to the reference fasta, the third is the path to the raw
#data directory. If 2 args are passed quality control is performed on the
#simulated data, if 3 are passed it is performed on the raw data.

QC() {

  #If there are two args, check that the STAR index exists and if not make it, then produce quality control statistics for simulated data
  if [ $# == 2 ]; then
    if [ ! "$(ls -A Simulation/indices/STAR)" ]; then
      Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir Simulation/indices/STAR --genomeFastaFiles $2  --sjdbGTFfile $1
    fi

    #If there is already a csv file of QC statistics, remove it, then create the header for the csv file
    if [ "$(ls -A Simulation/QC_stats/simulated)" ]; then
      rm Simulation/QC_stats/simulated/read_alignment_qc.csv
    fi
    echo "Filename,Unique,NonUnique,Unmapped,NumAlignments,NumReads" >> Simulation/QC_stats/simulated/read_alignment_qc.csv

    memory=`pwd`
    cd Simulation/data/simulated
    for i in $(find . -name '*.fq');
    do
      base=`echo $i |awk -F/ '{print $2}'`
      filename=`echo $base | rev | cut -d . -f2- | rev`
      end=`echo $base |awk -F. '{print $(NF-2)"."$(NF-1)"."$(NF)}' | awk -F. '{print $NF}'`
      cd $memory
      qualitycontrol $filename Simulation/data/simulated ${end%/} "simulated"

      num_jobs=`bjobs | wc -l`
      max_jobs=100

      #This prevents the number of queued jobs greatly exceeding 100.
      while [[ $num_jobs -gt $max_jobs ]];
      do
        sleep 100
        num_jobs=`bjobs | wc -l`
      done
    done

  #If there are three args, check that the STAR index exists and if not make it, then produce quality control statistics for raw data
  elif [ $# == 3 ]; then
    if [ ! "$(ls -A Simulation/indices/STAR)" ]; then
      Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir Simulation/indices/STAR --genomeFastaFiles $2  --sjdbGTFfile $1
    fi

    #If there is already a csv file of QC statistics, remove it, then create the header for the csv file
    if [ "$(ls -A Simulation/QC_stats/raw)" ]; then
      rm Simulation/QC_stats/raw/read_alignment_qc.csv
    fi
    echo "Filename,Unique,NonUnique,Unmapped,NumAlignments,NumReads" >> Simulation/QC_stats/raw/read_alignment_qc.csv

    memory=`pwd`
    cd $3
    for i in $(find . -name '*.fastq*');
    do
      base=`echo $i |awk -F/ '{print $2}'`
      filename=`echo $base | rev | cut -d . -f2- | rev`
      end=`echo $base |awk -F. '{print $(NF-2)"."$(NF-1)"."$(NF)}' | awk -F. '{print $NF}'`
      data_dir=${3%/}
      cd $memory
      qualitycontrol $filename $data_dir ${end%/} "raw"

      num_jobs=`bjobs | wc -l`
      max_jobs=100

      #This prevents the number of queued jobs greatly exceeding 100.
      while [[ $num_jobs -gt $max_jobs ]];
      do
        sleep 100
        num_jobs=`bjobs | wc -l`
      done
    done

  #Otherwise, print an error message and exit
  else
      echo "Incorrect number of arguments supplied. This function takes 2 or 3 arguments. The first argument is the path to the reference gtf file, the second is the path to the reference fasta, the third is the path to the raw data directory. If 2 args are passed quality control is performed on the simulated data, if 3 are passed it is performed on the raw data."
      exit 1
  fi
}


#Function which records quality statistics for reads
qualitycontrol() {

  #Name input arguments
  filename=$1
  data_dir=$2
  end=$3
  subdir=$4

  gz=`echo $end |awk -F. '{print $3}'`
  fast=`echo $end |awk -F. '{print $2}'`

  #If the file is gzipped, create an uncompressed copy, perform the analysis, then delete the uncompressed copy
  if [ "$gz" == gz ]; then
    gunzip -c $data_dir/$filename"_"$end > $data_dir/$filename'_1.fq'
    gunzip -c $data_dir/$filename"_2."$fast"."$gz > $data_dir/$filename'_2.fq'
    STAR_and_RSeQC $filename $data_dir '1.fq' $subdir
    rm $data_dir/$filename'_1.fq'
    rm $data_dir/$filename'_2.fq'

  else
    STAR_and_RSeQC $filename $data_dir $end $subdir
  fi
}

STAR_and_RSeQC() {

  #Name input arguments
  filename=$1
  data_dir=$2
  end=$3
  subdir=$4

  #Run STAR
  Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir Simulation/indices/STAR --readFilesIn $data_dir/$filename$end --outFileNamePrefix Simulation/bamfiles/$subdir/$filename --outSAMtype BAM SortedByCoordinate

  #Use bam_stat from the RSeQC package to find alignment statistics
  source Simulation/venv/bin/activate
  bam_stat.py -i Simulation/bamfiles/$subdir/$filename"Aligned.sortedByCoord.out.bam" >> Simulation/QC_stats/"temp_"$subdir/$filename"bam_stat"
  deactivate

  #Find number of unique mapping, non-unique mapping and unmapped reads
  Unique=`grep "mapq >= mapq_cut (unique):" Simulation/QC_stats/"temp_"$subdir/$filename"bam_stat" | awk '{print $NF}'`
  NonUnique=`grep "mapq < mapq_cut (non-unique)" Simulation/QC_stats/"temp_"$subdir/$filename"bam_stat" | awk '{print $NF}'`
  Unmapped=`grep "Unmapped reads:" Simulation/QC_stats/"temp_"$subdir/$filename"bam_stat" | awk '{print $NF}'`
  NumAlignments=`grep "Total records:" Simulation/QC_stats/"temp_"$subdir/$filename"bam_stat" | awk '{print $NF}'`

  lines="$(wc -l $data_dir/$filename$end | awk '{print $1}')"
  reads="$(echo $((lines / 4)))"

  #Append stats to csv file
  echo $filename","$Unique","$NonUnique","$Unmapped","$NumAlignments","$reads >> Simulation/QC_stats/$subdir/read_alignment_qc.csv

  #Tidy up
  rm -r Simulation/QC_stats/"temp_"$subdir/$filename*
  rm -r Simulation/bamfiles/$subdir/$filename*
}

export -f STAR_and_RSeQC
export -f qualitycontrol

"$@"
