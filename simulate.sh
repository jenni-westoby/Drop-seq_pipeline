#!/bin/bash
#RSEM simulation pipeline
#Note to self check dir names. Need to sort out index. Pass an arg to where the raw data is stored instead of requiring it to be moved?

#Function which submits RSEM simulations as LSF style jobs. Takes one arg, the directory where the data to be simulated is stored.
run_simulations(){
  if [ $# -ne 1 ]
    then
      echo "Incorrect number of arguments supplied. One argument should be passed to this function, the path to the directory in which the data is stored."
      exit 1
  fi
  memory=`pwd`
  cd $1
  for i in $(find . -name '*.fastq');
  do
    base=`echo $i |awk -F/ '{print $2}'`
    filename=`echo $base |awk -F_ '{print $1}'`
    cd $memory
    #The line below will need to be edited for your LSF job system.
    bsub -n8 -R"span[hosts=1]" -c 99999 -G team_hemberg -q normal -o $TEAM/temp.logs/output.$filename -e $TEAM/temp.logs/error.$filename -R"select[mem>100000] rusage[mem=100000]" -M100000 simulate $i $1

    num_jobs=`bjobs | wc -l`
    max_jobs=100

    #This prevents the number of queued jobs greatly exceeding 100.
    while [[ $num_jobs -gt $max_jobs ]];
    do
      sleep 100
      num_jobs=`bjobs | wc -l`
    done
  done
}

#Function which performs RSEM simulations. Takes 2 args, the filename of the cell and the directory in which it is stored.
simulate() {


 #Make filename strings
 filename=$1
 raw_data_dir=${2%/}

 lines="$(wc -l $raw_data_dir/$filename'.fastq' | awk '{print $1}')"
 reads="$(echo $((lines / 4)))"

 #Use RSEM to calculate expression
 ./Simulation/RSEM-1.3.0/rsem-calculate-expression --star\
       --star-path Simulation/STAR/bin/Linux_x86_64/ \
       -p 8 \
                   --estimate-rspd \
                   --append-names \
                   --output-genome-bam \
       --single-cell-prior --calc-pme \
                   $raw_data_dir/$filename'.fastq' \
                   Simulation/ref/reference Simulation/data/temp/$filename
                   #--fragment-length-mean 1600 --fragment-length-sd 50 \

  #extract first number of third line of filename.theta, which is an estimate of the portion of reads due to background noise
  background_noise=`sed '3q;d' Simulation/data/temp/$filename".stat"/$filename".theta" | awk '{print $1}'`


  #Simulate reads
  ./Simulation/RSEM-1.3.0/rsem-simulate-reads Simulation/ref/reference Simulation/data/temp/$filename".stat"/$filename".model" \
                                         Simulation/data/temp/$filename".isoforms.results" $background_noise $reads Simulation/data/simulated/$filename \
                                         --seed 0

  #Tidy up
  #rm -r Simulation/data/temp/$filename*
}

export -f simulate

"$@"
