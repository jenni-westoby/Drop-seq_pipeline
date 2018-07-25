#!/bin/bash
#Quantification pipeline

benchmark(){

  memory=`pwd`
  #If the wrong number of arguments was passed, print an error message
  if [ $# -eq 1 ]
    then
      #If an appropriate argument was passed
      if [ "$1" == "Kallisto" ] || [ "$1" == "Salmon" ] || [ "$1" == "eXpress" ] || [ "$1" == "RSEM" ] || [ "$1" == "Sailfish" ] || [ "$1" == "BRIE" ]; then

        #Make results directory
        mkdir Simulation/$1"_results"

        #If the argument is Salmon, make SMEM, quasi and alignment results directories
        if [ "$1" == "Salmon" ]; then
          cd Simulation/Salmon_results
          mkdir Salmon_Alignment_Results
          mkdir Salmon_SMEM_results
          mkdir Salmon_quasi_results
          cd ../..

          #If the SMEM index doesn't exist, make it
          if [ ! "$(ls -A Simulation/indices/Salmon_SMEM)" ]; then
            start_Salmon_SMEM_index=`date +%s`
            ./Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon index -t Simulation/ref/reference.transcripts.fa -i Simulation/indices/Salmon_SMEM/transcripts_index_SMEM --type fmd -p 8
            stop_Salmon_SMEM_index=`date +%s`
            printf $filename","$((stop_Salmon_SMEM_index-start_Salmon_SMEM_index)) >> Simulation/time_stats/time_Salmon_SMEM_index.csv
          fi

          #If the quasi index doesn't exist, make it
          if [ ! "$(ls -A Simulation/indices/Salmon_quasi)" ]; then
            start_Salmon_quasi_index=`date +%s`
            ./Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon index -t Simulation/ref/reference.transcripts.fa -i Simulation/indices/Salmon_quasi/transcripts_index_quasi --type quasi -k 31 -p 8
            stop_Salmon_quasi_index=`date +%s`
            printf $filename","$((stop_Salmon_quasi_index-start_Salmon_quasi_index)) >> Simulation/time_stats/time_Salmon_quasi_index.csv
          fi
        fi

        if [ "$1" == "Kallisto" ]; then
            #If there is no Kallisto index, make it
            if [ ! "$(ls -A Simulation/indices/Kallisto)" ]; then
                start_kallisto_index=`date +%s`
                ./Simulation/kallisto_linux-v0.43.1/kallisto index -i Simulation/indices/Kallisto/transcripts.idx Simulation/ref/reference.transcripts.fa
                stop_kallisto_index=`date +%s`
                printf $((stop_kallisto_index-start_kallisto_index)) >> Simulation/time_stats/time_kallisto_index.csv
            fi
        fi

        if [ "$1" == "Sailfish" ]; then
              #If there is no index for sailfish, make it
              if [ ! "$(ls -A Simulation/indices/Sailfish)" ]; then
                export LD_LIBRARY_PATH=`pwd`/SailfishBeta-0.10.0_CentOS5/lib:$LD_LIBRARY_PATH
                export PATH=`pwd`/SailfishBeta-0.10.0_CentOS5/bin:$PATH
                start_sailfish_index=`date +%s`
                LC_ALL=C ./Simulation/SailfishBeta-0.10.0_CentOS5/bin/sailfish index -p 8 -t Simulation/ref/reference.transcripts.fa -o Simulation/indices/Sailfish/ -k 31
                stop_sailfish_index=`date +%s`
                printf $filename","$((stop_sailfish_index-start_sailfish_index)) >> Simulation/time_stats/time_sailfish_index.csv
              fi
        fi

        #Run tool on simulated cells. Each cell is submitted as a seperate job.
        cd Simulation/data/simulated
        for i in $(find . -name '*.fq');
        do
          filename=`echo $i |awk -F/ '{print $2}'`
          cd $memory
          #The line below will need to be edited for your LSF job system.
          ./quantify.sh $1 $filename

          num_jobs=`bjobs | wc -l`
          max_jobs=100

          #This prevents the number of queued jobs greatly exceeding 100.
          while [[ $num_jobs -gt $max_jobs ]];
          do
            sleep 100
            num_jobs=`bjobs | wc -l`
          done

        done

        #If a data matrix hasn't been made for the ground_truth results, make it
        if [ ! -f Simulation/results_matrices/ground_truth.csv ]; then
          python ./generate.py ground_truth `pwd` Simulation/data/simulated
          chmod +x ground_truth_TPM.sh
          chmod +x ground_truth_FPKM.sh
          chmod +x ground_truth_Counts.sh
          ./ground_truth_Counts.sh
          ./ground_truth_TPM.sh
          ./ground_truth_FPKM.sh
          rm ground_truth_Counts.sh
          rm ground_truth_TPM.sh
          rm ground_truth_FPKM.sh
        fi

      else
        echo "Inappropriate argument passed. If one argument is passed, the following arguments are accepted: Kallisto, Salmon, eXpress and RSEM."
      fi

  else
    echo "Wrong number of arguments passed. One argument is accepted."
  fi



}

"$@"
