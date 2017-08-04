#!/bin/bash
#Quantification pipeline

benchmark(){

  memory=`pwd`
  #If the wrong number of arguments was passed, print an error message
  if [ $# -eq 1 ]
    then
      #If an appropriate argument was passed
      if [ "$1" == "Kallisto" ] || [ "$1" == "Salmon" ] || [ "$1" == "eXpress" ] || [ "$1" == "RSEM" ] || [ "$1" == "Sailfish" ]; then

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
            ./Simulation/Salmon-0.7.2_linux_x86_64/bin/salmon index -t Simulation/ref/reference.transcripts.fa -i Simulation/indices/Salmon_SMEM/transcripts_index_SMEM --type fmd -p 8
            stop_Salmon_SMEM_index=`date +%s`
            printf $filename","$((stop_Salmon_SMEM_index-start_Salmon_SMEM_index)) >> Simulation/time_stats/time_Salmon_SMEM_index.csv
          fi

          #If the quasi index doesn't exist, make it
          if [ ! "$(ls -A Simulation/indices/Salmon_quasi)" ]; then
            start_Salmon_quasi_index=`date +%s`
            ./Simulation/Salmon-0.7.2_linux_x86_64/bin/salmon index -t Simulation/ref/reference.transcripts.fa -i Simulation/indices/Salmon_quasi/transcripts_index_quasi --type quasi -k 31 -p 8
            stop_Salmon_quasi_index=`date +%s`
            printf $filename","$((stop_Salmon_quasi_index-start_Salmon_quasi_index)) >> Simulation/time_stats/time_Salmon_quasi_index.csv
          fi
        fi

        if [ "$1" == "Kallisto" ]; then
            #If there is no Kallisto index, make it
            if [ ! "$(ls -A Simulation/indices/Kallisto)" ]; then
                start_kallisto_index=`date +%s`
                ./Simulation/kallisto_linux-v0.43.0/kallisto index -i Simulation/indices/Kallisto/transcripts.idx Simulation/ref/reference.transcripts.fa
                stop_kallisto_index=`date +%s`
                printf $((stop_kallisto_index-start_kallisto_index)) >> Simulation/time_stats/time_kallisto_index.csv
            fi
        fi

        if [ "$1" == "Sailfish" ]; then
              #If there is no index for sailfish, make it
              if [ ! "$(ls -A Simulation/indices/Sailfish)" ]; then
                export LD_LIBRARY_PATH=`pwd`/Simulation/Sailfish-0.6.3-Linux_x86-64/lib:$LD_LIBRARY_PATH
                export PATH=`pwd`/Simulation/Sailfish-0.6.3-Linux_x86-64/bin:$PATH
                start_sailfish_index=`date +%s`
                LC_ALL=C ./Simulation/Sailfish-0.6.3-Linux_x86-64/bin/sailfish index -p 8 -t Simulation/ref/reference.transcripts.fa -o Simulation/indices/Sailfish/ -k 31
                stop_sailfish_index=`date +%s`
                printf $filename","$((stop_sailfish_index-start_sailfish_index)) >> Simulation/time_stats/time_sailfish_index.csv
              fi
        fi

        if [ "$1" == "BRIE" ]; then
              #If there mouse_features.csv doesn't exist, make it
              if [ ! "$(ls -A BRIE_splicing_events/mouse_features.csv)" ]; then
                #activate venv
                source Simulation/venv/bin/activate

                brie-event -a $TEAM/genomes/Mus_musculus.GRCm38.82.chr.gtf -o Simulation/BRIE_splicing_events

                brie-event-filter -a Simulation/BRIE_splicing_events/SE.gff3 --anno_ref=$TEAM/genomes/Mus_musculus.GRCm38.82.chr.gtf --reference=$TEAM/genomes/Mus_musculus.GRCm38.dna.toplevel.fa

                chmod +x ./bigWigSummary
                memory=`pwd`
                export PATH="$memory:$PATH"

                brie-factor -a Simulation/BRIE_splicing_events/SE.gold.gtf -r $TEAM/genomes/Mus_musculus.GRCm38.dna.toplevel.fa -c mm10.60way.phastCons.bw -o Simulation/BRIE_splicing_events/mouse_features.csv -p 10
              fi
        fi

        #Run tool on simulated cells. Each cell is submitted as a seperate job.
        cd Simulation/data/simulated
        for i in $(find . -name '*.fastq*' -o -name '*.fq*');
        do
          base=`echo $i |awk -F/ '{print $2}'`
          filename=`echo $base |awk -F_ '{print $1}'`
          cd $memory
          #The line below will need to be edited for your LSF job system.
          bsub -n8 -R"span[hosts=1]" -c 99999 -G team_hemberg -q normal -o $TEAM/temp.logs/"output."$filename$1 -e $TEAM/temp.logs/"error."$filename$1 -R"select[mem>200000] rusage[mem=200000]" -M 200000 ./quantify.sh $1 $filename
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