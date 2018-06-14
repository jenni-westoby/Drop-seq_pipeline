
#!/bin/bash

Kallisto() {

      #make a directory for the results of Kallisto for each cell
      filename=$1
      mkdir Simulation/Kallisto_real_results/$filename

      ./Simulation/kallisto_linux-v0.43.1/kallisto quant -i Simulation/indices/Kallisto/transcripts.idx --threads=8 --output-dir=Simulation/Kallisto_real_results/$filename Dropseq_Alignment_Cookbook/demultiplexed_fastqs/$1

      echo "target_id       length  eff_length      est_counts      tpm" >> Simulation/Kallisto_real_results/$filename/abundancesorted.tsv
      tail -n +2 Simulation/Kallisto_real_results/$filename/abundance.tsv | sort -n -k1.8 >> Simulation/Kallisto_real_results/$filename/abundancesorted.tsv
      mv Simulation/Kallisto_real_results/$filename/abundancesorted.tsv Simulation/Kallisto_real_results/$filename/abundance.tsv

}

"$@"
