#!/bin/bash
#Quantification pipeline
#Salmon results

RSEM(){

  filename=$1

  mkdir Simulation/RSEM_results

  #Start the clock for RSEM
  start_RSEM=`date +%s`

  #RSEM
  ./Simulation/RSEM-1.3.0/rsem-calculate-expression --star\
		    --star-path Simulation/STAR/bin/Linux_x86_64 \
		    -p 8 \
		    --fragment-length-mean 1600 --fragment-length-sd 50 \
        --append-names \
		    --single-cell-prior --calc-pme \
		    --time \
        Simulation/data/simulated/$filename'.fq' \
        Simulation/ref/reference Simulation/RSEM_results/$filename


	#Stop the clock for RSEM
	stop_RSEM=`date +%s`

	#Delete everything generated by RSEM except isoforms.results
	rm -r Simulation/RSEM_results/$filename'.stat'
	rm -r Simulation/RSEM_results/$filename'.transcript'*
	rm -r Simulation/RSEM_results/$filename'.genome'*
	rm -r Simulation/RSEM_results/$filename'.genes.results'

  #Trim the text added to the transcript names in the results file
  python ./trim.py `pwd` Simulation/RSEM_results/ $filename'.isoforms.results'
  mv Simulation/RSEM_results/'trimmed'$filename'.isoforms.results' Simulation/RSEM_results/$filename'.isoforms.results'

  #Sort the results file
  echo "transcript      gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct  posterior_mean_count    posterior_standard_deviation_of_count   pme_TPM pme_FPKM        IsoPct_from_pme_TPM" > Simulation/RSEM_results/$filename'.sortedisoforms.results'
  tail -n +2  Simulation/RSEM_results/$filename'.isoforms.results' | sort -n -k1.8 >> Simulation/RSEM_results/$filename'.sortedisoforms.results'
  mv Simulation/RSEM_results/$filename'.sortedisoforms.results' Simulation/RSEM_results/$filename'.isoforms.results'

  time_align=`grep "Aligning reads:" Simulation/RSEM_results/$filename".time" | awk '{print $3}'`
  time_expr=`grep "Estimating expression levels:" Simulation/RSEM_results/$filename".time" | awk '{print $4}'`
  time_cred=`grep "Calculating credibility intervals:" Simulation/RSEM_results/$filename".time" | awk '{print $4}'`

	printf $filename","$((stop_RSEM-start_RSEM))","$time_align","$time_expr","$time_cred"\n" >> Simulation/time_stats/time_RSEM.csv
}




Salmon(){

	#Rename/reformat input arguments
  	filename=$1

	mkdir Simulation/Salmon_results/Salmon_Alignment_Results/$filename
	mkdir Simulation/Salmon_results/Salmon_SMEM_results/$filename
	mkdir Simulation/Salmon_results/Salmon_quasi_results/$filename

  	if [ ! -f Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.out.bam' ]; then
    		STAR $filename
  	fi

	#Start the clock for Salmon alignment mode
	start_Salmon_align=`date +%s`

	#Salmon alignment mode
	Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon --no-version-check quant -a Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.out.bam' --fldMean 1600 --fldSD 50 -t Simulation/ref/reference.transcripts.fa -l A -o Simulation/Salmon_results/Salmon_Alignment_Results/$filename -p 8

	#Stop the clock for Salmon alignment mode

	stop_Salmon_align=`date +%s`

	rm -r Simulation/Salmon_results/Salmon_Alignment_Results/$filename/aux_info
	rm -r Simulation/Salmon_results/Salmon_Alignment_Results/$filename/cmd_info.json
	rm -r Simulation/Salmon_results/Salmon_Alignment_Results/$filename/lib_format_counts.json
	rm -r Simulation/Salmon_results/Salmon_Alignment_Results/$filename/libParams
	rm -r Simulation/Salmon_results/Salmon_Alignment_Results/$filename/logs

  	#Sort the results file
  	echo "Name    Length  EffectiveLength TPM     NumReads" > Simulation/Salmon_results/Salmon_Alignment_Results/$filename/quantsorted.sf
  	tail -n +2 Simulation/Salmon_results/Salmon_Alignment_Results/$filename/quant.sf | sort -n -k1.8 >> Simulation/Salmon_results/Salmon_Alignment_Results/$filename/quantsorted.sf
  	mv Simulation/Salmon_results/Salmon_Alignment_Results/$filename/quantsorted.sf Simulation/Salmon_results/Salmon_Alignment_Results/$filename/quant.sf

	#Start the clock for Salmon alignment free SMEM
	start_Salmon_SMEM=`date +%s`

	#Salmon alignment free SMEM
	Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon --no-version-check quant -i Simulation/indices/Salmon_SMEM/transcripts_index_SMEM -l A -r Simulation/data/simulated/$filename'.fq' --fldMean 1600 --fldSD 50 -o Simulation/Salmon_results/Salmon_SMEM_results/$filename -p 8

	#Stop the clock for Salmon SMEM
	stop_Salmon_SMEM=`date +%s`

	rm -r Simulation/Salmon_results/Salmon_SMEM_results/$filename/aux_info
	rm -r Simulation/Salmon_results/Salmon_SMEM_results/$filename/cmd_info.json
	rm -r Simulation/Salmon_results/Salmon_SMEM_results/$filename/lib_format_counts.json
	rm -r Simulation/Salmon_results/Salmon_SMEM_results/$filename/libParams
	rm -r Simulation/Salmon_results/Salmon_SMEM_results/$filename/logs

  	#Sort the results file
  	echo "Name    Length  EffectiveLength TPM     NumReads" > Simulation/Salmon_results/Salmon_SMEM_results/$filename/quantsorted.sf
  	tail -n +2 Simulation/Salmon_results/Salmon_SMEM_results/$filename/quant.sf | sort -n -k1.8 >> Simulation/Salmon_results/Salmon_SMEM_results/$filename/quantsorted.sf
  	mv Simulation/Salmon_results/Salmon_SMEM_results/$filename/quantsorted.sf Simulation/Salmon_results/Salmon_SMEM_results/$filename/quant.sf

	#Start the clock for alignment free quasi
	start_Salmon_quasi=`date +%s`

	#Salmon alignment free quasi
	Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon --no-version-check quant -i Simulation/indices/Salmon_quasi/transcripts_index_quasi -l A -r Simulation/data/simulated/$filename'.fq' --fldMean 1600 --fldSD 50 -o Simulation/Salmon_results/Salmon_quasi_results/$filename -p 8

	#Stop the clock for other alignment quasi
	stop_Salmon_quasi=`date +%s`

	rm -r Simulation/Salmon_results/Salmon_quasi_results/$filename/aux_info
	rm -r Simulation/Salmon_results/Salmon_quasi_results/$filename/cmd_info.json
	rm -r Simulation/Salmon_results/Salmon_quasi_results/$filename/lib_format_counts.json
	rm -r Simulation/Salmon_results/Salmon_quasi_results/$filename/libParams
	rm -r Simulation/Salmon_results/Salmon_quasi_results/$filename/logs

  	#Sort the results file
  	echo "Name    Length  EffectiveLength TPM     NumReads" > Simulation/Salmon_results/Salmon_quasi_results/$filename/quantsorted.sf
  	tail -n +2 Simulation/Salmon_results/Salmon_quasi_results/$filename/quant.sf | sort -n -k1.8 >> Simulation/Salmon_results/Salmon_quasi_results/$filename/quantsorted.sf
  	mv Simulation/Salmon_results/Salmon_quasi_results/$filename/quantsorted.sf Simulation/Salmon_results/Salmon_quasi_results/$filename/quant.sf

	#Output all times into Time.csv
	printf $filename","$((stop_Salmon_align-start_Salmon_align))","$((stop_Salmon_SMEM-start_Salmon_SMEM))","$((stop_Salmon_quasi-start_Salmon_quasi))"\n" >> Simulation/time_stats/time_Salmon.csv

}

eXpress () {

  	filename=$1

	#make a directory for the results of eXpress for each cell
	mkdir Simulation/eXpress_results/$filename

  	if [ ! -f Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.out.bam' ]; then
    		STAR $filename
  	fi

	#Start the clock for eXpress
	start_eXpress=`date +%s`

	#Run eXpress
	./Simulation/express-1.5.1-linux_x86_64/express Simulation/ref/reference.transcripts.fa Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.out.bam' -m 1600 -s 50 -o Simulation/eXpress_results/$filename

	#Stop the clock for eXpress
	stop_eXpress=`date +%s`

  	#Remove first column of results file then sort the file
  	cut -d$'\t' -f 2- Simulation/eXpress_results/$filename/results.xprs > Simulation/eXpress_results/$filename/resultscut.xprs
  	echo "target_id       length  eff_length      tot_counts      uniq_counts     est_counts      eff_counts      ambig_distr_alpha       ambig_distr_beta        fpkm    fpkm_conf_low   fpkm_conf_high  solvable        tpm" > Simulation/eXpress_results/$filename/results.xprs
  	tail -n +2 Simulation/eXpress_results/$filename/resultscut.xprs | sort -n -k1.8  >> Simulation/eXpress_results/$filename/results.xprs

	printf $filename","$((stop_eXpress-start_eXpress))"\n" >> Simulation/time_stats/time_eXpress.csv

}

Kallisto () {

	#make a directory for the results of Kallisto for each cell
	filename=$1
  mkdir Simulation/Kallisto_results
	mkdir Simulation/Kallisto_results/$filename

	#Start the clock for kallisto
	start_kallisto=`date +%s`

	./Simulation/kallisto_linux-v0.43.1/kallisto quant -i Simulation/indices/Kallisto/transcripts.idx --threads=8 --output-dir=Simulation/Kallisto_results/$filename --single -l 1600 -s 50 Simulation/data/simulated/$filename'.fq'

	#Stop the clock for kallisto
	stop_kallisto=`date +%s`

	printf $filename","$((stop_kallisto-start_kallisto))"\n" >> Simulation/time_stats/time_kallisto.csv

  	echo "target_id       length  eff_length      est_counts      tpm" >> Simulation/Kallisto_results/$filename/abundancesorted.tsv
  	tail -n +2 Simulation/Kallisto_results/$filename/abundance.tsv | sort -n -k1.8 >> Simulation/Kallisto_results/$filename/abundancesorted.tsv
 	mv Simulation/Kallisto_results/$filename/abundancesorted.tsv Simulation/Kallisto_results/$filename/abundance.tsv

}

Sailfish(){

  #make a directory for the results of sailfish for each cell
  filename=$1
  library_type=$2
  mkdir Simulation/Sailfish_results
  mkdir Simulation/Sailfish_results/$filename
  export LD_LIBRARY_PATH=`pwd`/SailfishBeta-0.10.0_CentOS5/lib:$LD_LIBRARY_PATH
  export PATH=`pwd`/SailfishBeta-0.10.0_CentOS5/bin:$PATH

  #Start the clock for sailfish
  start_sailfish=`date +%s`

  ./Simulation/SailfishBeta-0.10.0_CentOS5/bin/sailfish quant -p 8 -i Simulation/indices/Sailfish/ -l "IU"  --unmatedReads Simulation/data/simulated/$filename'.fq' -o Simulation/Sailfish_results/$filename
  stop_sailfish=`date +%s`

  printf $filename","$((stop_sailfish-start_sailfish))"\n" >> Simulation/time_stats/time_sailfish.csv

  rm Simulation/Sailfish_results/$filename/logs
  rm Simulation/Sailfish_results/$filename/quant_bias_corrected.sf
  rm Simulation/Sailfish_results/$filename/reads.count_info
  rm Simulation/Sailfish_results/$filename/reads.sfc

  echo "Name    Length  EffectiveLength TPM     NumReads" > Simulation/Sailfish_results/$filename/quantsorted.sf
  tail -n +2 Simulation/Sailfish_results/$filename/quant.sf | sort -n -k1.8 >> Simulation/Sailfish_results/$filename/quantsorted.sf
  mv Simulation/Sailfish_results/$filename/quantsorted.sf Simulation/Sailfish_results/$filename/quant.sf

}

STAR(){

  filename=$1

	#Start the clock for STAR
	start_STAR=`date +%s`

	#Make STAR reference
	./Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir Simulation/ref --readFilesIn Simulation/data/simulated/$filename'.fq' --outFileNamePrefix Simulation/bamfiles/simulated/$filename --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM

	#Stop the clock for STAR
	stop_STAR=`date +%s`

	#Start the clock for Samtools sort
	start_samtools_sort=`date +%s`

	#Sort STAR bamfiles
	./Simulation/samtools-1.5/samtools sort Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.out.bam' -T Simulation/bamfiles/simulated/$filename'temp' -o Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.sortedByCoord.out.bam'


	#Stop the clock for Samtools sort
	stop_samtools_sort=`date +%s`

	#Start the clock for Samtools index
	start_samtools_index=`date +%s`

	#Index STAR bamfiles
	./Simulation/samtools-1.5/samtools index Simulation/bamfiles/simulated/$filename'Aligned.toTranscriptome.sortedByCoord.out.bam'

	#Stop the clock for Samtools
	stop_samtools_index=`date +%s`

	printf $filename","$((stop_STAR-start_STAR))","$((stop_samtools_sort-start_samtools_sort))","$((stop_samtools_index-start_samtools_index))"\n" >> Simulation/time_stats/time_STAR_samtools.csv

}


#Export functions
export -f RSEM
export -f Salmon
export -f eXpress
export -f Kallisto
export -f STAR
export -f Sailfish

"$@"
