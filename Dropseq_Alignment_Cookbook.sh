#!/bin/bash

#script to carry out steps from Drop-seq alignment cookbook.

path_to_java=$1
path_to_ref_fasta=$2
path_to_ref_gtf=$3

cd Dropseq_Alignment_Cookbook

$path_to_java -jar picard.jar SortSam \
      I=SRR3587500.bam \
      O=SRR3587500_sorted.bam \
      SORT_ORDER=queryname

#Tag with cell barcodes
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar \
TagBamWithReadSequenceExtended INPUT=SRR3587500_sorted.bam \
OUTPUT=temp_files/SRR3587500_cell_barcode_tag.bam \
SUMMARY=temp_files/SRR3587500_cell_barcodes.txt \
BASE_RANGE=1-­12 BASE_QUALITY=10 BARCODED_READ=1 \
DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1

#Tag with UMIs
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar \
TagBamWithReadSequenceExtended INPUT=temp_files/SRR3587500_cell_barcode_tag.bam \
OUTPUT=temp_files/SRR3587500_unaligned_tagged_CellMolecular.bam \
SUMMARY=temp_files/unaligned_tagged_Molecular.bam_summary.txt \
BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True \
TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

#Remove reads where the cell barcode/UMI has low quality bases
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar \
FilterBAM TAG_REJECT=XQ \
INPUT=temp_files/SRR3587500_unaligned_tagged_CellMolecular.bam \
OUTPUT=temp_files/SRR3587500_unaligned_tagged_filtered.bam

#Trim SMART adapter bases
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar \
TrimStartingSequence INPUT=temp_files/SRR3587500_unaligned_tagged_filtered.bam \
OUTPUT=temp_files/SRR3587500_unaligned_tagged_trimmed_smart.bam \
OUTPUT_SUMMARY=temp_files/adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5

#Trim polyA
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar PolyATrimmer \
INPUT=temp_files/SRR3587500_unaligned_tagged_trimmed_smart.bam \
OUTPUT=temp_files/SRR3587500_unaligned_mc_tagged_polyA_filtered.bam \
OUTPUT_SUMMARY=$TEAM/Shekhar_retina_1_batch_1/polyA_trimming_report.txt \
MISMATCHES=0 NUM_BASES=6

#Convert bamfile to fastq
$path_to_java -XX:MaxHeapSize=1000m -jar picard.jar SamToFastq \
I=temp_files/SRR3587500_unaligned_mc_tagged_polyA_filtered.bam \
FASTQ=temp_files/SRR3587500_unaligned_mc_tagged_polyA_filtered.fastq \


#Generate STAR Index
../Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 \
--runMode genomeGenerate --genomeDir ../Simulation/indices/STAR \
--genomeFastaFiles $path_to_ref_fasta \
 --sjdbGTFfile $path_to_ref_gtf

#Run STAR
../Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 \
--genomeDir ../Simulation/indices/STAR \
--readFilesIn temp_files/SRR3587500_unaligned_mc_tagged_polyA_filtered.fastq \
--outFileNamePrefix temp_files/SRR3587500_star_ind

#Sort the sam file
$path_to_java -XX:MaxHeapSize=1000m -jar picard.jar SortSam \
I=temp_files/SRR3587500_star_indAligned.out.sam \
O=temp_files/SRR3587500_indaligned.sorted.bam \
SO=queryname

#Merge aligned + unaligned bams to recover cell barcode + UMI tags
$path_to_java -XX:MaxHeapSize=1000m -jar picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=$path_to_ref_fasta \
UNMAPPED_BAM=temp_files/SRR3587500_unaligned_mc_tagged_polyA_filtered.bam \
ALIGNED_BAM=temp_files/SRR3587500_indaligned.sorted.bam \
OUTPUT=temp_files/SRR3587500_ind_merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

#Tag reads with gene exons
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar \
TagReadWithGeneExon I=temp_files/SRR3587500_ind_merged.bam \
O=temp_files/SRR3587500_ind_merged_GE_tagged.bam \
ANNOTATIONS_FILE=$path_to_ref_gtf TAG=GE

#Bead synthesis error correction
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar DetectBeadSynthesisErrors \
I=temp_files/SRR3587500_ind_merged_GE_tagged.bam \
O=temp_files/SRR3587500_BSE.bam \
OUTPUT_STATS=temp_files/my_synthesis_stats_ind_BSE.txt \
SUMMARY=temp_files/my_synthesis_stats_ind_BSE.summary.txt \
NUM_BARCODES=20000 \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC


#Generate a digital gene expression matrix
$path_to_java -Xmx5000m -jar Drop-seq_tools-1.13/jar/dropseq.jar \
DigitalExpression I=temp_files/SRR3587500_BSE.bam \
O=temp_files/SRR3587500_BSE.bam.dge.txt.gz \
SUMMARY=temp_files/SRR3587500_BSE.bam.dge.summary.txt \
NUM_CORE_BARCODES=6000

#Put cell barcodes in a text file
gunzip temp_files/SRR3587500_BSE.bam.dge.txt.gz
head -n 1 temp_files/SRR3587500_BSE.bam.dge.txt | awk 'BEGIN{RS="\t"} 1' | tail -n +2 > \
temp_files/cell_barcodes.txt


# seeding adopted from https://stackoverflow.com/a/41962458/7820599
get_seeded_random()
{
  seed="$1";
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null;
}

#If you change the seed, watch out for blank lines in the output
seed=1;

#Randomly select 1000 cells
shuf --random-source=<(get_seeded_random $seed) -n 1000 temp_files/cell_barcodes.txt > temp_files/1000_cell_barcodes.txt

mkdir demultiplexed_bams

./demultiplex.sh temp_files/1000_cell_barcodes.txt

mkdir demultiplexed_fastqs

#Convert bams to fastqs
./convert_bam_to_fastq.sh Dropseq_Alignment_Cookbook/demultiplexed_bams /software/java/bin/java Dropseq_Alignment_Cookbook/picard.jar Dropseq_Alignment_Cookbook/demultiplexed_fastqs

#delete cell which triggered bugs in Salmon and RSEM
rm demultiplexed_fastqs/CTCGTTTTTTAA.fastq
