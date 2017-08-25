# Drop-seq_pipeline

Prerequisites:

-virtualenv

-github account

-reference genome gtf and fasta files. Note that if your data contains ERRCC spike-ins, you should concatenate the reference genome gtf file with the ERRCC gtf file, and concatenate the reference and ERRCC fastq files (see https://tools.thermofisher.com/content/sfs/manuals/cms_095048.txt)

-directory containing single cell RNA-seq data. This data should be demultiplexed, have any adaptors trimmed and should be in the format of gzipped fastq files.

To run the pipeline:

1. Follow the instructions from the Drop-seq Alignment Cookbook (http://mccarrolllab.com/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf) to generate a Digital Gene Expression (DGE) matrix. Make sure you keep the bamfile you used to generate the DGE matrix. Put the cell barcodes (the column names in the DGE matrix) into a text file with one barcode per line.

2. Use the following command to demultiplex the data:
  ```
  while read i;
  do
    samtools view path/to/bamfile/used/to/make/DGE | grep "XC:Z:"$i | cat path/to/header.sam - | samtools view -Sb - > path/to/output/directory/$i".bam"
  done < path/to/cell/barcodes/text/file
  ```
  Replace the paths with the correct paths for your data. If you are demultiplexing over 1000 cells, consider directing the   output into multiple directories. You can speed up the demultiplexing process by parallelising it. You will need a header for   your bamfiles - for convenience one is included in this directory (header.sam). The output of this command is one aligned   bamfile per cell.

4. Execute ./setup.sh setup. This will create a new directory called Simulation into which all the software required for this pipeline will be locally installed. In addition, empty directories are created within the Simulation directory which will eventually contain the RSEM references, various indices, the raw and simulated data, results matrices and graphs. This step will take ~30 minutes - 1 hour depending on your network speed.

5. Execute ./RSEM_ref.sh make_ref /path/to/gtf path/to/fasta, where the gtf and fasta files are the reference genome. This builds the RSEM reference.

6. Execute ./quality_control.sh QC path/to/gtf path/to/fasta path/to/raw/data. This creates a table of quality control statistics. Based on the results of this you can decide which cells you would like to simulate and which you are going to discard.

7. Once you have decided which cells to discard and have a directory containing only the gzipped cells you want to simulate, execute ./simulate.sh run_simulations path/to/raw/data. The simulated cells and their ground truth expression values are saved in Simulation/data/simulated.

8. If you wish, you can also perform quality control on your simulated cells based on read and alignment quality. This is probably wise, as RSEM sometimes generates cells with very few reads. Execute ./quality_control.sh QC path/to/gtf path/to/fasta and delete any problematic cells from the data/simulated directory.

9. Perform any further quality control you would like to perform prior to doing your benchmarking. For example, I use the scater package to filter based on patterns of expression, such as unusually high percentages of ERCCs. For this analysis I use the ground truth expression values produced in the simulation.

10. Before performing the quantification step, delete any cells that you don’t want to include in the benchmarking.

11. Execute ./benchmark.sh benchmark name_of_program_you_want_to_test. This will generate results matrices of expression values for the method you are interested in. Repeat for each method you want to test.

12. Execute ./make_matrix.sh make_matrix. This generates a compact results matrix for each method in results_matrices.
