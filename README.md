# Drop-seq_pipeline

Prerequisites:

-virtualenv

-github account

-reference genome gtf and fasta files. The Mus musculus Ensembl release 89 genome was used in this study

-java version 1.8

-R version 3.4.4

To run the pipeline:

Execute ./Dropseq_wrapper_script.sh path/to/java path/to/ref/fasta path/to/ref/gtf

In practice it is unlikely that your machine will have the resources to run the entire pipeline in one go, so you will probably need to split up the wrapper script and run it in bits.

As part of the pipeline, quality control steps are automatically carried out. For reference, these are the statistics used to filter the raw data:

| Statistic | Name of statistic in table | Threshold |
-------------|--------|---------
|No. uniquely mapping reads|Unique    | >20000 |
|No. of non-uniquely mapping reads|NonUnique|>3000|
|No. alignments|NumAlign|>25000|
|No. of reads|NumReads|>20000|

And the simulated data:

| Statistic | Name of statistic in table | Threshold |
-------------|--------|---------
|No. uniquely mapping reads|Unique    | <1000 |
|No. of non-uniquely mapping reads|NonUnique|>500|
|No. alignments|NumAlign|<1,000|

In addition, the scater package was used to filter cells in which more than 10% of reads mapped to mitochondrial genes in both the raw and simulated data.
