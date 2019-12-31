# RADAR
**RADAR** (**R**NA-editing **A**nalysis-pipeline to **D**ecode **A**ll twelve-types of **R**NA-editing events) is devised to detect and visualize all possible twelve-types of RNA editing events from RNA-seq datasets.
## Features
RADAR can be conveniently applied to identify RNA-editing from RNA-seq data with stringent filtering steps.
* All possible RNA-editing events from each given RNA-seq dataset are summarized into an Excel file.
* Numbers of all twelve-types of RNA editing events are plotted by histograms according to their genomic locations in Alu, repetitive non-Alu and non-repetitive regions.
* Manhattan plots are further used to illustrate RNA editing ratios of selected types of RNA-editing events, such as C-to-U or A-to-G.

## Schema
<img src="https://github.com/xiongyichun/RADAR/blob/master/docs/RADAR.jpg"  alt="RADAR pipeline" />

## Installation requirements
RADAR can be run directly without setup process after downloaded and unzipped, only if tools it depends on have been installed:

1. [HISAT2 (>=2.0.5)](https://ccb.jhu.edu/software/hisat2/index.shtml)
2. [BWA (>=0.7.9)](http://bio-bwa.sourceforge.net/)
3. [Samtools (>=1.7)](http://www.htslib.org/)
4. [GATK (>=4.0.1.0)](https://software.broadinstitute.org/gatk/)
5. [R (>=3.0.2)](https://www.r-project.org)<br/>
    Dependent R packages:
    * [ggplot2](https://ggplot2.tidyverse.org/index.html)
    * [dplyr](https://dplyr.tidyverse.org/index.html)
    
## Installation
As long as all **Installation requirements** have been fulfilled, RADAR can be run directly without setup process after downloaded and unzipped. 
```bash
git clone https://github.com/YangLab/RADAR
cd RADAR
chmod +x RADAR
./RADAR -h
```

## Configuration
Reference genome, genomic sequence index and genomic annotations should be provided to RADAR within the **RADAR.conf** file. Examples of how to get these annotations have been provided in **src/tools_preparation_of_RADAR_conf_annotation.sh**.

#### Reference and sequence index
1. Genome build version of the reference genome. <br />
    * Example: `genome_build_version=hg38`
2. Ribosomal DNA (rDNA) sequence index for BWA MEM, which can be created by command "bwa index \~/reference/Human/RNA_45S5/RNA45S5.fa" <br />
     * Example: `rDNA_idnex_bwa_mem=~/reference/Human/RNA_45S5/RNA45S5.fa`
3. Path to reference genome <br />
     * Example: `genome_fasta=~/reference/Human/hg38/hg38_all.fa`
4. Reference genome sequence index for HISAT2, which can be created by command "hisat2-build \~/reference/Human/hg38/hg38_all.fa \~/reference/Human/hg38/hg38_all.fa" <br />
     * Example: `genome_index_hisat2=~/reference/Human/hg38/hg38_all.fa`
5. Reference genome sequence index for BWA MEM, which can be created by command "bwa index \~/reference/Human/hg38/hg38_all.fa" <br />
     * Example: `genome_index_bwa_mem=~/reference/Human/hg38/hg38_all.fa`
6. Reference genome sequence index for Blat, which can be created by command "RADAR/tools/faToTwoBit \~/reference/Human/hg38/hg38_all.fa \~/reference/Human/hg38/hg38_all.fa.2bit" <br />
     * Example: `genome_index_blat=~/reference/Human/hg38/hg38_all.fa.2bit`
7. Reference genome sequence index for GATK in the same directory with reference genome, which can be created by command "gatk CreateSequenceDictionary -R \~/reference/Human/hg38/hg38_all.fa" <br />
     * Example: `genome_index_gatk=~/reference/Human/hg38/hg38_all.dict`

#### Variants annotation: dbSNP, 1000Genome, EVS
1. Total variants annotation from [NCBI dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) and GATK index in the same directory.<br />
     * Example of the total .vcf file: `dbSNP_all=~/annotation/Human/hg38/SNP/dbSNP_b151/NCBI_dbSNP_b151_all_hg38.vcf`<br />
     * Example of the GATK index for total .vcf (which can be created by command "gatk IndexFeatureFile -F \~/annotation/Human/hg38/SNP/dbSNP_b151/NCBI_dbSNP_b151_all_hg38.vcf"): `dbSNP_all_index_gatk=~/annotation/Human/hg38/SNP/dbSNP_b151/NCBI_dbSNP_b151_all_hg38.vcf.idx ` <br />
2. SNP annotation from NCBI dbSNP divided by chromosome  <br />
     * Example of the folder for NCBI dbSNP divided by chromosome: `SNP_dbSNP_divided_by_chromosome=~/annotation/Human/hg38/SNP/dbSNP_b151/split_chr`
3. SNP annotation from [The 1000 Genomes Project](https://www.internationalgenome.org/) divided by chromosome <br />
     * Example of the folder for SNP divided by chromosome: `SNP_1000Genome_divided_by_chromosome=~/annotation/Human/hg38/SNP/1000genomes/split_chr`
4. SNP annotation from [The University of Washington Exome Sequencing Project](http://evs.gs.washington.edu/EVS/) divided by chromosome <br />
     * Example of the folder for SNP divided by chromosome: `SNP_EVS_divided_by_chromosome=~/annotation/Human/hg38/SNP/EVS/split_chr`

#### Genome annotation
All genome annotations are in BED format:

| Field       | Description                                      |
| :---------- | :----------------------------------------------- |
| chrom       | Chromosome                                       |
| start       | Start position                                   |
| end         | End position                                     |
| name        | Repeat name or gene symbol/gene name             |
| score       | Smith Waterman alignment score for repeat region |
| strand      | + or - for strand                                |  

1. Annotation of Alu, repetitive non-Alu, non-repetitive genomic region in the BED format <br />
     * Example of Alu annotation: `annotation_Alu=~/annotation/Human/hg38/Alu.bed ` <br />
     * Example of repetitive non-Alu annotation: `annotation_Repetitive_non_Alu=~/annotation/Human/hg38/Repetitive_non-Alu.bed` <br />
     * Example of all repetitive annotation: `annotation_All_repetitive=~/annotation/Human/hg38/All_repetitive.bed` <br />
2. Annotation of RepeatMasker simple repeats from UCSC in BED format <br />
     * Example: `annotation_simple_repeats=~/annotation/Human/hg38/UCSC_RepeatMask_SimpleRepeats_hg38.bed`
3. Annotation of splice sites in BED format, which is used as the input of option "--known-splicesite-infile" during HISAT2 mapping. Official website of [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) has detailed how to create it.<br />
    * Example: `annotation_splice_sites=~/annotation/Human/hg38/ref_all_spsites_hg38.bed`
4. Annotation of intronic 4 bp flanking splice sites <br />
     * Example: `annotation_intronic_4site=~/annotation/Human/hg38/hg38_intronic_4site.bed`
5. Annotation of transcribed strands of genes <br />
     * Example: `annotation_gene_transcribed_strands=~/annotation/Human/hg38/ref_UCSC_refFlat.bed`
    
    


## Usage

RADAR pipeline can break down into three main steps, while read mapping and RNA-editing calling are integrated into one command during operation:

### STEP 1: Read mapping and RNA-editing calling
* For paired-end RNA-seq data: <br />
COMMAND: `./RADAR read_mapping_and_RNA_editing_calling -1 "full_path_of_fastq1" -2 "full_path_of_fastq2" --stranded "true/false"  -o "output_dir" -n "outname"  -t "maximum_threads" `
* For single-end RNA-seq data: <br />
COMMAND: `./RADAR read_mapping_and_RNA_editing_calling -s "full_path_of_fastq" --stranded "true/false"  -o "output_dir" -n "outname" -t "maximum_threads"  `
##### Options
`-s | --single | -single`: Fasta file for the single-end RNA-seq data. <br />
`-1 | --fq1 | -fq1`  and  `-2 | --fq2 | -fq2`: Fasta file for the paired-end RNA-seq data. <br />
`--stranded | -stranded`: If the RNA-seq is stranded or not. value: true or false. <br />
`-n | --outname | -outname`: The prefix of file name for the RNA-editing results. <br />
`-o | --outdir | -outdir`: Output directory of the results. <br />
`-t | --thread | -thread`: Maximum threads used for computation. <br />
`-h | --help | -help`: Print help information. <br />


### STEP 2: RNA-editing visualization
#### 1. Histogram plot for each treatment
COMMAND: `./RADAR histogram -i "outdir_of_read_mapping_and_RNA_editing_calling" -o "path_of_plot"  --outname_of_replicates "outname_of_replicates" `  <br />
##### Options
`-i | --inputdir | -inputdir`: The directory of the RNA-editing results.  <br />
`-o | --output | -output`: Full path of the pdf file for the histogram. <br />
`--outname_of_replicates | -outname_of_replicates`: The prefix of file name for the RNA-editing results for each replicates from the same treatment. The separator between outnames should be comma, for example, "s1_rep1,s1_rep2,s1_rep3". <br />
`-h | --help | -help`: Print help information.  <br />

#### 2. Manhattan plot of specific RNA-editing type 
COMMAND: `./RADAR Manhattan_plot -i "outdir_of_read_mapping_and_RNA_editing_calling" -o "path_of_plot" --RNA_editing_type "RNA_editing_type" --outname_of_samples "outname_of_samples_to_plot" --color_of_samples "colors_of_samples_in_the_plot" `  <br />
##### Options
`-i | --inputdir | -inputdir`: The directory of the RNA-editing results.  <br />
`-o | --output | -output`: Full path of the pdf file for the Manhattan plot.  <br />
`--RNA_editing_type | -RNA_editing_type`: Interested RNA-editing type for the Manhattan plot, which was selected from all twelve-types RNA-editing, including A-to-C, A-to-G, A-to-U, C-to-A, C-to-G, C-to-U, G-to-A, G-to-C, G-to-U, U-to-A, U-to-C, U-to-G. <br />
`--outname_of_samples | -outname_of_samples`: Outname of samples to plot. The separator between outnames should be comma, for example, "s1_rep1,s1_rep2,s1_rep3,s2_rep1,s2_rep2,s2_rep3". <br />
`--color_of_samples | -color_of_samples`: Color of hex RGB format for the dot of samples in the plot. Colors should be within double quotations and seperated by comma (,). For example, "#919191,#919191,#FF3F00,#FF3F00,#FF3F00". For multiple samples, provide matched colors and samples; for one sample, provide two colors to distinguish adjacent chromosomes. <br />
`-h | --help | -help`: Print help information.  <br />

## Citation
Wang X#, Ding C#, Yu W#, Xiong YC#, He S#, Yang B#, Wang Y, Li J, Lu Z, Zhu W, Wu J, Wei J, Huang X, Liu Z*, Yang L* and Chen J*. Cas12a base editors induce efficient and specific editing with low DNA damage response. in revision.

## License
Copyright (C) 2020 YangLab. [Licensed GPLv3](https://github.com/xiongyichun/RADAR/blob/master/LICENSE) for open source use or contact YangLab (yanglab@@picb.ac.cn) for commercial use.
