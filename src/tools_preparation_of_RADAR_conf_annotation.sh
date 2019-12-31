#--------------------------------------------------------------------------------------------------------------------------------------------------
#### Example of generating SNP annotation from NCBI dbSNP divided by chromosome 

## directory
dir_anno=~

## input: reference genome fasta file
genome_fasta=${dir_anno}/reference/Human/hg38/hg38_all.fa
## input: variant downloaded from NCBI dbSNP
dbSNP_all=${dir_anno}/annotation/Human/hg38/SNP/dbSNP_b151/NCBI_dbSNP_b151_all_hg38.vcf

## output: SNP from dbSNP
SNP_dbSNP=${dir_anno}/annotation/Human/hg38/SNP/NCBI_dbSNP_b151_all_hg38_SNP.vcf
## output: SNP annotation from NCBI dbSNP divided by chromosome
SNP_dbSNP_divided_by_chromosome=${dir_anno}/annotation/Human/hg38/SNP/dbSNP_b151/split_chr
test -d ${SNP_dbSNP_divided_by_chromosome} || mkdir -p ${SNP_dbSNP_divided_by_chromosome}

gatk SelectVariants -select-type SNP -R ${genome_fasta} -V  ${dbSNP_all} -O ${SNP_dbSNP}
awk -v dir_SNP=${SNP_dbSNP_divided_by_chromosome} 'BEGIN{OFS="\t"}$0 !~/^#/{$1=substr($1,4,length($1));print $0 >dir_SNP"/chr"$1}' ${SNP_dbSNP}
gzip ${SNP_dbSNP_divided_by_chromosome}/chr*


#--------------------------------------------------------------------------------------------------------------------------------------------------
#### Example of generating annotation_Alu, annotation_Repetitive_non_Alu, annotation_All_repetitive, annotation_simple_repeats

## input: annotation of RepeatMasker downloaded from UCSC (http://genome.ucsc.edu/)
UCSC_RepeatMask=UCSC_RepeatMask_hg38.txt

cat ${UCSC_RepeatMask} | awk '{print $6"\t"$7"\t"$8"\t"$11"\t"$2"\t"$10}' >  ${annotation_All_repetitive}
name_Alu="FAM_FRAM_FLAM_A_FLAM_C"
awk -v nA=${name_Alu} '{if(index(nA,$4)>0 || $4 ~ "Alu"){print $0}}'  ${annotation_All_repetitive} > ${annotation_Alu}
awk 'FILENAME==ARGV[1]{a[$0]++}FILENAME==ARGV[2]&& !a[$0]'  ${annotation_Alu} ${annotation_All_repetitive}  > ${annotation_Repetitive_non_Alu}
awk '{ if( $4 ~/\([ACTGN]?[ACTGN]?[ACTGN]?[ACTGN]?[ACTGN]?\)/ ){print} }' ${annotation_All_repetitive} > ${annotation_simple_repeats}
 

#### Example of generating annotation_intronic_4site, annotation_gene_transcribed_strands
## input: UCSC refFlat annotation 
UCSC_refFlat=UCSC_refFlat.txt

cat ${UCSC_refFlat} | awk '{split($10,exonStarts,","); for(idx=1;idx<length(exonStarts);idx++){print $3"\t"(exonStarts[idx]-4)"\t"(exonStarts[idx])"\t"$2"\t0\t"$4};split($11,exonEnds,",");for(idx=1;idx<length(exonEnds);idx++){print $3"\t"(exonEnds[idx])"\t"(exonEnds[idx]+4)"\t"$2"\t0\t"$4}; }' >  ${annotation_intronic_4site}
cat ${UCSC_refFlat} | awk '{print $3"\t"$5"\t"$6"\t"$2"\t""0\t"$4}' > ${annotation_gene_transcribed_strands}

