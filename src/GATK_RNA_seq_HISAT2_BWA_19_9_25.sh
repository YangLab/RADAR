#!/usr/bin/env bash
set -eo pipefail
echo -e `date` dir:$PWD  command:$0 "$*" >>$0\.arg
######!!!!!!!!!Warining!!!!!!!!!!############
##USE bash instead of sh to run this shell script!
######!!!!!!!!!Warining!!!!!!!!!!############

################################
MYDIR=`dirname $0`
dir_scripts=${MYDIR}/../src
dir_tools=${MYDIR}/../tools


RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller(){
    #RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller $workspace $fq_path $outname $layout $genome_build_version

    ###### Wed Sep 25 20:00:00 CST 2019
    ####################init#################################################
    ##Usage:
    #RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller $outdir $fq_path $outname $layout $genome_build_version

    #Parameters example:
    #outdir="${ABE_m_path}/SameBam_TwoCaller/HISAT2_6mis_bam"
    #fq_path="${ABE_m_path}/dep_files/fastq"
    #outname="NT-Cas9_rep1"
    #layout="paired"
    #genome_build_version=hg38
    ##
    #used for RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller,\
    #assemble subfunctions:mapping, sam fine tune, caliing and variants filtering.

    ##required subfunctions:
    #HISAT2_2mismatch_following_BWA_6mismatch_mapping $outname $layout $fq_path $HISAT2_mapping_path $BWA_mapping_path $fine_tune_path
    #sam_fine_tune $bam_file $fine_tune_path
    #gatk_Variant_filter $bam_file "_2pass" $outdir ${outname}

    ###update:
    ###### Sat Nov 9 14:34:48 CST 2019
    #add parameter "$genome_build_version"
    ####################init END#################################################
	fq0=""
	fq1=""
	fq2=""
	stranded=false
	thread=10
	outdir="./"
	outname=""
	genome_build_version=""
	genome_index_hisat2=""
	genome_index_bwa_mem=""
	genome_index_blat=""
	genome_fasta=""
	SNP_dbSNP_divided_by_chromosome=""
	SNP_1000Genome_divided_by_chromosome=""
	SNP_EVS_divided_by_chromosome=""
	rDNA_index_bwa_mem=""
	dbSNP_all=""
	annotation_Alu=""
	annotation_Repetitive_non_Alu=""
	annotation_All_repetitive=""
	annotation_RepeatMasker_simple_repeats=""
	annotation_splice_sites=""
	annotation_intronic_4site=""
	annotation_gene_transcribed_strands=""
	
	GETOPT_ARGS=`getopt  -o s:1:2:t:o:n:g:  -al single:,fq1:,fq2:,stranded:,thread:,outdir:,outname:,genome_build_version:,genome_index_hisat2:,genome_index_bwa_mem:,genome_index_blat:,genome_fasta:,SNP_dbSNP_divided_by_chromosome:,SNP_1000Genome_divided_by_chromosome:,SNP_EVS_divided_by_chromosome:,rDNA_index_bwa_mem:,dbSNP_all:,annotation_Alu:,annotation_Repetitive_non_Alu:,annotation_All_repetitive:,annotation_RepeatMasker_simple_repeats:,annotation_intronic_4site:,annotation_splice_sites:,annotation_gene_transcribed_strands: -- "$@"`
echo ${GETOPT_ARGS}
	eval set -- "$GETOPT_ARGS"
	while [ -n "$1" ]
	do
        	case "$1" in
                	-s|--single | -single) fq0=$2; shift 2;;
                	-1|--fq1 | -fq1) fq1=$2; shift 2;;
                	-2|--fq2 | -fq2) fq2=$2; shift 2;;
			--stranded | -stranded) stranded=$2; shift 2;;
                	-t|--thread | -thread) thread=$2; shift 2;;
                	-o|--outdir | -outdir) outdir=$2; shift 2;;
               		-n|--outname | -outname) outname=$2; shift 2;;
                	-g|--genome_build_version | -genome_build_version) genome_build_version=$2; shift 2;;
			--rDNA_index_bwa_mem | -rDNA_index_bwa_mem) rDNA_index_bwa_mem=$2; shift 2;;
                	--genome_index_hisat2 | -genome_index_hisat2) genome_index_hisat2=$2; shift 2;;
                	--genome_index_bwa_mem | -genome_index_bwa_mem) genome_index_bwa_mem=$2; shift 2;;
			--genome_index_blat | -genome_index_blat) genome_index_blat=$2; shift 2;;
                	--genome_fasta | -genome_fasta) genome_fasta=$2; shift 2;;
			--dbSNP_all | -dbSNP_all) dbSNP_all=$2; shift 2;;
                	--SNP_dbSNP_divided_by_chromosome | -SNP_dbSNP_divided_by_chromosome) SNP_dbSNP_divided_by_chromosome=$2; shift 2;;
                	--SNP_1000Genome_divided_by_chromosome | -SNP_1000Genome_divided_by_chromosome) SNP_1000Genome_divided_by_chromosome=$2; shift 2;;
                	--SNP_EVS_divided_by_chromosome | -SNP_EVS_divided_by_chromosome) SNP_EVS_divided_by_chromosome=$2; shift 2;;
			--annotation_Alu | -annotation_Alu) annotation_Alu=$2; shift 2;;
			--annotation_Repetitive_non_Alu |annotation_Repetitive_non_Alu) annotation_Repetitive_non_Alu=$2; shift 2;;
			--annotation_All_repetitive | -annotation_All_repetitive) annotation_All_repetitive=$2; shift 2;;
			--annotation_RepeatMasker_simple_repeats | -annotation_RepeatMasker_simple_repeats) annotation_RepeatMasker_simple_repeats=$2; shift 2;;
			--annotation_intronic_4site | -annotation_intronic_4site) annotation_intronic_4site=$2; shift 2;;
                        --annotation_splice_sites | -annotation_splice_sites) annotation_splice_sites=$2; shift 2;;
			--annotation_gene_transcribed_strands | -annotation_gene_transcribed_strands) annotation_gene_transcribed_strands=$2; shift 2;;
                	--) break ;;
                	*) echo $1,$2,$show_usage; break ;;
        	esac
	done	




    tmp_a="bb${fq0}bb"
    test "$tmp_a" == "bbbb" && {
        layout=paired
        fq_path=`dirname ${fq1}`
        }|| {
            layout=single
            fq_path=`dirname ${fq0}`
        }


    #genome_build_version=$5  ###### Sat Nov 9 14:35:52 CST 2019 can serve as global value when this function is called.


    rRNAdeplete_path="${outdir}/${outname}_tmp/0-remove_rRNA"
    HISAT2_mapping_path="${outdir}/${outname}_tmp/1-1-0HISAT_map"
    BWA_mapping_path="${outdir}/${outname}_tmp/2-0-0BWA_map"



    fine_tune_path="${outdir}/${outname}_tmp/3-0-0Combine_bam"
    ###1. mapping, HISAT2_2mismatch_following_BWA_6mismatch_mapping
    HISAT2_2mismatch_following_BWA_6mismatch_mapping $outname $layout $fq_path $HISAT2_mapping_path $BWA_mapping_path $fine_tune_path ${thread} ${rRNAdeplete_path} ${rDNA_index_bwa_mem} ${annotation_splice_sites} 
    
    ###2. sam_fine_tune, Markduplicate, BQSR
    bam_file=${fine_tune_path}/${outname}_combine.bam
    sam_fine_tune $bam_file $fine_tune_path ${dbSNP_all}

    ###3. GATK_HaplotypeCaller Calling
    bam_file_prefix=${outname}_combine
    vcf_ready_bam_file=$fine_tune_path/${bam_file_prefix}_rgadd_dedupped_split_recal.bam
    gatk_Variant_filter $vcf_ready_bam_file "_2pass" $outdir ${outname} $genome_build_version ${genome_fasta} ${genome_index_blat} ${SNP_dbSNP_divided_by_chromosome} ${SNP_1000Genome_divided_by_chromosome} ${SNP_EVS_divided_by_chromosome} ${annotation_Alu} ${annotation_Repetitive_non_Alu} ${annotation_All_repetitive} ${annotation_RepeatMasker_simple_repeats}  ${annotation_intronic_4site} ${annotation_gene_transcribed_strands} ${stranded}
    
    ###4. HPB, editing ratio and clear temp file
    HPB_editing_ratio_clear ${bam_file}  $outdir ${outname} $genome_build_version 

}


HISAT2_2mismatch_following_BWA_6mismatch_mapping(){
    #HISAT2_2mismatch_following_BWA_6mismatch_mapping $outname $layout $fq_path $HISAT2_mapping_path $BWA_mapping_path $combine_bam
    ###### Wed Sep 25 19:04:01 CST 2019
    #used for HISAT2 and BWA two pass mapping.
    #First, HISAT2 2 mismatches
    #HISAT2 unmapped reads are extracted and send to BWA 6 mismatches mapping.
    #for paired ends reads, mismatches are all count at read level not fragment level, for example, HISAT2 2 mismatches means read1's mismatches maximum count is 2,read2's mismatches maximum count is 2, read1 mismatches plus read2 mismatches can up to 4.

    outname=$1
    layout=$2
    fq_path=$3
    HISAT2_mapping_path=$4
    BWA_mapping_path=$5
    combine_bam=$6
    thread=$7
    rRNAdeplete_path=$8   
    rDNA_index_bwa_mem=$9
    annotation_splice_sites=${10}
 
    test -d ${rRNAdeplete_path} || mkdir -p ${rRNAdeplete_path}
    test -d $HISAT2_mapping_path || mkdir -p $HISAT2_mapping_path
    test -d $BWA_mapping_path || mkdir -p $BWA_mapping_path
    test -d ${combine_bam} ||mkdir -p ${combine_bam}
    HISAT_map=$HISAT2_mapping_path
    bwa_map=${BWA_mapping_path}

    #ref_rRNA=reference/Human/RNA_45S5/RNA45S5.fa

    if [ "$layout" == "paired" ];then
	## remove rRNA
	bwa  mem -t ${thread} ${rDNA_index_bwa_mem}  $fq1  $fq2 > ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam 2>${rRNAdeplete_path}/log_BWA_rRNA_`date +%Y_%m_%d`.log
	samtools view -bh -f 4 ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam > ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
	samtools sort -n ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
	samtools fastq -1 ${rRNAdeplete_path}/${outname}_R1.fastq.gz -2 ${rRNAdeplete_path}/${outname}_R2.fastq.gz -s ${rRNAdeplete_path}/${outname}_singleton.fq ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
	rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
	rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
	fq1=${rRNAdeplete_path}/${outname}_R1.fastq.gz
	fq2=${rRNAdeplete_path}/${outname}_R2.fastq.gz

        ### 1. HISAT2 2 mismatches mapping
        ##################need confirm: --rna-strandness RF update:###### Wed Sep 25 19:56:57 CST 2019 confirmed, 0.12 VS 0.88
        #hisat2  --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${dep_path}/${genome_build_version}/${genome_build_version}_annotation/ref_all_spsites.txt --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p 10 -x /picb/rnomics1/database/Human/${genome_build_version}/genome/${genome_build_version}_all -1  ${fq_path}/${outname}${fq_suffix_1}  -2 ${fq_path}/${outname}${fq_suffix_2}  --un-conc-gz ${HISAT_map}/${outname}_un_conc_%.fastq.gz -S ${HISAT_map}/${outname}_HISAT2_mapped.sam  2>${HISAT_map}/log_HISAT2_2mismatch_${outname}_`date +%Y_%m_%d`.log 
        hisat2  --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${thread} -x   ${genome_index_hisat2} -1  $fq1  -2 $fq2  --un-conc-gz ${HISAT_map}/${outname}_un_conc_%.fastq.gz -S ${HISAT_map}/${outname}_HISAT2_mapped.sam  2>${HISAT_map}/log_HISAT2_2mismatch_${outname}_`date +%Y_%m_%d`.log ###### Wed Nov 20 08:06:39 CST 2019 fzc
        
        samtools view -h -F 4 ${HISAT_map}/${outname}_HISAT2_mapped.sam|awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) print $0 } else print $0 " not have XM tag"}}'|awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag"  }}' >${HISAT_map}/${outname}_unique_mismatch2.sam &

        #samtools view -f 4 -S ${HISAT_map}/${outname}_HISAT2_mapped.sam |awk 'BEGIN{FS="\t"}{print $1}'|sort|uniq >${HISAT_map}/${outname}_hisat2_unmap.readid 
        #$seqtk subseq ${HISAT_map}/${outname}_un_conc_1.fastq.gz ${HISAT_map}/${outname}_hisat2_unmap.readid |gzip > ${HISAT_map}/${outname}_unmapped_1.fastq.gz  &
        #$seqtk subseq ${HISAT_map}/${outname}_un_conc_2.fastq.gz ${HISAT_map}/${outname}_hisat2_unmap.readid |gzip > ${HISAT_map}/${outname}_unmapped_2.fastq.gz  
        #wait
        samtools view -bh -f 4 ${HISAT_map}/${outname}_HISAT2_mapped.sam > ${HISAT_map}/${outname}_HISAT2_mapped-unmapped.bam
        samtools sort -n ${HISAT_map}/${outname}_HISAT2_mapped-unmapped.bam  -o ${HISAT_map}/${outname}_HISAT2_mapped-unmapped_sorted.bam
        samtools fastq -1 ${HISAT_map}/${outname}_unmapped_1.fastq.gz -2 ${HISAT_map}/${outname}_unmapped_2.fastq.gz -s ${HISAT_map}/${outname}_unmapped_singleton.fq  ${HISAT_map}/${outname}_HISAT2_mapped-unmapped_sorted.bam


        ### 2. BWA 6 mismatches mapping
        
        test -d $bwa_map||mkdir -p $bwa_map
        bwa mem -t ${thread}  -A 1 -B 4  ${genome_index_bwa_mem}  ${HISAT_map}/${outname}_unmapped_1.fastq.gz  ${HISAT_map}/${outname}_unmapped_2.fastq.gz > ${bwa_map}/${outname}_bwa_mapped.sam 2>${bwa_map}/log_BWA_6mismatch_`date +%Y_%m_%d`.log


    elif [ "$layout" == "single" ];then
	## remove rRNA

	bwa  mem -t ${thread} ${rDNA_index_bwa_mem}  $fq0  > ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam 2>${rRNAdeplete_path}/log_BWA_rRNA_`date +%Y_%m_%d`.log
        samtools view -bh -f 4 ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam > ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
        samtools sort -n ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
        samtools fastq  ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam | gzip > ${rRNAdeplete_path}/${outname}.fastq.gz
        rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
        rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
        fq0=${rRNAdeplete_path}/${outname}.fastq.gz
	

	## two round mapping
        #hisat2 --secondary --no-temp-splicesite --known-splicesite-infile ${dep_path}/${genome_build_version}/${genome_build_version}_annotation/ref_all_spsites.txt --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p 10 -x /picb/rnomics1/database/Human/${genome_build_version}/genome/${genome_build_version}_all -U ${fq_path}/${outname}.fastq.gz -S ${HISAT_map}/${outname}_HISAT2_mapped.sam 2>${log_path}/${genome_build_version}/bmc/HISAT2/log_hisat2_${outname}.log 

        hisat2 --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${thread} -x ${genome_index_hisat2} -U $fq0 -S ${HISAT_map}/${outname}_HISAT2_mapped.sam 2>${HISAT_map}/log_HISAT2_2mismatch_${outname}_`date +%Y_%m_%d`.log 

        samtools view -h -F 4 ${HISAT_map}/${outname}_HISAT2_mapped.sam|awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) print $0 } else print $0 " not have XM tag" }}'|awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag" }}' >${HISAT_map}/${outname}_unique_mismatch2.sam &
        samtools view -bS -f 4 -o ${HISAT_map}/${outname}_HISAT2_unmapped.bam ${HISAT_map}/${outname}_HISAT2_mapped.sam
	samtools sort -n ${HISAT_map}/${outname}_HISAT2_unmapped.bam  -o ${HISAT_map}/${outname}_HISAT2_unmapped_sort.bam
        samtools fastq  ${HISAT_map}/${outname}_HISAT2_unmapped_sort.bam | gzip  > ${HISAT_map}/${outname}_HISAT2_unmapped.fastq.gz
 
        #/picb/rnomics4/rotation/fuzhican/software/conda/envs/circ/bin/bamToFastq -i  ${HISAT_map}/${outname}_HISAT2_unmapped.bam  -fq /dev/stdout | gzip >${HISAT_map}/${outname}_HISAT2_unmapped.fastq.gz

        ### 2. BWA 6 mismatches mapping
        bwa mem -t ${thread} ${genome_index_bwa_mem}  ${HISAT_map}/${outname}_HISAT2_unmapped.fastq.gz > ${bwa_map}/${outname}_bwa_mapped.sam
    fi
    python ${dir_scripts}/bwa_unique_mismatch6.py ${bwa_map}/${outname}_bwa_mapped.sam ${bwa_map}/${outname}_bwa_unique_mis6_mapq0.sam 
    
    samtools  view -bT ${genome_fasta} -o ${bwa_map}/${outname}_unmapped.nbam ${bwa_map}/${outname}_bwa_unique_mis6_mapq0.sam
    samtools  sort ${bwa_map}/${outname}_unmapped.nbam -o ${bwa_map}/${outname}_unmapped.sort.bam
    samtools  view -H ${bwa_map}/${outname}_unmapped.sort.bam > ${bwa_map}/${outname}_bwa.header
    
    wait
    cat ${bwa_map}/${outname}_bwa.header ${HISAT_map}/${outname}_unique_mismatch2.sam > ${combine_bam}/${outname}_accepted_hits.nsam
    samtools view -bT ${genome_fasta} -o  ${combine_bam}/${outname}_accepted_hits.nbam ${combine_bam}/${outname}_accepted_hits.nsam
    samtools sort ${combine_bam}/${outname}_accepted_hits.nbam -o ${combine_bam}/${outname}_accepted_hits.sort.bam

    samtools merge -f ${combine_bam}/${outname}_combine.bam ${combine_bam}/${outname}_accepted_hits.sort.bam ${bwa_map}/${outname}_unmapped.sort.bam
    samtools index ${combine_bam}/${outname}_combine.bam
#    need_rm_intermediate_file_list=(${combine_bam}/${outname}_accepted_hits.nsam ${combine_bam}/${outname}_accepted_hits.nbam ${combine_bam}/${outname}_accepted_hits.sort.bam ${bwa_map}/${outname}_unmapped.nbam ${bwa_map}/${outname}_unmapped.nbam ${bwa_map}/${outname}_unmapped.sort.bam ${bwa_map}/${outname}_bwa_unique_mis6_mapq0.sam ${bwa_map}/${outname}_bwa_mapped.sam ${HISAT_map}/${outname}_HISAT2_mapped.sam ${HISAT_map}/${outname}_unique_mismatch2.sam ${HISAT_map}/${outname}_unmapped_1.fastq.gz ${HISAT_map}/${outname}_unmapped_2.fastq.gz ${HISAT_map}/${outname}_hisat2_unmap.readid ${HISAT_map}/${outname}_un_conc_2.fastq.gz ${HISAT_map}/${outname}_un_conc_1.fastq.gz)
#    for need_rm_intermediate_file in ${need_rm_intermediate_file_list[@]}
#    do
#        test -e $need_rm_intermediate_file &&rm $need_rm_intermediate_file
#    done
}


sam_fine_tune(){
    #sam_fine_tune $bam_file $fine_tune_path $bam_file_basename_prefix
    bam_file=$1
    fine_tune_path=$2
    dbSNP_all=$3
    bam_file_basename_prefix=""  ###### Tue Oct 8 19:27:03 CST 2019 fzc added
    test -d $fine_tune_path || mkdir -p $fine_tune_path
    if [ "$bam_file_basename_prefix" == "" ];then
        bam_file_basename_prefix=`basename $bam_file|awk -F".bam" '{print $1}'`
    fi
    bam_file_prefix=${fine_tune_path}/${bam_file_basename_prefix}


    picard=${dir_tools}/picard.jar
    knownSNP_for_BQSR=${dbSNP_all}
	#${dep_path}/${genome_build_version}/SNP/dbSNP_b151/NCBI_dbSNP_b151_all_${genome_build_version}.vcf

    java -jar ${picard} AddOrReplaceReadGroups I=${bam_file} O=${bam_file_prefix}_rgadd.bam SO=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 

    java -jar $picard MarkDuplicates I=${bam_file_prefix}_rgadd.bam O=${bam_file_prefix}_rgadd_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${bam_file_prefix}_rgadd_MarkDuplicates_output.metrics 

    gatk SplitNCigarReads -R ${genome_fasta} -I ${bam_file_prefix}_rgadd_dedupped.bam -O ${bam_file_prefix}_rgadd_dedupped_split.bam

    gatk  BaseRecalibrator -R ${genome_fasta} -I ${bam_file_prefix}_rgadd_dedupped_split.bam --known-sites ${knownSNP_for_BQSR}  -O ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk4.grv 
    gatk ApplyBQSR  -R ${genome_fasta} -I ${bam_file_prefix}_rgadd_dedupped_split.bam  --bqsr-recal-file ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk4.grv -O ${bam_file_prefix}_rgadd_dedupped_split_recal.bam

    tmp_remove_file_list=(${bam_file_prefix}_rgadd.bam ${bam_file_prefix}_rgadd_dedupped.bam ${bam_file_prefix}_rgadd_dedupped_split.bam ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk4.grv ${bam_file_prefix}_rgadd_MarkDuplicates_output.metrics)

    for file1 in ${tmp_remove_file_list[@]}
    do
    test -e $file1 && rm $file1
    done
}


gatk_Variant_filter(){
    #gatk_Variant_filter $bam_file $tag $outdir
    #tag="_2pass"

    local from_bam=$1
    local tag=$2
    local outdir=$3
    local outname=$4
    local genome_build_version=$5
    local genome_fasta=$6
    local genome_index_blat=$7
    local SNP_dbSNP_divided_by_chromosome=$8
    local SNP_1000Genome_divided_by_chromosome=$9
    local SNP_EVS_divided_by_chromosome=${10}
    local annotation_Alu=${11} 
    local annotation_Repetitive_non_Alu=${12} 
    local annotation_All_repetitive=${13} 
    local annotation_RepeatMasker_simple_repeats=${14} 
    local annotation_intronic_4site=${15} 
    local annotation_gene_transcribed_strands=${16}   
    local stranded=${17}
 
    m1_path=${outdir}/${outname}_tmp/gatk
    local var_call=${m1_path}/4-0-0var_calling
    test -d ${var_call}||mkdir -p ${var_call}

    #--base-quality-score-threshold default 18
    gatk HaplotypeCaller -R ${genome_fasta} -I ${from_bam} --minimum-mapping-quality 0 --dont-use-soft-clipped-bases true  -stand-call-conf 0 -O ${var_call}/${outname}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0.vcf ###### Sun Nov 10 16:37:06 CST 2019 remove "--min-base-quality-score 20"

    echo -e "[`date`]BEGIN filter..."
    local vcf_filter_path=${m1_path}/5-0-0vcf_filter/${outname}
    test -d ${vcf_filter_path}||mkdir -p ${vcf_filter_path}
    #extract SNP
    gatk SelectVariants -select-type SNP -R ${genome_fasta} -V ${var_call}/${outname}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0.vcf -O ${vcf_filter_path}/${outname}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP.vcf
    inter_name="_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP"
    
    bash ${dir_scripts}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m gatk -n $outname -o ${outdir} -i ${inter_name} -g ${genome_build_version} --SNP_dbSNP_divided_by_chromosome ${SNP_dbSNP_divided_by_chromosome} --SNP_1000Genome_divided_by_chromosome ${SNP_1000Genome_divided_by_chromosome} --SNP_EVS_divided_by_chromosome ${SNP_EVS_divided_by_chromosome}
    bash ${dir_scripts}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m gatk -n $outname -o ${outdir} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${genome_build_version} --genome_fasta ${genome_fasta} --genome_index_blat ${genome_index_blat}  --annotation_Alu ${annotation_Alu} --annotation_Repetitive_non_Alu ${annotation_Repetitive_non_Alu} --annotation_All_repetitive ${annotation_All_repetitive} --annotation_RepeatMasker_simple_repeats ${annotation_RepeatMasker_simple_repeats}  --annotation_intronic_4site ${annotation_intronic_4site} --annotation_gene_transcribed_strands ${annotation_gene_transcribed_strands} --stranded ${stranded}


}



HPB_editing_ratio_clear(){

    local combined_bam=$1
    local outdir=$2
    local outname=$3
    local genome_build_version=$4


dir_vcf=${outdir}/${outname}_tmp/gatk/5-0-0vcf_filter/${outname}

## compute HPB:
## HPB=1000,000,000*DP/unique_mapping_read_base_count

#### 1. calculate total unique_mapping_read base count
samtools stats ${combined_bam} > ${outdir}/${outname}_tmp/unique_mapped-stats.txt
total_base_number=` grep 'total length:' ${outdir}/${outname}_tmp/unique_mapped-stats.txt | awk '{print $4}' `

#### 2. calculate HPB based on the formula and HPB > 3
#### calculate HPB and filter HPB>3, editing_ratio>0.05, ( without QD>2 )
awk -F '\t' -v totalBaseCount=${total_base_number} '{split($8,a,";");for(col in a){if(index(a[col],"DP")>0){split( a[col],dp,"=");print  $0"\t"(1000000000*dp[2]/totalBaseCount);break; }}}'  ${dir_vcf}/${outname}_recommended_final_result_Alu.vcf  | awk '{if($NF>=3){print $0}}' | awk -F '\t' '{split($9,a,":");split($10,b,":");for(idx in a){if(index(a[idx],"AD")>0){ split(b[idx],ref_alt,","); if(ref_alt[1]+ref_alt[2]>0){print $0"\t"ref_alt[2]/(ref_alt[1]+ref_alt[2])   }}}}' | awk '{if($NF>=0.05){print $0}}' | sed 's/\tT\t/\tU\t/g'   > ${outdir}/${outname}_Alu.vcf

awk -F '\t' -v totalBaseCount=${total_base_number} '{split($8,a,";");for(col in a){if(index(a[col],"DP")>0){split( a[col],dp,"=");print  $0"\t"(1000000000*dp[2]/totalBaseCount);break; }}}' ${dir_vcf}/${outname}_recommended_final_result_Repetitive_non-Alu.vcf   | awk '{if($NF>=3){print $0}}' | awk -F '\t' '{split($9,a,":");split($10,b,":");for(idx in a){if(index(a[idx],"AD")>0){ split(b[idx],ref_alt,","); if(ref_alt[1]+ref_alt[2]>0){print $0"\t"ref_alt[2]/(ref_alt[1]+ref_alt[2])   }}}}' | awk '{if($NF>=0.05){print $0}}' | sed 's/\tT\t/\tU\t/g'  > ${outdir}/${outname}_Repetitive_non_Alu.vcf

awk -F '\t' -v totalBaseCount=${total_base_number} '{split($8,a,";");for(col in a){if(index(a[col],"DP")>0){split( a[col],dp,"=");print  $0"\t"(1000000000*dp[2]/totalBaseCount);break; }}}'  ${dir_vcf}/${outname}_recommended_final_result_Nonrepetitive.vcf  | awk '{if($NF>=3){print $0}}' | awk -F '\t' '{split($9,a,":");split($10,b,":");for(idx in a){if(index(a[idx],"AD")>0){ split(b[idx],ref_alt,","); if(ref_alt[1]+ref_alt[2]>0){print $0"\t"ref_alt[2]/(ref_alt[1]+ref_alt[2])   }}}}' | awk '{if($NF>=0.05){print $0}}' | sed 's/\tT\t/\tU\t/g'   > ${outdir}/${outname}_Non_repetitive.vcf





}

preparation_for_plotting(){
	echo "Preparation for histogram."
        GETOPT_ARGS=`getopt  -o i:h  -al inputdir:,help -- "$@"`
	#echo ${GETOPT_ARGS}
        inputdir=""

        eval set -- "$GETOPT_ARGS"
        while [ -n "$1" ]
        do
                case "$1" in
                        -o|--inputdir | -inputdir) inputdir=$2; shift 2;;
                        -h|--help | -help) usage_preparation_for_plotting; exit 1;;
                        --) break ;;
                        *) echo $1,$2, usage_preparation_for_plotting; break ;;
                esac
        done


outnames=` ls ${inputdir}/*.vcf | grep -v 'result_merged_twelve_types_RNA_editing_events_of_all_samples.vcf' | awk -F '/' '{print $NF}' | sed -e 's/_Non_repetitive.vcf//g' -e 's/_Repetitive_non_Alu.vcf//g' -e 's/_Alu.vcf//g' | sort | uniq  `


tmp_file_12_kinds_variant=${inputdir}/tmp_summarized_number_of_twelve_types_RNA_editing_events_of_all_samples.txt
true > ${tmp_file_12_kinds_variant}
final_RNA_variant_result_stats_file=${inputdir}/result_summarized_number_of_twelve_types_RNA_editing_events_of_all_samples.txt
#merged_vcf=${inputdir}/result_merged_twelve_types_RNA_editing_events_of_all_samples.vcf
#true > ${merged_vcf}

while read one_outname
do
	## stats 12 kind variants count in Alu, Repetitive_non_Alu, and Nonrepetitive region
	echo -e "${one_outname}: Alu" >> ${tmp_file_12_kinds_variant}
	cat ${inputdir}/${one_outname}_Alu.vcf | awk '{print $4">"$5}' | sed 's/\*,//g' | sed 's/,.*//g' | sort | uniq -c >> ${tmp_file_12_kinds_variant}
	echo -e "${one_outname}: Nonrepetitive" >> ${tmp_file_12_kinds_variant}
	cat ${inputdir}/${one_outname}_Non_repetitive.vcf | awk '{print $4">"$5}' | sed 's/\*,//g' | sed 's/,.*//g' | sort | uniq -c >> ${tmp_file_12_kinds_variant}
	echo -e "${one_outname}: Repetitive_non_Alu" >> ${tmp_file_12_kinds_variant}
	cat ${inputdir}/${one_outname}_Repetitive_non_Alu.vcf | awk '{print $4">"$5}' | sed 's/\*,//g' | sed 's/,.*//g' | sort | uniq -c >> ${tmp_file_12_kinds_variant}

	## merge all vcf result into one which was used in Manhattan plot
	#cat ${inputdir}/${one_outname}_Alu.vcf | sed 's/$/\t'${one_outname}'\tAlu/g' >> ${merged_vcf}
	#cat ${inputdir}/${one_outname}_Non_repetitive.vcf | sed 's/$/\t'${one_outname}'\tNon_repetitive/g' >> ${merged_vcf} 
	#cat ${inputdir}/${one_outname}_Repetitive_non_Alu.vcf | sed 's/$/\t'${one_outname}'\tRepetitive_non_Alu/g' >> ${merged_vcf}
done <<< "${outnames}"

echo -e "\n\n\n"  >> ${tmp_file_12_kinds_variant}

## transform statistics of 12 types of RNA-editing into excel form
python ${dir_scripts}/summarize_all_12_types_of_RNA_editing.py ${tmp_file_12_kinds_variant} ${final_RNA_variant_result_stats_file}
rm  ${tmp_file_12_kinds_variant}


}



"$@"
