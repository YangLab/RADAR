user_bin=/picb/rnomics4/rotation/fuzhican/bin
dep_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files
trim_fq_path=${dep_path}/fastq

m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
all_srr_list=(SRR8570461 SRR8570463 SRR8570465 TREAT-BE3_rep1 NT-BE3_rep1 TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240 SRR4235528 SRR4235527 20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
all_control_srr_list=(TREAT-Cas9_rep1 NT-Cas9_rep1 20190619_LSQ10_19302 20141031_293FT_ADAR_scr_polyAplus SRR8570460 SRR8570462 SRR8832240)
all_treat_srr_list=(SRR8570461 SRR8570463 SRR8570465 TREAT-BE3_rep1 NT-BE3_rep1)
Six_testdata_srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
#skip: 20190619_LSQ9_19301 SRR8570464
NatMethods_50_brain_srr_list=( ERR030890 ERR030882 SRR039628 SRR039629 SRR039630 SRR039631 SRR039632 SRR039633 SRR085726 SRR085473 SRR087416 SRR309137 SRR309138 SRR309133 SRR309134 SRR309135 SRR309136 SRR036966 SRR014262 SRR085474 SRR309139 SRR309140 SRR085471 SRR309141 SRR309142 SRR309143 SRR309144 SRR085725 SRR090440 SRR090441 SRR090442 SRR107727 SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900 SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306839 SRR306840 SRR306841 SRR306842 SRR306844 SRR309262 )

 
 NatMethods_50_brain_single_end_srr_all=(SRR309139 SRR014262 SRR309136 SRR085474 SRR111899 SRR309142 SRR111900 SRR309143 SRR111902 SRR039630 SRR111906 SRR111905 SRR111903 SRR111895 SRR111898 SRR085471 SRR039632 SRR309135 SRR087416 SRR111904 SRR111897 SRR111901 SRR107727 SRR309140 SRR306839 SRR085725 SRR039633 SRR309141 SRR309134 SRR306844 SRR309138 SRR306841 SRR085473 SRR111907 SRR309133 SRR111896 SRR309144 SRR309137 SRR039631 ERR030890 SRR085726 SRR036966 SRR039629 SRR039628) #44



 NatMethods_50_brain_pair_end_srr_all=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306840 SRR306842 SRR309262 ERR030882) #16
Time_m_path=/data/rnomics6/fuzhican/project/RNA_editing_Time_flow_fine_tune
ABE_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
sh_path=${ABE_m_path}/sh
ref_genome=hg38
sam_markduplicate(){
    srr=TREAT-BE3_rep1
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/bmc_markdup
    # The first sort can be omitted if the file is already name ordered
    samtools sort -@ 3 -n -o ${srr}_combine_readgroup_sort_recal_gatk4_namesort.bam ${srr}_combine_readgroup_sort_recal_gatk4.bam

    # Add ms and MC tags for markdup to use later
    samtools fixmate -m ${srr}_combine_readgroup_sort_recal_gatk4_namesort.bam ${srr}_combine_readgroup_sort_recal_gatk4_namesort_fixmate.bam

    # Markdup needs position order
    samtools sort -@ 3 -o ${srr}_combine_readgroup_sort_recal_gatk4_namesort_fixmate_sort.bam ${srr}_combine_readgroup_sort_recal_gatk4_namesort_fixmate.bam

    # Finally mark duplicates
    samtools markdup ${srr}_combine_readgroup_sort_recal_gatk4_namesort_fixmate_sort.bam ${srr}_combine_readgroup_sort_recal_gatk4_namesort_fixmate_sort_markdup.bam
}
trim_run(){
    srr_list=(SRR4235528 SRR4235527)
    dep_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files
    trim_fq_path=${dep_path}/fastq
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    out_path=${m_path}/fastqc_out
    for p_fq_prefix in ${srr_list[@]}
    do
    {
        #p_fq=${p_fq_prefix}_R1.fastq.gz
        p_trimed_fq_basename=${trim_fq_path}/${p_fq_prefix}_trimed.fastq.gz
        java -jar ${user_bin}/trimmomatic-0.38.jar  PE -threads  5 ${dep_path}/fastq/${p_fq_prefix}_1.fastq.gz ${dep_path}/fastq/${p_fq_prefix}_2.fastq.gz -baseout ${p_trimed_fq_basename} ILLUMINACLIP:/picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36 
        fastqc -o $out_path -d ~/tmp -q ${trim_fq_path}/${p_fq_prefix}_trimed*P.fastq.gz
    }&
    done
    wait
    }
split_bam_by_chr(){
    #m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand
    #vcf_file=$1 #${m_path}/TREAT-BE3_rep1_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf
    bam_file=$1 #${m_path}/TREAT-BE3_rep1_dedupped_split_recal.bam
    bam_file_prefix=$(basename $bam_file|cut -d. -f1)
    out_path=$2 #${m_path}/split_bam_by_chr
    #chr_list=($(awk '{a[$1]++}END{for (i in a){print i}}' $vcf_file))
    for chr in ${chr_list[@]}
    do
    {
        #echo split_bam_by_chr $chr
        samtools1 view -@ 3 -b $bam_file $chr |samtools1 sort -@ 3 -o ${out_path}/${bam_file_prefix}_${chr}.bam
        samtools1 index -@ 3 ${out_path}/${bam_file_prefix}_${chr}.bam
    }
    done
}

merge_num(){
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc/3-0-0Combine_bam
    list1=($(ls *_accepted_hits.sort.bam|awk -F "_accepted" '{print $1}'))
    for srr in ${list1[@]}
    do
    {
        hisat2=$(samtools1 view -@ 3 ${srr}_accepted_hits.sort.bam|wc -l ) &
        bwa=$(samtools1 view -@ 3 ../2-0-0BWA_map/${srr}_unmapped.sort.bam|wc -l) &
        com=$(samtools1 view -@ 3 ${srr}_combine.bam|wc -l)&
        wait
        com_c=`expr $hisat2 + $bwa`;echo $srr $hisat2 $bwa  $com_c $com  
    }
    done
}
extract_plus_readid(){
    bamf=$1
    plus_readid=$2
    tmp_dir=`dirname $plus_readid`
    test -d $tmp_dir|| mkdir -p $tmp_dir
    samtools1 view -@ 3 -F 1024 -f 80 $bamf |awk '{print $1 >"'${plus_readid}'_1_"$3}' &
    samtools1 view -@ 3 -F 1024 -f 128 -F 16 $bamf |awk '{print $1 >"'${plus_readid}'_2_"$3}'
    wait 
    list1=($(ls ${plus_readid}_1_chr*|awk -F"_1_" '{print $2}'))
    list1+=($(ls ${plus_readid}_2_chr*|awk -F"_2_" '{print $2}'))
    chr_list=($(echo ${list1[@]}|awk 'BEGIN{RS=" ";ORS="\n"}{print }'|sort|uniq))
    for chr in ${chr_list[@]}
    do
    {
        cat ${plus_readid}_1_$chr  ${plus_readid}_2_$chr >${plus_readid}_$chr
    }&
    done
}
extract_minus_readid(){
    bamf=$1
    minus_readid=$2
    tmp_dir=`dirname $minus_readid`
    test -d $tmp_dir|| mkdir -p $tmp_dir
    samtools1 view -@ 3 -F 1024 -f 64 -F 16 $bamf |awk '{print $1 >"'${minus_readid}'_1_"$3}' &
    samtools1 view -@ 3 -F 1024 -f 144 $bamf |awk '{print $1 >"'${minus_readid}'_2_"$3}'
    wait
    list1=($(ls ${minus_readid}_1_chr*|awk -F"_1_" '{print $2}'))
    list1+=($(ls ${minus_readid}_2_chr*|awk -F"_2_" '{print $2}'))
    chr_list=($(echo ${list1[@]}|awk 'BEGIN{RS=" ";ORS="\n"}{print }'|sort|uniq))
    for chr in ${chr_list[@]}
    do
    {
        test -e ${minus_readid}_1_$chr||touch ${minus_readid}_1_$chr
        test -e ${minus_readid}_2_$chr||touch ${minus_readid}_2_$chr
        cat ${minus_readid}_1_$chr  ${minus_readid}_2_$chr >${minus_readid}_$chr
    }&
    done
}

pileup_intersection(){
    #need ${minus_readid}_$chr from function plus_minus_extract_readid
    chr=$1
    pos=$2
    bamf=$3
    bamf_basename=`basename ${bamf}`
    tmp_path=$(dirname `dirname ${bamf}`)
    minus_readid=${tmp_path}/minus_readid_path/minus_readid
    pileup_reads=${tmp_path}/`echo ${bamf_basename}|cut -d. -f1`_${pos}_pileup_reads_`random /1..102/`
    #/_256729split_bam_by_chr/TREAT-BE3_rep1_dedupped_split_recal_chr7.bam__pileup_reads_84: No such file or directory
    samtools1 view -F 1024 -@ 3 $bamf ${chr}:${pos}-${pos}|cut -f1 > $pileup_reads
    total_reads_count=`echo $pileup_reads|wc -l `
    plus_mapped_reads_count=`grep -Fxc -f $pileup_reads ${minus_readid}_$chr`
    plus_mapped_fraction=`echo "scale=0;$plus_mapped_reads_count*1000/$total_reads_count"|bc`
    rm $pileup_reads
    echo $plus_mapped_fraction >>${tmp_path}/`echo ${bamf_basename}|awk -F"_chr" '{print $1}'`_pileup_intersection.log
    test $plus_mapped_fraction -gt 700 && exit 1 || exit 0
    #exit 1: $plus_mapped_fraction great than 700, change the base
    
    #the threshold 700/1000 (0.7) come from python ~/bin/REDItoolDenovo.py
    #2013_Bioinformatics_REDItools_ high-throughput RNA editing detection made easy
}


gatk_chr_level_v0(){
    #bash shell solution
    chr=""
    vcf_file_chr=""
    bamf_chr=""
    #while IFS= read -r line
    #do
    # {
     #   pos=echo
      #  #system("fun_test4 10")
        
    #}
    #done >$input
}
gatk_chr_level_v1(){
    #awk solution
    #use function pileup_intersection
    vcf_file_chr=$1
    chr=`basename $vcf_file_chr`
    vcf_file_chr_fine_tune=$tmp_path_vcf_fine_tune/$chr
    bamf_chr=$2
    export -f pileup_intersection
    awk  '
    BEGIN{
        BaseChange["C"]="G";
        BaseChange["G"]="C";
        BaseChange["A"]="T";
        BaseChange["T"]="A";
        FS="\t";
        OFS="\t"
    }
    {
        pos=$2
        change_OrNot=system("pileup_intersection '$chr' " pos " '$bamf_chr'" )
        #change_OrNot = 1 , need change
        #change_OrNot = 0 , need not change
        if (change_OrNot =="1"){
            change_log[$4" "$5" to "BaseChange[$4]" "BaseChange[$5]]++
            $4=BaseChange[$4]
            $5=BaseChange[$5]
            $8=$8";SC=1"
        }else{
            $8=$8";SC=0"
        }
        print $0
    }END{
        for (i in change_log){
            printf("%s: %'\''d\n",i,change_log[i]) >>"'$tmp_path_vcf_fine_tune'/change_log_'$chr'"
        }
    }
    '  $vcf_file_chr > $vcf_file_chr_fine_tune
}
gatk_distinguish_plus_minus_main_flow(){
    test_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk
    srr_list=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240) 
    srr_list=(SRR4235528 SRR4235527)
    srr_list=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    srr_list=(SRR4235528)
    srr_list=(TREAT-Cas9_rep1)
    srr_list=(SRR8570461 SRR8570463 SRR8570465)
    srr_list+=(TREAT-BE3_rep1 NT-BE3_rep1)

    for srr in ${srr_list[@]}
    do
    {
        vcf_file=${m_path}/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf
        bam_file=${m_path}/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
        eval tmp_path=${m_path}/5-0-0vcf_filter/${srr}/distinguish_plus_minus
        bam_file_prefix=$(basename $bam_file|cut -d. -f1)
        #################!Warning!###############
        #Due to the following "eval"s, this function\
        # (gatk_distinguish_plus_minus_main_flow)\
        # can not use multi-threads in out-most loop.
        ##################!Warning!##############
        eval tmp_path_vcf=${tmp_path}/split_vcf_by_chr
        test -d $tmp_path_vcf||mkdir -p $tmp_path_vcf
        eval tmp_path_vcf_fine_tune=${tmp_path}/split_vcf_by_chr_fine_tuned
        test -d $tmp_path_vcf_fine_tune||mkdir -p $tmp_path_vcf_fine_tune
        eval tmp_path_bam=${tmp_path}/split_bam_by_chr
        test -d $tmp_path_bam||mkdir -p $tmp_path_bam
        eval chr_list=($(awk '$0 !~/^#/{a[$1]++}END{for (i in a){print i}}' $vcf_file ))
        eval minus_readid=${tmp_path}/minus_readid_path/minus_readid
        test -d ${tmp_path}/minus_readid_path||mkdir ${tmp_path}/minus_readid_path
        extract_minus_readid $bam_file $minus_readid &
        awk '{print $0 >"'${tmp_path_vcf}'/"$1}' <(awk '$0 !~/^#/' $vcf_file) &
        split_bam_by_chr  $bam_file $tmp_path_bam
        wait
        chr_list=($(ls ${tmp_path_vcf}))
        for chr in ${chr_list[@]}
        do
        {
            test -e ${minus_readid}_$chr||touch ${minus_readid}_$chr
            gatk_chr_level_v1 ${tmp_path_vcf}/$chr  ${tmp_path_bam}/${bam_file_prefix}_${chr}.bam
        }&
        done
        wait
        cat $tmp_path_vcf_fine_tune/chr* >${m_path}/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_strand.vcf
    }
    done


}

gatk_distinguish_plus_minus_main_flow_test(){
    test_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand
    m_path=$test_m_path
    vcf_file=${m_path}/TREAT-BE3_rep1_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf
    bam_file=${m_path}/TREAT-BE3_rep1_dedupped_split_recal.bam
    eval tmp_path=${m_path}/tmp_distinguish_plus_minus
    bam_file_prefix=$(basename $bam_file|cut -d. -f1)
    #################!Warning!###############
    #Due to the following "eval"s, this function\
    # (gatk_distinguish_plus_minus_main_flow)\
    # can not use multi-threads in out-most loop.
    ##################!Warning!##############
    eval tmp_path_vcf=${tmp_path}/split_vcf_by_chr
    test -d $tmp_path_vcf||mkdir -p $tmp_path_vcf
    eval tmp_path_vcf_fine_tune=${tmp_path}/split_vcf_by_chr_fine_tuned
    test -d $tmp_path_vcf_fine_tune||mkdir -p $tmp_path_vcf_fine_tune
    eval tmp_path_bam=${tmp_path}/split_bam_by_chr
    test -d $tmp_path_bam||mkdir -p $tmp_path_bam
    eval chr_list=($(awk '{a[$1]++}END{for (i in a){print i}}' $vcf_file))
    eval minus_readid=${tmp_path}/minus_readid_path/minus_readid
    test -d ${tmp_path}/minus_readid_path||mkdir ${tmp_path}/minus_readid_path

    extract_minus_readid $bam_file $minus_readid &
    awk '{print $0 >"'${tmp_path_vcf}'/"$1}' $vcf_file &
    split_bam_by_chr  $bam_file $tmp_path_bam
    wait
    chr_list=($(ls ${tmp_path_vcf}))
    for chr in ${chr_list[@]}
    do
    {
        gatk_chr_level_v1 ${tmp_path_vcf}/$chr  ${tmp_path_bam}/${bam_file_prefix}_${chr}.bam
    }&
    done
    wait

}
bmc_chr_level_v1(){
    #subfunction of bmc_distinguish_plus_minus_main_flow
    #awk solution
    #use function pileup_intersection
    vcf_file_chr=$1
    chr=`basename $vcf_file_chr`
    vcf_file_chr_fine_tune=$tmp_path_vcf_fine_tune/$chr
    bamf_chr=$2
    export -f pileup_intersection
    awk  '
    BEGIN{
        BaseChange["c"]="G";
        BaseChange["g"]="C";
        BaseChange["a"]="T";
        BaseChange["t"]="A";
        FS="\t";
        OFS="\t"
    }
    {
        split($1,a_tmp_1,":")
        pos=a_tmp_1[2]
        change_OrNot=system("pileup_intersection '$chr' " pos " '$bamf_chr'" )
        #change_OrNot = 1 , need change
        #change_OrNot = 0 , need not change
        if (change_OrNot =="1"){
            change_log[$2" "$17" to "BaseChange[$2]" "BaseChange[$17]]++
            $2=BaseChange[$2]
            $17=BaseChange[$17]
        }
        print $0
    }END{
        for (i in change_log){
            printf("%s: %'\''d\n",i,change_log[i]) >>"'$tmp_path_vcf_fine_tune'/change_log_'$chr'"
        }
    }
    '  $vcf_file_chr > $vcf_file_chr_fine_tune
}
bmc_distinguish_plus_minus_main_flow(){
    test_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc
    srr_list=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240) 
    srr_list+=(SRR4235528 SRR4235527)
    srr_list+=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    for srr in ${srr_list[@]}
    do
    {
        vcf_file=${m_path}/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_recal_1
        bam_file=${m_path}/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
        eval tmp_path=${m_path}/4-0-0Editing_sites/${srr}/distinguish_plus_minus
        bam_file_prefix=$(basename $bam_file|cut -d. -f1)
        #################!Warning!###############
        #Due to the following "eval"s, this function\
        # (gatk_distinguish_plus_minus_main_flow)\
        # can not use multi-threads in out-most loop.
        ##################!Warning!##############
        eval tmp_path_vcf=${tmp_path}/split_vcf_by_chr
        test -d $tmp_path_vcf||mkdir -p $tmp_path_vcf
        eval tmp_path_vcf_fine_tune=${tmp_path}/split_vcf_by_chr_fine_tuned
        test -d $tmp_path_vcf_fine_tune||mkdir -p $tmp_path_vcf_fine_tune
        eval tmp_path_bam=${tmp_path}/split_bam_by_chr
        test -d $tmp_path_bam||mkdir -p $tmp_path_bam
        eval chr_list=($(awk '{split($1,a,":");chr=a[1];tmp_a[chr]++}END{for (i in tmp_a){print i}}' $vcf_file))
        eval minus_readid=${tmp_path}/minus_readid_path/minus_readid
        test -d ${tmp_path}/minus_readid_path||mkdir ${tmp_path}/minus_readid_path
        extract_minus_readid $bam_file $minus_readid &
        awk '{split($1,a,":");chr=a[1];print $0 >"'${tmp_path_vcf}'/"chr}' $vcf_file &
        split_bam_by_chr  $bam_file $tmp_path_bam
        wait
        chr_list=($(ls ${tmp_path_vcf}))
        for chr in ${chr_list[@]}
        do
        {
            test -e ${minus_readid}_$chr||touch ${minus_readid}_$chr
            bmc_chr_level_v1 ${tmp_path_vcf}/$chr  ${tmp_path_bam}/${bam_file_prefix}_${chr}.bam
        }&
        done
        wait
        cat $tmp_path_vcf_fine_tune/chr* >${m_path}/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_recal_strand
    }
    done
}
awk_base_change_percentage(){
    #2019_07_20
    #base_change_percentage.sh
    method=$1 #distinguish input file format, support option: bmc, gatk, avimimic
    input_file=$2
    awk 'BEGIN{OFS="\t"}$0!~/^#/{
        if ("'$method'" == "bmc"){
            base_change_array[$2"_"$(NF-1)]++
        }else if ("'$method'" == "gatk"){
            base_change_array[$4"_"$5]++
        }} END{
            for (i in base_change_array){
                print i,base_change_array[i]/NR
            }
        }
        ' $input_file|sort -k2
}
geneanno(){
    convert2annovar=/picb/rnomics4/rotation/fuzhican/bin/convert2annovar.pl
ANNOVAR_humandb=/picb/rnomics4/rotation/fuzhican/download/annovar/humandb
filter_region_type_list=(Repetitive_non-Alu Alu Nonrepetitive)
#awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17)}' ${bmc_input} >${bmc_output}
srr=SRR8832240
srr=SRR8570460
srr=20141031_293FT_ADAR_scr_polyAplus
filter_region_type=Repetitive_non-Alu
filter_region_type=Nonrepetitive

ref_genome=hg38
vcf_name=${srr}_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_${filter_region_type}.vcf
vcf_name=${srr}_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_Alu.vcf
vcf_name=${srr}_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS.vcf

vcf_filter_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk/5-0-0vcf_filter/${srr}

avinput_name=`echo $vcf_name|cut -d. -f1 `.avinput
$convert2annovar -format vcf4 ${vcf_filter_path}/${vcf_name} > ${vcf_filter_path}/${avinput_name}
ANNOVAR_out=`echo $vcf_name|cut -d. -f1 `_ANNOVAR_annatation_knownGene
annotate_variation.pl --geneanno -dbtype refGene --buildver ${ref_genome} -out ${vcf_filter_path}/${ANNOVAR_out} ${vcf_filter_path}/${avinput_name} ${ANNOVAR_humandb}  
#annotate_variation.pl --geneanno  -exonsort --splicing_threshold 5 -exonicsplicing -dbtype knownGene --buildver ${ref_genome} -out ${ANNOVAR_out} ${bmc_output} ${ANNOVAR_humandb}  
}
test_knwon_SNP(){
    convert2annovar=/picb/rnomics4/rotation/fuzhican/bin/convert2annovar.pl
    ref_genome=hg38
    ANNOVAR_humandb=/picb/rnomics4/rotation/fuzhican/download/annovar/humandb
    test_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/known_SNP_intronic
    vcf_filter_path=$test_path
    srr=chr21
    vcf_name=$srr
    avinput_name=${vcf_name}.avinput
    $convert2annovar -format vcf4 ${vcf_filter_path}/${vcf_name} > ${vcf_filter_path}/${avinput_name}
    ANNOVAR_out=${vcf_filter_path}/${srr}_ANNOVAR_annatation_knownGene
    annotate_variation.pl --geneanno -dbtype refGene --buildver ${ref_genome} -out ${ANNOVAR_out} ${vcf_filter_path}/${avinput_name} ${ANNOVAR_humandb}  
}
awk_annovar_region_percentage(){
    input_file=NT-Cas9_rep1.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5.avinput_recal_ANNOVAR_annatation_knownGene.variant_function
    awk '{
        region_array[$1]++
        }END{
            for (i in region_array){
                print i,region_array[i]*100/NR
            }
        }
        ' $input_file|sort -n -k2
}

distinguish_plus_minus_main_annotation_version(){
    #cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand_annotation
    ref_genome=hg38
    method=$1
    local input_file=$2
    local out_final_part_name=$3
    local tmp_work_path=`dirname $input_file`
    local input_file_name=`basename $input_file`
    local srr=`basename $tmp_work_path`
    test $method == "bmc" && {
        result_path=4-0-0Editing_sites 
        local inter_name=`echo $input_file_name|awk -F ".BQ20o6ES95v2" '{print $2}'`
        }
    test $method == "gatk" && { 
        result_path=5-0-0vcf_filter
        local inter_name=`echo $input_file_name|awk -F "_Haplo" '{print $2}'`
        }   
    #inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf";out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_annotationStrand.vcf"
    #/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc/4-0-0Editing_sites/SRR111896_trimed/SRR111896_trimed.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_annotationStrand_recal

    local input_result=$input_file
    local out_aviput="${tmp_work_path}/${srr}${inter_name}_recal"
    in_bed="${dep_path}/${ref_genome}/ref_all_6.bed"
    out_bed_minus="${dep_path}/${ref_genome}/ref_all_6_minus.bed"
    out_bed_plus="${dep_path}/${ref_genome}/ref_all_6_plus.bed"
    #out_1baseBed="${tmp_work_path}/${srr}${inter_name}_bed"
    #out_1baseBed_bothStrand="${tmp_work_path}/${srr}${inter_name}_bed_plus_minus"
    out_1baseBed_minus="${tmp_work_path}/${srr}${inter_name}_bed_minus"
    out_1baseBed_plus="${tmp_work_path}/${srr}${inter_name}_bed_plus"
    out_strand_result="${out_final_part_name}"
    out_strand_result_minus="${tmp_work_path}/${srr}${inter_name}_annotationStrandminus_recal"
    out_strand_result_plus="${tmp_work_path}/${srr}${inter_name}_annotationStrandplus_recal"
    both_strands_site="${tmp_work_path}/${srr}${inter_name}_both_strand_sites"
    test -e $out_bed_minus||awk '$NF ~"-"'  $in_bed >$out_bed_minus
    test -e $out_bed_plus|| awk '$NF ~"+"'  $in_bed >$out_bed_plus
    test $method == "bmc" && ( test -e $out_aviput || awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $input_result >${out_aviput} )
    test $method == "gatk" && ( test -e $out_aviput ||${user_bin}/convert2annovar.pl -format vcf4 $input_result > ${out_aviput})
    #annotate_variation.pl ${out_aviput} `dirname $in_bed` --buildver ${ref_genome}  -bedfile `basename $in_bed` -dbtype bed -regionanno -out ${out_1baseBed_bothStrand}
    annotate_variation.pl ${out_aviput} `dirname $out_bed_minus` --buildver ${ref_genome}  -bedfile `basename $out_bed_minus` -dbtype bed -regionanno -out ${out_1baseBed_minus}
    annotate_variation.pl ${out_aviput} `dirname $out_bed_plus` --buildver ${ref_genome}  -bedfile `basename $out_bed_plus` -dbtype bed -regionanno -out ${out_1baseBed_plus}
    if [ $method == "bmc" ];then
    {
        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if ( a[$1]){print}}' ${out_1baseBed_plus}.${ref_genome}_bed $input_result >${out_strand_result_plus}
        awk 'BEGIN{
            BaseChange["C"]="G";
            BaseChange["G"]="C";
            BaseChange["A"]="T";
            BaseChange["T"]="A";
            BaseChange["c"]="g";
            BaseChange["g"]="c";
            BaseChange["a"]="t";
            BaseChange["t"]="a";
            OFS="\t"}
            FILENAME==ARGV[1]{a[$3":"$4]++}
            FILENAME==ARGV[2]{
                if ( a[$1]){
                    $2=BaseChange[$2]
                    $(NF-1)=BaseChange[$(NF-1)]
                    print $0
                    }}' ${out_1baseBed_minus}.${ref_genome}_bed $input_result >${out_strand_result_minus}
        awk 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{a[$1]++}FILENAME==ARGV[2]{if (a[$1]){print $1}}' ${out_strand_result_minus} ${out_strand_result_plus} >${both_strands_site} 
        awk 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{a[$0]++}FILENAME==ARGV[2]{if (!a[$1]){print}}FILENAME==ARGV[3]{if (!a[$1]){print}}' ${both_strands_site}  ${out_strand_result_minus} ${out_strand_result_plus} >${out_strand_result}
    }
    elif [ $method == "gatk" ];then
    {
        awk 'BEGIN{FS="\t"}FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if ( a[$1":"$2]){print}}' ${out_1baseBed_plus}.${ref_genome}_bed $input_result >${out_strand_result_plus}
        awk 'BEGIN{
            BaseChange["C"]="G";
            BaseChange["G"]="C";
            BaseChange["A"]="T";
            BaseChange["T"]="A";
            BaseChange["c"]="g";
            BaseChange["g"]="c";
            BaseChange["a"]="t";
            BaseChange["t"]="a";
            FS="\t"
            OFS="\t"}
            FILENAME==ARGV[1]{a[$3":"$4]++}
            FILENAME==ARGV[2]{
                if ( a[$1":"$2]){
                    $4=BaseChange[$4]
                    $5=BaseChange[$5]
                    print $0
                    }}' ${out_1baseBed_minus}.${ref_genome}_bed $input_result >${out_strand_result_minus}
    
        awk 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{a[$1":"$2]++}FILENAME==ARGV[2]{if (a[$1":"$2]){print $1":"$2}}' ${out_strand_result_minus} ${out_strand_result_plus} >${both_strands_site} 
        awk 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{a[$0]++}FILENAME==ARGV[2]{if (!a[$1":"$2]){print}}FILENAME==ARGV[3]{if (!a[$1":"$2]){print}}' ${both_strands_site}  ${out_strand_result_minus} ${out_strand_result_plus} >${out_strand_result}

    }
    fi
    echo $method $srr
    #awk_base_change_percentage $method ${out_strand_result}
    rm  ${out_1baseBed_minus}.${ref_genome}_bed ${out_1baseBed_plus}.${ref_genome}_bed ${out_strand_result_plus} ${out_strand_result_minus} ${both_strands_site} $out_aviput
}
run_2013_Natmethods(){
    method_list=(gatk bmc)
    for srr in ${all_control_srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            #distinguish_plus_minus_main_annotation_version $method $srr
            2013_Natmethods_non-Alu_deSimpleRepeat $method $srr
            
        }
        done
    }
    done

}


 mpileup_mutation_bmc(){
    ref_genome=hg38
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/bmc"
    combine_bam=${m_path}/3-0-0Combine_bam
    log_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/log/${ref_genome}/flow_summary
    dep_path='/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files'
    bwa_index_path=${dep_path}/bwa_mem_index_${ref_genome}
    ref_genome_path=${bwa_index_path}/${ref_genome}_all.fa
    a=(TREAT-Cas9_rep1)
    a+=(NT-BE3_rep1 NT-Cas9_rep1)
    a+=(SRR8570460 SRR8570461 SRR8570462 SRR8570463 SRR8570464 SRR8570465 SRR8832240)
    a=(NT-BE3_rep1 SRR8570461 SRR8570463 SRR8570465)#need do 2019_07_12
    a=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240) 
    a=(SRR4235528 SRR4235527)
    a=(20190619_LSQ10_19302)
    a=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8832240) 
    n=0
    for srr in ${a[@]}
    do
    {
        if [ $n -eq 3 ];then
        n=0
        wait
        fi
        from_bam=${combine_bam}/${srr}_combine_readgroup_sort_recal_gatk4.bam
        #rm ${combine_bam}/${srr}_mpileup_direct_out.gz
        {
            samtools1 mpileup -Q 20 -d 10000000 -OI -f ${ref_genome_path}  $from_bam |awk '{gsub(/\^./,"",$5);if ($5 ~ /[ATCGNactgn]/){print }}'|gzip >${combine_bam}/${srr}_mpileup_mutation.gz
            token=$(printf "%'d" `zcat ${combine_bam}/${srr}_mpileup_mutation.gz|wc -l `)
            echo $srr $token mpileup_mutation_count|tee -a ${log_path}/mpileup_mutation_count_`date +%Y_%m_%d`.log
        }&
        let n+=1
    }
    done
    wait
}
mpileup_mutation_gatk(){
    ref_genome=hg38
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/gatk"
    sam_fine_tune=${m_path}/3-0-0sam_fine-tune_2pass
    log_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/log/${ref_genome}/flow_summary
    dep_path='/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files'
    bwa_index_path=${dep_path}/bwa_mem_index_${ref_genome}
    ref_genome_path=${bwa_index_path}/${ref_genome}_all.fa
    a=(TREAT-Cas9_rep1)
    a+=(NT-BE3_rep1 NT-Cas9_rep1)
    a+=(SRR8570460 SRR8570461 SRR8570462 SRR8570463 SRR8570464 SRR8570465 SRR8832240)
    a=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240)
    a=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    a=(SRR4235528 SRR4235527)
    n=0
    for srr in ${a[@]}
    do
    {
        if [ $n -eq 3 ];then
        n=0
        wait
        fi
        from_bam=${sam_fine_tune}/${srr}_dedupped_split_recal.bam
        #rm ${combine_bam}/${srr}_mpileup_direct_out.gz
        {
            samtools1 mpileup -Q 20 -d 10000000  -OI -f ${ref_genome_path}  $from_bam |awk '{gsub(/\^./,"",$5);if ($5 ~ /[ATCGNactgn]/){print }}'|gzip >${sam_fine_tune}/${srr}_mpileup_mutation.gz
            token=$(printf "%'d" `zcat ${sam_fine_tune}/${srr}_mpileup_mutation.gz|wc -l `)
            echo $srr $token mpileup_mutation_count|tee -a ${log_path}/gatk_mpileup_mutation_count_`date +%Y_%m_%d`.log
        }&
        let n+=1
    }
    done
    wait
}
read_position_gatk(){
    a=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240)
    a+=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    a=(SRR4235528 SRR4235527 SRR8832240) ##need do
    n=0
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk"
    for srr in ${a[@]}
    do 
    {
        mpileup_mutation_out=${m_path}/3-0-0sam_fine-tune_2pass/${srr}_mpileup_mutation.gz
        vcf_file=${m_path}/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf
        tmp_pos_count=${m_path}/6-0-0Summary/read_position_${srr}.txt
        awk 'FILENAME==ARGV[1]{a[$1":"$2]++}FILENAME==ARGV[2]{if ($1":"$2 in a){print $0}}' $vcf_file <(zcat $mpileup_mutation_out)|awk '{split($NF,a,",");for (i in a){pos[a[i]]++}}END{for (i in pos){print i"\t"pos[i]}}' |tee $tmp_pos_count
    }
    done
}
read_position_bmc(){
    a=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240)
    a+=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    a=(SRR4235528 SRR4235527)
    n=0
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc"
    for srr in ${a[@]}
    do
    {
        mpileup_mutation_out=${m_path}/3-0-0Combine_bam/${srr}_mpileup_mutation.gz
        vcf_file=${m_path}/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_recal_1
        tmp_pos_count=${m_path}/6-0-0Summary/read_position_${srr}.txt
        awk 'FILENAME==ARGV[1]{a[$1]++}FILENAME==ARGV[2]{if ($1":"$2 in a){print $0}}' $vcf_file <(zcat $mpileup_mutation_out)|awk '{split($NF,a,",");for (i in a){pos[a[i]]++}}END{for (i in pos){print i"\t"pos[i]}}' > $tmp_pos_count
    }
    done
    wait
}
REDItool_test(){
    ref_genome=hg38
    dep_path='/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files'
    bwa_index_path=${dep_path}/bwa_mem_index_${ref_genome}
    ref_genome_path=${bwa_index_path}/${ref_genome}_all.fa
    test_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/REDItool_test
    bam_file=${test_m_path}/TREAT-BE3_rep1_Aligned.out.bam
    python ~/bin/REDItoolDenovo.py -I $bam_file -f $ref_genome_path -o REDItool_output_3 -t 8 -e -d -l -U AG,TC,CT,GA -p -u -m255 -T6-0 -W -v 1 -n 0.0
}
test_2(){
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand
    samtools1 view -@ 3 -h  -F 1024 -f 64 -F 16 $bamf |tee TREAT-BE3_rep1_dedupped_split_recal_minus.sam|awk '$0 !~/^@/{print $1 >"minus_readid_F16/minus_1_"$3}' &
    samtools1 view -@ 3  -F 1024 -f 144 $bamf |tee -a TREAT-BE3_rep1_dedupped_split_recal_minus.sam|awk '$0 !~/^@/{print $1 >"minus_readid_F16/minus_2_"$3}' &
    wait
    featureCounts -s 0 -p  --fraction -O -T 10 -a /picb/rnomics4/rotation/fuzhican/download/human/hg38/annotation/hg38_gencode.v29.annotation.gtf -t exon -g gene_id -o TREAT-BE3_rep1_dedupped_split_recal_minus_featureCount.txt  TREAT-BE3_rep1_dedupped_split_recal_minus.sam 
    awk '$5 ~"-"{a["-"]+=$NF}$5 ~"+"{a["+"]+=$NF}END{for (i in a){print i,a[i]}}' <(head TREAT-BE3_rep1_dedupped_split_recal_minus_featureCount.txt)
}
test_3(){
    a=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8570464 SRR8832240)
    a+=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    #a+=(SRR4235528 SRR4235527)
    for srr in ${a[@]}
    do
    {
        cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk/5-0-0vcf_filter/${srr}/distinguish_plus_minus/split_vcf_by_chr_fine_tuned
        echo -e "\n\n" $srr
        cat change_log_*|awk -F":" '{a[$1]+=$2}END{for (i in a){print i,a[i]}}'
    }
    done
}
ACTG_averge_dp(){

    awk -F"\t" '$5!~/,/{change_base[$4"-to-"$5]++;split($10,bb,":");dp=bb[3];split(bb[2],cc,",");if (dp!=0){alt_dp=cc[2];er=alt_dp/dp;change_base_er[$4"-to-"$5]+=er}}END{for (i in change_base){print i,change_base[i],change_base_er[i]/change_base[i]}}' <(awk '$0!~/^#/' 20190619_LSQ9_19301_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf)|sort
}
fastqc_1(){
    srr_list=(SRR4235528 SRR4235527)
    srr_list+=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    out_path=${m_path}/fastqc_out
    dep_path=${m_path}/dep_files
    test -d $out_path||mkdir $out_path
    for fq_prefix in ${srr_list[@]};do
    {
        
        fastqc -o $out_path -d ~/tmp -q ${dep_path}/fastq/${fq_prefix}*.fastq.gz 
    }&
    done
}
gatk_mapping_summary(){
    #sh/gatk_flow_summary_CBE.sh
    ref_genome=hg38
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/gatk"
    sam_fine_tune=${m_path}/3-0-0sam_fine-tune_2pass
    echo srr after_dedup mismatch_contained
    srr_list=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    srr_list=(SRR8570460 SRR8570461 SRR8570462 SRR8570463 SRR8570464 SRR8570465 SRR8832240)
    srr_list+=(TREAT-Cas9_rep1 NT-Cas9_rep1 NT-BE3_rep1 TREAT-BE3_rep1)
    #a=(20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302 20190619_LSQ9_19301)
    srr_list+=(SRR4235528 SRR4235527)
    for srr in ${srr_list[@]}
    do
    {
        bma_file=${sam_fine_tune}/${srr}_dedupped_split.bam
        after_dedup=`samtools1 view -@ 4 -F 1024 $bma_file|wc -l`
        mismatch_contained=`samtools1 view -@ 4 -F 1024 $bma_file|awk -F"\t" '$0 ~"nM:i:0"{a[1]=1;if ($1 in a){method_list=(gatk bmc)
    for srr in ${all_control_srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            #distinguish_plus_minus_main_annotation_version $method $srr
            2013_Natmethods_non-Alu_deSimpleRepeat $method $srr
            
        }
        done
    }}else{a[$1]++;mis_contained+=1}}END{print mis_contained}'`
        echo $srr $after_dedup $mismatch_contained

    }
    done
    wait

}
cluster_Alu_nonAlu_nonRepeat_bmc(){
    ref_genome=hg38
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/bmc"
    
    deSNPtag=de_commonSNP
    tag=_recal
    filter_region_type_list=(Repetitive_non-Alu Alu)
    srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    for srr in ${srr_list[@]}
    do
    {
        srr=${srr}_trimed
        #out_path=${m_path}/6-0-0Summary/cluster_Alu_nonAlu_NonRepetitive/${srr}
        #test -d $out_path||mkdir $out_path 
        #20141031_293FT_ADAR_scr_polyAplus.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        editing_sites=${m_path}/4-0-0Editing_sites/${srr}
        out_path=$editing_sites
        bmc_result_file=${editing_sites}/${srr}.BQ20o6ES95v2.allvariants_${deSNPtag}_HPB3_ER5${tag}
        #bmc_result_file=${editing_sites}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
        avi_input=${editing_sites}/${srr}.BQ20o6ES95v2.allvariants_${deSNPtag}_HPB3_ER5.avinput${tag}
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $bmc_result_file >${avi_input}
        for filter_region_type in ${filter_region_type_list[@]}
        do
        {
            out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}${tag}
            out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
            annotate_variation.pl ${avi_input} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_1baseBed}
            awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >${out_resultFile}
            rm ${out_1baseBed}.${ref_genome}_bed
        }
        done
        filter_region_type=All_repetitive
        out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_Nonrepetitive${tag}
        out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
        annotate_variation.pl ${avi_input} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_1baseBed}

        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >$out_resultFile
        rm ${out_1baseBed}.${ref_genome}_bed
    }
    done
}

cluster_Alu_nonAlu_nonRepeat_gatk(){
    ref_genome=hg38
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/gatk"
    out_path=${m_path}/6-0-0Summary/cluster_Alu_nonAlu_NonRepetitive
    test -d $out_path||mkdir $out_path 
    deSNPtag=de_commonSNP
    tag=_2pass
    convert2annovar=/picb/rnomics4/rotation/fuzhican/bin/convert2annovar.pl
    annotate_variation=/picb/rnomics4/rotation/fuzhican/bin/annotate_variation.pl
    filter_region_type_list=(Repetitive_non-Alu Alu)
    srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    for srr in ${srr_list[@]}
    do
    {
        srr=${srr}_trimed
        #20141031_293FT_ADAR_scr_polyAplus_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf
        local vcf_name=${srr}_HaplotypeCaller_Variants${tag}_SNP_${deSNPtag}_hits10.vcf
        #vcf_name=${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
        local avinput_name=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}.avinput`
        vcf_filter_path=${m_path}/5-0-0vcf_filter/${srr}
        out_path=$vcf_filter_path
        $convert2annovar -format vcf4 ${vcf_filter_path}/${vcf_name} > ${out_path}/${avinput_name}
        for filter_region_type in ${filter_region_type_list[@]}
        do
        {
            alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
            alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}.vcf`
            $annotate_variation ${out_path}/${avinput_name} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_path}/${alu_contain_bed}

            awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if ( a[$1"_"$2]){print $0}}' ${out_path}/${alu_contain_bed}.${ref_genome}_bed ${vcf_filter_path}/${vcf_name} >${out_path}/${alu_filtered_vcf}
            rm ${out_path}/$alu_contain_bed.${ref_genome}_bed
        }
        done
        filter_region_type=All_repetitive
        alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
        alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_Nonrepetitive.vcf`
        $annotate_variation ${out_path}/${avinput_name} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_path}/${alu_contain_bed}
        awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if (! a[$1"_"$2]){print $0}}' ${out_path}/${alu_contain_bed}.${ref_genome}_bed ${vcf_filter_path}/${vcf_name} >${out_path}/${alu_filtered_vcf}
        rm ${out_path}/$alu_contain_bed.${ref_genome}_bed
    }
    done
}

core_of_split_non_alu_to_nonAlu_nonRepeat(){
    ref_genome=hg38
    method=$1
    input_file=$2
    deSNPtag=de_commonSNP
    filter_region_type_list=(Repetitive_non-Alu Alu)
    if [ $method == "bmc" ];then
    {
        bmc_result_file=$input_file
        editing_sites=`dirname $bmc_result_file`
        srr=`basename $editing_sites`
        tag=_recal
        avi_input=${editing_sites}/${srr}.BQ20o6ES95v2.allvariants_${deSNPtag}_HPB3_ER5.avinput${tag}
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $bmc_result_file >${avi_input}
        for filter_region_type in ${filter_region_type_list[@]}
        do
        {
            out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}${tag}
            out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
            annotate_variation.pl ${avi_input} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_1baseBed}
            awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >${out_resultFile}
            rm ${out_1baseBed}.${ref_genome}_bed 
        }
        done
        filter_region_type=All_repetitive
        out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_Nonrepetitive${tag}
        out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
        annotate_variation.pl ${avi_input} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_1baseBed}
        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >$out_resultFile
        rm ${out_1baseBed}.${ref_genome}_bed $avi_input
    }
    elif [ $method == "gatk" ];then
    {
        convert2annovar=/picb/rnomics4/rotation/fuzhican/bin/convert2annovar.pl
        annotate_variation=/picb/rnomics4/rotation/fuzhican/bin/annotate_variation.pl
        vcf_name=$input_file
        avinput_name=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}.avinput`
        $convert2annovar -format vcf4 ${vcf_name} > ${avinput_name}
        for filter_region_type in ${filter_region_type_list[@]}
        do
        {
            alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
            alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}.vcf`
            $annotate_variation ${avinput_name} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out  ${alu_contain_bed}

            awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if ( a[$1"_"$2]){print $0}}'  ${alu_contain_bed}.${ref_genome}_bed  ${vcf_name} > ${alu_filtered_vcf}
            rm  ${alu_contain_bed}.${ref_genome}_bed
        }
        done
        filter_region_type=All_repetitive
        alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
        alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_Nonrepetitive.vcf`
        $annotate_variation ${avinput_name} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out  ${alu_contain_bed}
        awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if (! a[$1"_"$2]){print $0}}' ${alu_contain_bed}.${ref_genome}_bed  ${vcf_name} > ${alu_filtered_vcf}
        rm  ${alu_contain_bed}.${ref_genome}_bed ${avinput_name}
    }
    fi
}
tmp_run(){
    a=(EMX1-BE3_rep2 EMX1-nCas9_rep2 NT-BE3_rep2 NT-nCas9_rep2)
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    for i in ${a[@]}
    do
    {
        bash ${m_path}/sh/bmc_HISAT2_BWA_flow_6_27_begin_at_deFPS.sh -q $i -R hg38 -p >${m_path}/log/hg38/bmc/bmc_HISAT2_BWA_flow_6_27_begin_at_deFPS_${i}.log 2>&1
    }
    done
}
pileup_extract_readseq_sam_format(){
    #sub function of prepare_for_blat_2013_Natmethods_version
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    pos=$2
    chr1=`echo $pos|cut -d":" -f1`
    pos0=`echo $pos|cut -d":" -f2|cut -d"-" -f1`
    bam_file=$1
    tmp_path=~/tmp/blat_filter
    input=$3
    readid_file=${tmp_path}/${input}/tmp_mpileup_readid_`echo $pos|cut -d"-" -f1`_$(random)
    sample_sam_file=${tmp_path}/${input}/view_cluster_sample_`echo $pos|cut -d"-" -f1`_$(random).sam

    samtools1 view -F 1024 $bam_file $pos > ${sample_sam_file} &
    #samtools1 mpileup -Q 20 -O -r $pos -f $ref_genome_path --output-QNAME $bam_file |awk -F"\t" '{
    grep "${chr1}"$'\t'"${pos0}"$'\t' ${tmp_path}/${input}/mpilup_output/${chr1}_mpileup_output|awk -F"\t" '{
        gsub(/\$/,"",$5)
        gsub(/\^./,"",$5)
        split($7,readpos_array,",")
        split($8,readid_array,",")
        for ( i in readpos_array){
            if (readpos_array[i]>6){readpos_preserved[i]=i}
        }
        readbase=$5
        while  (readbase ~ /[AaTtGgCcNn]/  ){
            match(readbase,/[AaTtGgCcNn]/)
            readbase_preserved[RSTART]=RSTART
            sub(/[AaTtGgCcNn]/,"Z",readbase)
             }
        for (i in readbase_preserved){
            if (i in readpos_preserved){
                print readid_array[i]
            }
        }
        }
        ' > ${readid_file} 
    
    wait
    #awk 'FILENAME==ARGV[1]{readid[$0]++}FILENAME==ARGV[2]{if ($1 in readid){seq_fore7=substr($10,1,7);print ">"$1"_"seq_fore7"\n" $10}}' ${readid_file} ${sample_sam_file} 
    awk 'FILENAME==ARGV[1]{readid[$0]++}FILENAME==ARGV[2] &&readid[$1]' ${readid_file} ${sample_sam_file} 
    rm $readid_file $sample_sam_file
}
pileup_extract_readseq(){
    #sub function of prepare_for_blat_2013_Natmethods_version
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    pos=$2
    chr1=`echo $pos|cut -d":" -f1`
    pos0=`echo $pos|cut -d":" -f2|cut -d"-" -f1`
    bam_file=$1
    tmp_path=~/tmp/blat_filter
    input=$3
    readid_file=${tmp_path}/${input}/tmp_mpileup_readid_`echo $pos|cut -d"-" -f1`_$(random)
    sample_sam_file=${tmp_path}/${input}/view_cluster_sample_`echo $pos|cut -d"-" -f1`_$(random).sam

    samtools1 view -F 1024 $bam_file $pos > ${sample_sam_file} &
    #samtools1 mpileup -Q 20 -O -r $pos -f $ref_genome_path --output-QNAME $bam_file |awk -F"\t" '{
    grep "${chr1}"$'\t'"${pos0}"$'\t' ${tmp_path}/${input}/mpilup_output/${chr1}_mpileup_output|awk -F"\t" '{
        gsub(/\$/,"",$5)
        gsub(/\^./,"",$5)
        split($7,readpos_array,",")
        split($8,readid_array,",")
        for ( i in readpos_array){
            if (readpos_array[i]>6){readpos_preserved[i]=i}
        }
        readbase=$5
        while  (readbase ~ /[AaTtGgCcNn]/  ){
            match(readbase,/[AaTtGgCcNn]/)
            readbase_preserved[RSTART]=RSTART
            sub(/[AaTtGgCcNn]/,"Z",readbase)
             }
        for (i in readbase_preserved){
            if (i in readpos_preserved){
                print readid_array[i]
            }
        }
        }
        ' > ${readid_file} 
    
    wait
    awk 'FILENAME==ARGV[1]{readid[$0]++}FILENAME==ARGV[2]{if ($1 in readid){seq_fore7=substr($10,1,7);print ">"$1"_"seq_fore7"\n" $10}}' ${readid_file} ${sample_sam_file} 
    #awk 'FILENAME==ARGV[1]{readid[$0]++}FILENAME==ARGV[2] &&readid[$1]' ${readid_file} ${sample_sam_file} 
    rm $readid_file $sample_sam_file
}
prepare_for_blat_science_comment_version(){
    #blat:Standalone BLAT v. 36x2
    #blat  -makeOoc=11.ooc /picb/rnomics1/database/Human/hg38/genome/hg38_all.fa  /picb/rnomics1/database/Human/hg38/genome/hg38_all.fa out.psl
    input=$1
    ooc=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/11.ooc
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    ref_genome_path_blat=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.2bit
    output_fa=tmp_`date +%s`.fa
    #faToTwoBit $ref_genome_path  `echo $ref_genome_path|cut -d. -f1`.2bit
    while IFS= read -r line
    do
    {
        pos=`echo $line|cut -d $'\t' -f1`
        chr1=`echo $pos|cut -d":" -f1`
        pos_0=`echo $pos|cut -d":" -f2`
        pos_1=`expr $pos_0 - 25`
        pos_2=`expr $pos_0 + 24`
        samtools1 faidx $ref_genome_path ${chr1}:${pos_1}-${pos_2} >>$output_fa
    }
    done < $input
    
    blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -noHead $ref_genome_path $output_fa intern_out.psl
    #gfServer start 127.0.0.1 49252 -repMatch=2253 -stepSize=5 -log=untrans.log $ref_genome_path_blat
    pslScore.pl intern_out.psl >out_with_score.psl
}
psl_filter(){
    #sub function of prepare_for_blat_2013_Natmethods_version
    ###per site in bmc/gatk
    psl_with_score_file=$1
    pos_ori=$2 #chr1:367182
    chr1=`echo $pos_ori|cut -d":" -f1`
    pos_0_0=$(expr `echo $pos_ori|cut -d":" -f2` - 1)
    tmp_a=`awk '{
        count[$4]++
        if ($4 in fa_title_array && count[$4]==2){
            ratio=$5/fa_title_array[$4]
            if (ratio < 0.95){
                blat_survive[$4]++
            }
        }
        else{
            if ("'$chr1'"==$1 && $2 <= '$pos_0_0' &&'$pos_0_0' < $3){
                fa_title_array[$4]=$5
            } 
        }
        }END{
            #print length(blat_survive)
            for (i in count){ ##for read only hit one time! 
                if (count[i]==1 && "'$chr1'"==$1 && $2 <= '$pos_0_0' &&'$pos_0_0' < $3){
                    blat_survive[i]++
                }
            }
            #print length(blat_survive),length(count)
            if (2*length(blat_survive)>length(count)){
                print 1
            }else{print 0}
            }' $psl_with_score_file`
    echo $tmp_a
}
base_oversee(){
    #require function: psl_filter pileup_extract_readseq
    Time_m_path=/data/rnomics6/fuzhican/project/RNA_editing_Time_flow_fine_tune
    ABE_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    
    sh_path=${ABE_m_path}/sh
    ref_genome=hg38
    srr="Human_Forebrain_10Week_PostConception_Female_rep0"
    bam_file=${Time_m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
    input=$1
    method=bmc

    input_basename=`basename $input`
    
    tmp_path=~/tmp/blat_filter
    test -e $out_blat_filter && rm $out_blat_filter
    test -d ${tmp_path}/${input_basename}/split_chr && rm ${tmp_path}/${input_basename}/split_chr/* || mkdir -p ${tmp_path}/${input_basename}/split_chr
    test -d ${tmp_path}/${input_basename}/mpilup_output && rm ${tmp_path}/${input_basename}/mpilup_output/* || mkdir -p ${tmp_path}/${input_basename}/mpilup_output
    test $method = "bmc" && {
        awk '{split($1,aa,":");print aa[1]" "aa[2] >"'${tmp_path}'/'${input_basename}'/split_chr/"aa[1]}' $input 
    }
    test $method = "gatk" && {
        awk '{print $1" "$2 >"'${tmp_path}'/'${input_basename}'/split_chr/"$1}' $input 
    }
    
    chr_list=($(ls ${tmp_path}/${input_basename}/split_chr/))
    n=0
    for chr in ${chr_list[@]}
    do
    {
        if [ "$n" -eq "5" ];then
        n=0
        wait
        fi
        let n+=1
        samtools1 mpileup -Q 20 -O -f $ref_genome_path --output-QNAME -r $chr -l ${tmp_path}/${input_basename}/split_chr/$chr $bam_file > ${tmp_path}/${input_basename}/mpilup_output/${chr}_mpileup_output &
    }
    done
    ############
    if [ $method == "bmc" ];then
    {
        pos_ori=`echo $line|cut -d $' ' -f1` #chr1:367182 ###gatk\bmc
    }
    elif [ $method == "gatk" ];then
    {
        pos_ori=`echo $line|awk '{print $1":"$2}'` #chr1:367182 ###gatk\bmc
    }
    fi
    chr1=`echo $pos_ori|cut -d":" -f1`  ###gatk\bmc
    pos_0=`echo $pos_ori|cut -d":" -f2`  ###gatk\bmc
    output_fa=${tmp_path}/${input_basename}/tmp_${pos_ori}_$(random).fa
    out_client_psl=${tmp_path}/${input_basename}/tmp_out_client_${pos_ori}_$(random).psl
    out_with_score_psl=${tmp_path}/${input_basename}/tmp_out_with_score_${pos_ori}_$(random).psl

    pileup_extract_readseq  $bam_file ${pos_ori}-${pos_0} $input_basename >$output_fa 
    pileup_extract_readseq_sam_format  $bam_file ${pos_ori}-${pos_0} $input_basename

    gfClient 127.0.0.1 49252 ""  -nohead $output_fa $out_client_psl >/dev/null

    pslScore.pl $out_client_psl |awk 'BEGIN{OFS="\t"}{sub(/:[0-9]+-[0-9]+/,"",$4);print}'|sort  -k4,4 -k5nr,5 >$out_with_score_psl

    filter_or_not=$(psl_filter $out_with_score_psl $pos_ori)

}

prepare_for_blat_2013_Natmethods_version(){
    #blat:Standalone BLAT v. 36x2
    #blat  -makeOoc=11.ooc /picb/rnomics1/database/Human/hg38/genome/hg38_all.fa  /picb/rnomics1/database/Human/hg38/genome/hg38_all.fa out.psl
    input=$1
    input_basename=`basename $input`
    #ooc=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/11.ooc
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    ref_genome_path_blat=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.2bit
    
    tmp_path=~/tmp/blat_filter
    out_blat_filter=$3
    bam_file=$2
    method=$4
    test -e $out_blat_filter && rm $out_blat_filter
    test -d ${tmp_path}/${input_basename}/split_chr && rm ${tmp_path}/${input_basename}/split_chr/* || mkdir -p ${tmp_path}/${input_basename}/split_chr
    test -d ${tmp_path}/${input_basename}/mpilup_output && rm ${tmp_path}/${input_basename}/mpilup_output/* || mkdir -p ${tmp_path}/${input_basename}/mpilup_output
    test $method == "bmc" && {
        awk '{split($1,aa,":");print aa[1]" "aa[2] >"'${tmp_path}'/'${input_basename}'/split_chr/"aa[1]}' $input 
    }
    test $method == "gatk" && {
        awk '{print $1" "$2 >"'${tmp_path}'/'${input_basename}'/split_chr/"$1}' $input 
    }
    
    chr_list=($(ls ${tmp_path}/${input_basename}/split_chr/))
    n=0
    for chr in ${chr_list[@]}
    do
    {
        if [ "$n" -eq "5" ];then
        n=0
        wait
        fi
        let n+=1
        samtools1 mpileup -Q 20 -O -f $ref_genome_path --output-QNAME -r $chr -l ${tmp_path}/${input_basename}/split_chr/$chr $bam_file > ${tmp_path}/${input_basename}/mpilup_output/${chr}_mpileup_output &
    }
    done

    #faToTwoBit $ref_genome_path  `echo $ref_genome_path|cut -d. -f1`.2bit
    n=0
    while IFS= read -r line
    do
    if [ "$n" -eq "10" ];then
    {
        n=0
        wait
    }
    fi
    let n+=1
    {
        #pos=chr1:187497
        #echo "`date`_token_0"
        if [ $method == "bmc" ];then
        {
            pos_ori=`echo $line|cut -d $' ' -f1` #chr1:367182 ###gatk\bmc
        }
        elif [ $method == "gatk" ];then
        {
            pos_ori=`echo $line|awk '{print $1":"$2}'` #chr1:367182 ###gatk\bmc
        }
        fi
        chr1=`echo $pos_ori|cut -d":" -f1`  ###gatk\bmc
        pos_0=`echo $pos_ori|cut -d":" -f2`  ###gatk\bmc
        output_fa=${tmp_path}/${input_basename}/tmp_${pos_ori}_$(random).fa
        out_client_psl=${tmp_path}/${input_basename}/tmp_out_client_${pos_ori}_$(random).psl
        out_with_score_psl=${tmp_path}/${input_basename}/tmp_out_with_score_${pos_ori}_$(random).psl
        #echo "pileup_extract_readseq  $bam_file ${pos_ori}-${pos_0}"
        #echo "`date`_token_1"
        pileup_extract_readseq  $bam_file ${pos_ori}-${pos_0} $input_basename >$output_fa 
        #echo "`date`_token_2"
        gfClient 127.0.0.1 49252 ""  -nohead $output_fa $out_client_psl >/dev/null
        #echo "`date`_token_3"
        #blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -noHead $ref_genome_path_blat $output_fa intern_out.psl
        pslScore.pl $out_client_psl |awk 'BEGIN{OFS="\t"}{sub(/:[0-9]+-[0-9]+/,"",$4);print}'|sort  -k4,4 -k5nr,5 >$out_with_score_psl
        #echo "`date`_token_4"
        filter_or_not=$(psl_filter $out_with_score_psl $pos_ori)
        #echo "`date`_token_5"
        #echo filter_or_not,$filter_or_not
        if [ "$filter_or_not" == "1" ];then
        {
            #echo 'yyy'
            echo "$line" >>$out_blat_filter
        }
        fi
        rm $output_fa $out_client_psl $out_with_score_psl
        #echo "`date`_token_6"
    } &
    done < "$input"
    wait
    rm -r ${tmp_path}/${input_basename}

}

2013_Natmethods_variantsReads_misFrequency(){
    #to all: Alu , rep non-alu ,non-rep
    ##we required at least three variant reads and mismatch frequency ≥0.1
    method=$1
    srr=$2
    if [ $method == "bmc" ];then
    {
        input_result_file="20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        awk '$NF*$3 >3 && $NF>.1' $input_result_file
    }
    fi
}
tmp_IFS_test(){
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/test_IFS
    input_file=tmp_1
    while IFS= read -r line
    do
    {
        echo "$line" >>tmp_2 ##this syntax is ok!
    }
    done <"$input_file"
}

2013_NatMethods_filter_regions_in_bed(){
    ###this function filtering sites localed in regions specified within in_bed.
    #in detail, this function can be use to filter Alu, Simple Repeat, intronic.
    ref_genome=hg38
    method=$1 #bmc or gatk
    in_bed=$2 #"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    input_result=$3 #absolute path of input result file
    output_result=$4 #absolute path of output result file
    filter_or_not=$5
    tmp_work_path=`dirname $input_result`
    srr=`basename $tmp_work_path` #####only suit for result_file deposited at director named by srr
    inter_name=`echo $input_result|awk -F "$srr" '{print $3}'`
    out_aviput="${tmp_work_path}/${srr}${inter_name}_recal" #temp internal file
    out_1baseBed="${tmp_work_path}/${srr}${inter_name}_1baseBed_recal" #temp internal file
    test $method == "bmc" && ( test -e $out_aviput || awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $input_result >${out_aviput} )
    test $method == "gatk" && ( test -e $out_aviput ||${user_bin}/convert2annovar.pl -format vcf4 $input_result > ${out_aviput})
    annotate_variation.pl ${out_aviput} `dirname $in_bed` --buildver ${ref_genome}  -bedfile `basename $in_bed` -dbtype bed -regionanno -out ${out_1baseBed}
    if [ "$filter_or_not" == "" ];then
    {
        if [ $method == "bmc" ];then
        {
            awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print}}' ${out_1baseBed}.${ref_genome}_bed $input_result >${output_result}
        }
        elif [ $method == "gatk" ];then
        {
            awk 'BEGIN{FS="\t"}FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1":"$2]){print}}' ${out_1baseBed}.${ref_genome}_bed $input_result >${output_result}
        }
        fi
    }
    else
    {
        if [ $method == "bmc" ];then
        {
            awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (  a[$1]){print}}' ${out_1baseBed}.${ref_genome}_bed $input_result >${output_result}
        }
        elif [ $method == "gatk" ];then
        {
            awk 'BEGIN{FS="\t"}FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (  a[$1":"$2]){print}}' ${out_1baseBed}.${ref_genome}_bed $input_result >${output_result}
        }
        fi
    }
    fi
    rm ${out_1baseBed}.${ref_genome}_bed $out_aviput
    echo $method $srr
    awk_base_change_percentage $method ${output_result}
}

run_2013_NatMethods_filter_regions_in_bed(){
    method_list=(gatk bmc)
    ref_genome=hg38
    m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target"
    for srr in ${Six_testdata_srr_list[@]}
    do
    {
        #srr=${srr}_trimed
        for method in ${method_list[@]}
        do
        {
            run_2013_NatMethods_filter_regions_in_bed_core $method $srr $m_path
        }
        done
    } &
    done
    wait 
}
run_2013_NatMethods_filter_regions_in_bed_core_start_form_hits10(){ 
    ref_genome=hg38  
    method=$1 
    srr=$2
    m_path=$3
    in_bed=${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed
    test $method == "bmc" && {
        #20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal"
        }
    test $method == "gatk" && { 
        #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf";
        out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu.vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result
    in_bed="${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    test $method == "bmc" && {
        #20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_recal" #.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal"
        }
    test $method == "gatk" && { 
        #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf";
        out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu.vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result

    in_bed=${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_intronic_4site.bed #${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed #"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    test $method == "bmc" && {
        #20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_recal" #.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal"
        }
    test $method == "gatk" && { 
        #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf";
        out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu.vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result



    if [ $method == "bmc" ];then
    {
        bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
        
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_recal
        output_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
        2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
        echo -e "$method\t$srr"
        awk_base_change_percentage $method $output_filter
    }
    elif [ $method == "gatk" ];then
    {
        bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp.vcf
        output_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
        2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
        echo "$method\t$srr"
        awk_base_change_percentage $method $output_filter
    }
    fi

    ################################ filter blat ################################
    if [ $method == "bmc" ];then
    {
        bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
        
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
        out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_recal
        prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter $method
    }
    elif [ $method == "gatk" ];then
    {
        bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
        out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat.vcf
        prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter $method
    }
    fi
    ################################strandness################################
    if [ $method == "bmc" ];then
    {
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_recal
        out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_recal
        distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter
    }
    elif [ $method == "gatk" ];then
    {
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat.vcf
        out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.vcf
        distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter
    }
    fi
    ################################split_non-alu_to_nonAlu_nonRepeat#################################
    if [ $method == "bmc" ];then
    {
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_recal
        core_of_split_non_alu_to_nonAlu_nonRepeat $method $bmc_file #output name is auto determined by input 
    }
    elif [ $method == "gatk" ];then
    {
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.vcf
        core_of_split_non_alu_to_nonAlu_nonRepeat $method $bmc_file #output name is auto determined by input 
    }
    fi
    ################################reshape to avimimic################################
    #python solution
    #be_off_target_huge.py function:reshape_avimimic_filter_region_type
}

2013_NatMethods_filter_SNP(){
    #download cat 20170504_GRCh38_positions_manifest.txt|grep '_sites'|cut -f1|xargs -I{} wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/\{\}
    #link
    ref_genome=hg38  
    method=$1 
    input_file=$2
    output_file=$3
    SNP_path=$4 #/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/dbSNP_b151/split_chr
    work_path=`dirname $input_file`
    input_file_basename=`basename $input_file`
    n=0
    test -d ${work_path}/tmp_deknownSNP/chr  ||mkdir -p ${work_path}/tmp_deknownSNP/chr 
    test -d ${work_path}/tmp_deknownSNP/tmp_result  ||mkdir -p ${work_path}/tmp_deknownSNP/tmp_result
    test $method == "bmc" && {
        awk -F ":" '{print $0 > "'${work_path}'/tmp_deknownSNP/chr/"$1}' ${input_file}
        chr_list=($(ls ${work_path}/tmp_deknownSNP/chr|sort -k1.4n))
        for chr in ${chr_list[@]}
        do
        if [ $n -eq 5 ];then
        n=0
        wait
        fi
        let n+=1
        {
            echo $chr| grep -q "_"  &&continue
            if [ -e "${SNP_path}/${chr}.gz" ];then
            awk 'FILENAME==ARGV[1]{array_tmp[$1":"$2]=1}FILENAME==ARGV[2]{chrn=substr($1,4,length($1));if (! array_tmp[chrn]){print}}'  <(zcat ${SNP_path}/${chr}.gz) ${work_path}/tmp_deknownSNP/chr/${chr} >${work_path}/tmp_deknownSNP/tmp_result/${input_file_basename}_${chr}  
            else
            echo "$chr in $input_file do not exists in $SNP_path"
            fi
        }&
        done
        
    }
    test $method == "gatk" && {
        awk '$0 !~/^#/ && $0 !~/^$/{print  > "'${work_path}'/tmp_deknownSNP/chr/"$1}' ${input_file}
        chr_list=($(ls ${work_path}/tmp_deknownSNP/chr|sort -k1.4n))
        for chr in ${chr_list[@]}
        do
        if [ $n -eq 5 ];then
        n=0
        wait
        fi
        let n+=1
        {
            echo $chr| grep -q "_"  &&continue
            if [ -e "${SNP_path}/${chr}.gz" ];then
            awk 'FILENAME==ARGV[1]{array_tmp[$1":"$2]=1}FILENAME==ARGV[2]{chrn=substr($1,4,length($1));if (! array_tmp[chrn":"$2]){print}}'  <(zcat ${SNP_path}/${chr}.gz) ${work_path}/tmp_deknownSNP/chr/${chr} >${work_path}/tmp_deknownSNP/tmp_result/${input_file_basename}_${chr} 
            else
            echo "$chr in $input_file do not exists in $SNP_path"
            fi
        }&
        done
        
    }
    wait
    cat ${work_path}/tmp_deknownSNP/tmp_result/* >$output_file
    rm -r ${work_path}/tmp_deknownSNP

     
}

run_2013_NatMethods_filter_SNP_begin_at_SNP(){
    method_list=(gatk)
    filter_region_type_list_part=(Repetitive_non-Alu Nonrepetitive)
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    sh_path=`pwd`/sh
    srr_list=(SRR8570462 SRR8832240 SRR8570460 TREAT-Cas9_rep1 NT-Cas9_rep1 20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302)
    srr_list=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306840 SRR306842 SRR309262)
    srr_list=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937)
    srr_list=(SRR6191044 SRR6191045 SRR6191046 20190619_LSQ10_19302_derRNA 20190619_LSQ9_19301_derRNA)
    srr_list=(SRR6191044 SRR6191045)
    need_combine=(ERR030890 SRR039628 SRR039630 SRR039632 SRR309137 SRR309133 SRR309135 SRR309139 SRR309141 SRR309143 ERR030882 SRR039629 SRR039631 SRR039633 SRR309138 SRR309134 SRR309136 SRR309140 SRR309142 SRR309144)
    #trimed_or_not="_shorted"
    #for srr in ${NatMethods_50_brain_single_end_srr_all[@]}
    #for srr in ${NatMethods_50_brain_pair_end_srr_all[@]}
    for srr in ${srr_list[@]}
    do
    {
        echo ${need_combine[@]}|grep -q $srr && continue
        srr=${srr}${trimed_or_not}
        for method in ${method_list[@]}
        do
        {
            test $method == "gatk" && {
                suffix=".vcf"
                tmp_work_path=${m_path}/hg38/${method}/5-0-0vcf_filter/${srr}           
                inter_name="_HaplotypeCaller_Variants_2pass_SNP"                
            }
            test $method == "bmc" && {
                suffix="_recal"
                tmp_work_path=${m_path}/hg38/${method}/4-0-0Editing_sites/${srr}           
                inter_name=".BQ20o6ES95v2.allvariants_deFPS"
                ln -s ${tmp_work_path}/${srr}.BQ20o6ES95v2.allvariants_recal_deFPS $input_file
            }
            input_file=${tmp_work_path}/${srr}${inter_name}${suffix}
            output_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151${suffix}
            test -e $input_file ||{
                echo "$input_file do not exists,skip it!"
                continue
            }
            test -e ${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS.vcf && {
                echo "${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS.vcf exists,skip it!"
                continue
            }
            2013_NatMethods_filter_SNP gatk $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/dbSNP_b151/split_chr

            input_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151.vcf
            output_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes.vcf
            2013_NatMethods_filter_SNP gatk $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/1000genomes/split_chr

            input_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes.vcf
            output_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS.vcf
            2013_NatMethods_filter_SNP gatk $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/EVS/split_chr
        }
        done
    }
    done

}
tmp_run_2013_NatMethods_filter_regions_in_bed_core_flexible(){
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    sh_path=${m_path}/sh
    srr_list=(SRR8570462 SRR8832240 SRR8570460 TREAT-Cas9_rep1 NT-Cas9_rep1 20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302)
    need_combine=(ERR030890 SRR039628 SRR039630 SRR039632 SRR309137 SRR309133 SRR309135 SRR309139 SRR309141 SRR309143 ERR030882 SRR039629 SRR039631 SRR039633 SRR309138 SRR309134 SRR309136 SRR309140 SRR309142 SRR309144)
    #srr_list=(SRR8570462)
    #srr_list=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306840 SRR306842 SRR309262)
    srr_list=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937)
    srr_list=(20190619_LSQ10_19302 20141031_293FT_ADAR_scr_polyAplus SRR8832240 TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8832240)
    srr_list=(20190619_LSQ9_19301)
    srr_list=(SRR6191044 SRR6191045)
    #trimed_or_not="_shorted"
    n=0
    #for srr in ${NatMethods_50_brain_single_end_srr_all[@]}
    #for srr in ${NatMethods_50_brain_pair_end_srr_all[@]}
    for srr in ${srr_list[@]}
      
    do
    {
        if [ $n == 2 ];then
        n=0
        wait
        fi
        let n+=1
        echo ${need_combine[@]}|grep -q $srr && continue
        srr=${srr}${trimed_or_not}
        tmp_work_path=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}
        #awk -F ":" '$7>=10' ${tmp_work_path}/${srr}_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS.vcf >${tmp_work_path}/${srr}_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10.vcf
        bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr "/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target" "_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS" >log/hg38/gatk/${srr}_run_2013_NatMethods_filter_regions_in_bed_core_flexible_`date +%Y_%m_%d`.log 2>&1 &
        inter_name="_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand"
        suffix=".vcf"
        input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
        #core_of_split_non_alu_to_nonAlu_nonRepeat gatk $input_result
    }
    done
    wait
}
run_2013_NatMethods_filter_SNP_flexible(){
    ref_genome=hg38  
    method=$1 
    srr=$2
    m_path=$3
    tmp_inter_name=$4
    sh_path=${ABE_m_path}/sh
    test $method == "gatk" && {
        suffix=".vcf"
        tmp_work_path=${m_path}/hg38/${method}/5-0-0vcf_filter/${srr}           
        inter_name="_HaplotypeCaller_Variants_2pass_SNP"                
    }
    test $method == "bmc" && {
        suffix="_recal"
        tmp_work_path=${m_path}/hg38/${method}/4-0-0Editing_sites/${srr}           
        inter_name=".BQ20o6ES95v2.allvariants_deFPS"
        
    }
    if [ "$tmp_inter_name" != "" ];then
    inter_name=$tmp_inter_name
    fi
    input_file=${tmp_work_path}/${srr}${inter_name}${suffix}
    #ln -s ${tmp_work_path}/${srr}.BQ20o6ES95v2.allvariants_recal_deFPS $input_file
    output_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151${suffix}
    test -e $input_file ||{
        echo "$input_file do not exists,skip it!"
        continue
    }
    test -e ${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix} && tmp_a=`head ${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix}`
    test "$tmp_a" != ""  && {
        echo "${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix} exists,skip it!"
        continue
    }


    bash ${sh_path}/support_flow_19_7_8.sh 2013_NatMethods_filter_SNP $method $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/dbSNP_b151/split_chr

    input_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151${suffix}
    output_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes${suffix}
    bash ${sh_path}/support_flow_19_7_8.sh 2013_NatMethods_filter_SNP $method $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/1000genomes/split_chr

    input_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes${suffix}
    output_file=${tmp_work_path}/${srr}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix}
    bash ${sh_path}/support_flow_19_7_8.sh 2013_NatMethods_filter_SNP $method $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/EVS/split_chr
}
run_2013_NatMethods_filter_commonSNP_flexible(){
    ref_genome=hg38  
    method=$1 
    srr=$2
    m_path=$3
    tmp_inter_name=$4
    sh_path=${ABE_m_path}/sh
    test $method == "gatk" && {
        suffix=".vcf"
        tmp_work_path=${m_path}/hg38/${method}/5-0-0vcf_filter/${srr}           
        inter_name="_HaplotypeCaller_Variants_2pass_SNP"                
    }
    test $method == "bmc" && {
        suffix="_recal"
        tmp_work_path=${m_path}/hg38/${method}/4-0-0Editing_sites/${srr}           
        inter_name=".BQ20o6ES95v2.allvariants_deFPS"
    }
    if [ "$tmp_inter_name" != "" ];then
    inter_name=$tmp_inter_name

    fi
    input_file=${tmp_work_path}/${srr}${inter_name}${suffix}
    #ln -s ${tmp_work_path}/${srr}.BQ20o6ES95v2.allvariants_recal_deFPS $input_file
    output_file=${tmp_work_path}/${srr}${inter_name}_decommonSNP${suffix}
    test -e $input_file ||{
        echo "$input_file do not exists,skip it!"
        exit
    }
    test -e ${output_file} && tmp_a=`head ${output_file}`
    test "$tmp_a" != ""  && {
        echo "${output_file} exists,skip it!"
        exit
    }
    test $method == "gatk" && {
        awk 'FILENAME==ARGV[1]{array_tmp[$1":"$2]=1}FILENAME==ARGV[2]{if (! array_tmp[$1":"$2]){print}}' <(bcftools view -H ${dep_path}/${ref_genome}/NCBI_dbSNP_b151_common_${ref_genome}.vcf)  $input_file >$output_file
    }
    test $method == "bmc" && {
        awk 'FILENAME==ARGV[1]{array_tmp[$1":"$2]=1}FILENAME==ARGV[2]{if (! array_tmp[$1]){print}}' <(bcftools view -H ${dep_path}/${ref_genome}/NCBI_dbSNP_b151_common_${ref_genome}.vcf) $input_file >$output_file
    }

}
run_2013_NatMethods_filter_SNP(){
    method_list=(gatk)
    filter_region_type_list_part=(Repetitive_non-Alu Nonrepetitive)
    filter_region_type_list_part=(1)
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    sh_path=`pwd`/sh
    srr_list=(SRR8570462 SRR8832240 SRR8570460 TREAT-Cas9_rep1 NT-Cas9_rep1 20141031_293FT_ADAR_scr_polyAplus 20190619_LSQ10_19302)
    #srr_list=(SRR8570462)
    #for srr in ${srr_list[@]}
    #do
    # {
    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr "/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target" "_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS" >log/hg38/gatk/${srr}_run_2013_NatMethods_filter_regions_in_bed_core_flexible_`date +%Y_%m_%d`.log 2>&1
    #}
    #done

    for srr in ${srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            test $method == "gatk" && tmp_work_path=${m_path}/hg38/${method}/5-0-0vcf_filter/${srr}
            for filter_region_type in ${filter_region_type_list_part[@]}
            do
            {
                inter_name="_HaplotypeCaller_Variants_2pass_SNP"
                input_file=${tmp_work_path}/${srr}${inter_name}_${filter_region_type}.vcf
                output_file=${tmp_work_path}/${srr}${inter_name}_${filter_region_type}_deAllSNP_dbSNP_b151.vcf
                2013_NatMethods_filter_SNP gatk $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/dbSNP_b151/split_chr

                input_file=${tmp_work_path}/${srr}${inter_name}_${filter_region_type}_deAllSNP_dbSNP_b151.vcf
                output_file=${tmp_work_path}/${srr}${inter_name}_${filter_region_type}_deAllSNP_dbSNP_b151_1000genomes.vcf
                2013_NatMethods_filter_SNP gatk $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/1000genomes/split_chr

                input_file=${tmp_work_path}/${srr}${inter_name}_${filter_region_type}_deAllSNP_dbSNP_b151_1000genomes.vcf
                output_file=${tmp_work_path}/${srr}${inter_name}_${filter_region_type}_deAllSNP_dbSNP_b151_1000genomes_EVS.vcf
                2013_NatMethods_filter_SNP gatk $input_file $output_file /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/SNP/EVS/split_chr    
            }
            done
            
        }
        done
    }
    done
}
2013_NatMethods_filter_variants3(){
    method=$1
    input_file=$2
    output_file=$3
    variants_thredhold=$4
    if [ "$variants_thredhold" == "" ];then
    variants_thredhold=2.9
    fi
    test $method == "bmc" &&{
        awk '$3*$NF>='$variants_thredhold'' $input_file >$output_file
    }
    test $method == "gatk" &&{
        awk '$0!~/^#/{split($10,tmp1_a,":");dp=tmp1_a[3];split(tmp1_a[2],tmp1_c,",");v=tmp1_c[2];if (v>='$variants_thredhold'){print}}' $input_file >$output_file
    }
}
run_2013_NatMethods_filter_regions_in_bed_core_flexible(){ 
    #$method $srr ${Time_m_path}/SameBam_TwoCaller/STAR_bam "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS"
    ref_genome=hg38  
    method=$1 
    srr=$2
    m_path=$3
    tmp_inter_name=$4
    tmp_work_path_input=$5
    tmp_bam_file=$6
    need_nonAlu_v3=$7
    echo $1 $2 $3 $4
    gatk_inter_name="_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS"
    bmc_inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5"
    ################################ prepare #################################

    test $method == "bmc" && {
        result_path=4-0-0Editing_sites
        inter_name=$bmc_inter_name
        suffix="_recal"
        }
    test $method == "gatk" && { 
        result_path=5-0-0vcf_filter
        inter_name=${gatk_inter_name}
        suffix=".vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    if [ "$tmp_inter_name" != "" ];then
    inter_name=$tmp_inter_name
    fi
    if [ "$tmp_work_path_input" != "" ];then
    tmp_work_path=$tmp_work_path_input
    fi
    if [ "$need_nonAlu_v3" == "" ];then
    need_nonAlu_v3=True
    fi
    ##############################variants>=2####################################
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    test $method == "gatk" && {
        out_final_part_name=${inter_name}"_v2"
        output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
        awk '$0!~/^#/{split($10,tmp1_a,":");dp=tmp1_a[3];split(tmp1_a[2],tmp1_c,",");v=tmp1_c[2];if (v>=2){print}}' $input_result >$output_result
    }
    test $method == "bmc" && {
        out_final_part_name=${inter_name}
    }

    


    ################################split all to Alu and non-Alu#################################
    in_bed=${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed
    inter_name=$out_final_part_name
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    out_final_part_name_1=${inter_name}"_non-Alu"
    out_final_part_name_2=${inter_name}"_Alu"
    output_result_1="${tmp_work_path}/${srr}${out_final_part_name_1}${suffix}"
    output_result_2="${tmp_work_path}/${srr}${out_final_part_name_2}${suffix}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result_1"

    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result_1
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result_2 reverse

    ##############################non Alu variants>=3####################################
    if [ $need_nonAlu_v3 == "True" ];then
    inter_name=$out_final_part_name_1 
    out_final_part_name=${inter_name}"_v3" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "2013_NatMethods_filter_variants3 $method $input_result $output_result"
    2013_NatMethods_filter_variants3 $method $input_result $output_result
    else
    inter_name=$out_final_part_name_1 
    out_final_part_name=${inter_name}
    fi



    ################################ filter SimpleRepeat ################################

    in_bed="${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    inter_name=$out_final_part_name
    out_final_part_name=${inter_name}"_deSimpleRepeat" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result

    ################################ filter intronic4bp ################################

    in_bed=${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_intronic_4site.bed #${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed #"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_intronic4bp" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result

    ################################ filter Homopolymer ################################

    if [ $method == "bmc" ];then
    {
        bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_dedupped_recal_gatk4.bam
    }
    elif [ $method == "gatk" ];then
    {
        bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
    }
    fi
    if [ "$tmp_bam_file" != "" ];then
    bam_file=$tmp_bam_file
    fi
    test -e $bam_file || {
        echo "$bam_file not exists!"
        exit
    }
    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_deHomopolymer" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "2013_Natmethods_deHomopolymer $input_result $bam_file $output_result $method"
    2013_Natmethods_deHomopolymer $input_result $bam_file $output_result $method

    ################################ filter blat ################################

    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_blat" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "prepare_for_blat_2013_Natmethods_version $input_result $bam_file $output_result $method"
    prepare_for_blat_2013_Natmethods_version $input_result $bam_file $output_result $method

    ################################strandness################################

    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_annotationStrand" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "distinguish_plus_minus_main_annotation_version $method $input_result $output_result"
    distinguish_plus_minus_main_annotation_version $method $input_result $output_result

    ################################split_non-alu_to_nonAlu_nonRepeat#################################

    #_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.
    inter_name=$out_final_part_name 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    #echo $input_result
    echo "core_of_split_non_alu_to_nonAlu_nonRepeat $method $input_result #output name is auto determined by input "
    core_of_split_non_alu_to_nonAlu_nonRepeat $method $input_result #output name is auto determined by input 

    ################################reshape to avimimic################################
    #python solution
    #be_off_target_huge.py function:reshape_avimimic_filter_region_type
}

run_2013_NatMethods_filter_regions_in_bed_core(){ 
    ref_genome=hg38  
    method=$1 
    srr=$2
    m_path=$3
    
    in_bed=${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed
    test $method == "bmc" && {
        #20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal"
        }
    test $method == "gatk" && { 
        #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf";
        out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu.vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result
    in_bed="${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    test $method == "bmc" && {
        #20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_recal" #.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal"
        }
    test $method == "gatk" && { 
        #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf";
        out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu.vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result

    in_bed=${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_intronic_4site.bed #${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed #"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    test $method == "bmc" && {
        #20190619_LSQ10_19302.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal"
        out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_recal" #".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_recal" #.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_recal"
        }
    test $method == "gatk" && { 
        #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP.vcf";
        out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat.vcf" #"_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu.vcf"
        }  
    tmp_work_path=${m_path}/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result



    if [ $method == "bmc" ];then
    {
        bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
        
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_recal
        output_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
        2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
        echo -e "$method\t$srr"
        awk_base_change_percentage $method $output_filter
    }
    elif [ $method == "gatk" ];then
    {
        bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp.vcf
        output_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
        2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
        echo "$method\t$srr"
        awk_base_change_percentage $method $output_filter
    }
    fi

    ################################ filter blat ################################
    if [ $method == "bmc" ];then
    {
        bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
        
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
        out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_recal
        prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter $method
    }
    elif [ $method == "gatk" ];then
    {
        bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
        out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat.vcf
        prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter $method
    }
    fi
    ################################strandness################################
    if [ $method == "bmc" ];then
    {
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_recal
        out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_recal
        distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter
    }
    elif [ $method == "gatk" ];then
    {
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat.vcf
        out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.vcf
        distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter
    }
    fi
    ################################split_non-alu_to_nonAlu_nonRepeat#################################
    if [ $method == "bmc" ];then
    {
        bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_recal
        core_of_split_non_alu_to_nonAlu_nonRepeat $method $bmc_file #output name is auto determined by input 
    }
    elif [ $method == "gatk" ];then
    {
        bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.vcf
        core_of_split_non_alu_to_nonAlu_nonRepeat $method $bmc_file #output name is auto determined by input 
    }
    fi
    ################################reshape to avimimic################################
    #python solution
    #be_off_target_huge.py function:reshape_avimimic_filter_region_type
}



2013_Natmethods_run_core(){
    method_list=(bmc)
    ref_genome=hg38
    for srr in ${Six_testdata_srr_list[@]}
    do
    {
        #srr=${srr}_trimed
        for method in ${method_list[@]}
        do
        {
            if [ $method == "bmc" ];then
            {
                bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
                out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_annotationStrand_recal
                #distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter
                

                bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_Alu_recal
                out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_Alu_annotationStrand_recal
                distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter


                bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_annotationStrand_recal
                
                #core_of_split_non_alu_to_nonAlu_nonRepeat $method $bmc_file #output name is auto determined by input 
            }
            elif [ $method == "gatk" ];then
            {
                bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
                out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_annotationStrand.vcf
                #distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter


                bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_Alu.vcf
                out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_Alu_annotationStrand.vcf
                distinguish_plus_minus_main_annotation_version $method $bmc_file $out_blat_filter

                bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_annotationStrand.vcf
                #core_of_split_non_alu_to_nonAlu_nonRepeat $method $bmc_file #output name is auto determined by input 

            }
            fi
        }
        done
    } &
    done
    wait 
}

2013_Natmethods_non-Alu_deSimpleRepeat(){
    ref_genome=hg38
    method=$1
    srr=$2
    test $method == "bmc" && {
        result_path=4-0-0Editing_sites ;inter_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_annotationStrand_recal";out_final_part_name=".BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_annotationStrand_deSimpleRepeat_recal"
        }
    test $method == "gatk" && { 
        result_path=5-0-0vcf_filter;inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_annotationStrand.vcf";out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_annotationStrand_deSimpleRepeat.vcf"
        }  
    tmp_work_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/${method}/${result_path}/${srr}
    input_result="${tmp_work_path}/${srr}${inter_name}"
    out_aviput="${tmp_work_path}/${srr}${inter_name}_recal"
    in_bed="${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_hg38.bed"
    out_1baseBed="${tmp_work_path}/${srr}${inter_name}_1baseBed_recal"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}"
    test $method == "bmc" && ( test -e $out_aviput || awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $input_result >${out_aviput} )
    test $method == "gatk" && ( test -e $out_aviput ||${user_bin}/convert2annovar.pl -format vcf4 $input_result > ${out_aviput})
    annotate_variation.pl ${out_aviput} `dirname $in_bed` --buildver ${ref_genome}  -bedfile `basename $in_bed` -dbtype bed -regionanno -out ${out_1baseBed}
    if [ $method == "bmc" ];then
    {
        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print}}' ${out_1baseBed}.${ref_genome}_bed $input_result >${output_result}
    }
    elif [ $method == "gatk" ];then
    {
        awk 'BEGIN{FS="\t"}FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1":"$2]){print}}' ${out_1baseBed}.${ref_genome}_bed $input_result >${output_result}
    }
    fi
    rm ${out_1baseBed}.${ref_genome}_bed $out_aviput
    echo $method $srr
    awk_base_change_percentage $method ${output_result}
}
2013_Natmethods_prepare_intronic_bed(){
    #preare ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_intronic_4site.bed
    ref_genome=hg38
    test -e ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite.bed || awk -F "\t" '$3=="exon"{match($9,/transcript_id \"([^;])+\"/);geneid=substr($9,RSTART+15,RLENGTH-16);match($9,/exon_number \"([^;])+\"/);exonnum=substr($9,RSTART+13,RLENGTH-14);print $1"\t"$4-1"\t"$5"\t"geneid"_"exonnum"\t0\t"$7}' ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all.gtf >${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite.bed  

    awk 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{chromsize_array[$1]=$2}FILENAME==ARGV[2] {if ($2!="0"){print $1,$2-4,$2,$4,$5,$6}if ($3+4<=chromsize_array[$1]){print $1,$3,$3+4,$4,$5,$6}}' ${dep_path}/${ref_genome}/${ref_genome}_annotation/${ref_genome}.chromsize ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite.bed  >${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_intronic_4site.bed
    #bedtools slop -i ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite.bed -g ${dep_path}/${ref_genome}/${ref_genome}_annotation/${ref_genome}.chromsize -b -4 >${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite_minus4.bed
    #bedtools subtract -a ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite.bed -b ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_junctionsite_minus4.bed |bedtools sort -g ${dep_path}/${ref_genome}/${ref_genome}_annotation/${ref_genome}.chromsize -i /dev/stdin > ${dep_path}/${ref_genome}/${ref_genome}_junctionsite/${ref_genome}_intronic_4site.bed

}

2013_Natmethods_deHomopolymer(){
    input=$1
    #ooc=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/11.ooc
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    ref_genome_path_blat=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.2bit
    out_Homopolymer_filter=$3
    bam_file=$2
    method=$4
    test -e $out_Homopolymer_filter && rm $out_Homopolymer_filter
    #faToTwoBit $ref_genome_path  `echo $ref_genome_path|cut -d. -f1`.2bit
    n=0
    while IFS="\t" read -r line
    do
    if [ "$n" -eq "20" ];then
    {
        n=0
        wait
    }
    fi
    let n+=1
    {
        #pos_ori=chr1:187497
        #echo "`date`_token_0"
        if [ $method == "bmc" ];then
        {
            pos_ori=`echo $line|cut -d $' ' -f1` #chr1:367182 ###gatk\bmc
        }
        elif [ $method == "gatk" ];then
        {
            pos_ori=`echo $line|awk '{print $1":"$2}'` #chr1:367182 ###gatk\bmc
        }
        fi
        chr1=`echo $pos_ori|cut -d":" -f1`  
        pos_0=`echo $pos_ori|cut -d":" -f2`  
        pos_1=`expr $pos_0 - 4` 
        pos_2=`expr $pos_0 + 4` 
        fa=$(samtools1 faidx $ref_genome_path ${chr1}:${pos_1}-${pos_2}|awk 'END{print $0}')
        filter_or_not=$(echo $fa|awk '{m_nul=substr($0,5,1);re=m_nul m_nul m_nul m_nul m_nul;if ($0 ~ re){print 0}else{print 1}}')
        #echo filter_or_not,$filter_or_not
        if [ "$filter_or_not" == "1" ];then
        {
            #echo 'yyy'
            echo "$line" >>$out_Homopolymer_filter
        }
        fi
    } &
    done < "$input"
    #echo $method $srr
    #awk_base_change_percentage $method ${output_result}

}
patch_blank_tab_split(){
    m_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target
    #filter_region_type_list=(Repetitive_non-Alu Alu Nonrepetitive)
    filter_region_type_list=(Repetitive_non-Alu Nonrepetitive)
    method_list=(bmc )
    Six_testdata_srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    #all_treat_srr_list
    #all_control_srr_list
    #all_control_srr_list=(SRR8832240)
    #all_control_srr_list
    for srr in ${Six_testdata_srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            if [ $method == "bmc" ];then
            {
                bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
                
                bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_recal
                output_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
                cp $output_filter ${output_filter}_backed
                awk 'BEGIN{RS=" ";ORS="\t"}{print }' ${output_filter}_backed >${output_filter}
                output_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_blat_Nonrepetitive_recal
                #cp $output_filter ${output_filter}_backed
                #awk 'BEGIN{RS=" ";ORS="\t"}{print }' ${output_filter}_backed >${output_filter}

                output_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_blat_Repetitive_non-Alu_recal
                #cp $output_filter ${output_filter}_backed
                #awk 'BEGIN{RS=" ";ORS="\t"}{print }' ${output_filter}_backed >${output_filter}

                #2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
                #echo -e "$method\t$srr"
                #awk_base_change_percentage $method $output_filter
            }
            elif [ $method == "gatk" ];then
            {
                bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
                
                bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp.vcf
                output_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
                cp $output_filter ${output_filter}_backed
                awk 'BEGIN{RS=" ";ORS="\t"}{print }' ${output_filter}_backed >${output_filter}
                output_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_Nonrepetitive_blat.vcf
                #cp $output_filter ${output_filter}_backed
                #awk 'BEGIN{RS=" ";ORS="\t"}{print }' ${output_filter}_backed >${output_filter}
                output_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_Repetitive_non-Alu_blat.vcf
                #cp $output_filter ${output_filter}_backed
                #awk 'BEGIN{RS=" ";ORS="\t"}{print }' ${output_filter}_backed >${output_filter}

                #2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
                #echo -e "$method\t$srr"
                #awk_base_change_percentage $method $output_filter
            }
            fi
        }
        done
        #20141031_293FT_ADAR_scr_polyAplus.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_Repetitive_non-Alu_recal
        
    }
    done
}

run_2013_Natmethods_deHomopolymer(){
    m_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target
    #filter_region_type_list=(Repetitive_non-Alu Alu Nonrepetitive)
    filter_region_type_list=(Repetitive_non-Alu Nonrepetitive)
    method_list=(bmc)
    #all_treat_srr_list
    #all_control_srr_list
    #all_control_srr_list=(SRR8832240)
    for srr in ${all_control_srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            if [ $method == "bmc" ];then
            {
                bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
                
                bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_recal
                output_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
                #2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
                echo -e "$method\t$srr"
                awk_base_change_percentage $method $output_filter
            }
            elif [ $method == "gatk" ];then
            {
                bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
                
                bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp.vcf
                output_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
                #2013_Natmethods_deHomopolymer $bmc_file $bam_file $output_filter $method
                echo $method $srr
                awk_base_change_percentage $method $output_filter
            }
            fi
        }
        done
        #20141031_293FT_ADAR_scr_polyAplus.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_Repetitive_non-Alu_recal
        
    }
    done
}
tmp_interscetion_blat(){
    m_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target
    filter_region_type_list=(Repetitive_non-Alu Nonrepetitive)
    method_list=(gatk bmc)
    #srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    for srr in ${all_control_srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            if [ $method == "bmc" ];then
            {
                
                #bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
                for filter_region_type in  ${filter_region_type_list[@]}
                do
                {
                    input_blat_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_blat_${filter_region_type}_recal
                    input_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_${filter_region_type}_recal
                    output_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_${filter_region_type}_recal
                    awk 'FILENAME==ARGV[1]{array[$1]++}FILENAME==ARGV[2] && array[$1]{print}' $input_blat_file $input_file >$output_file
                    
                }
                done
            }
            
            elif [ $method == "gatk" ];then
            {
               # bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
                for filter_region_type in  ${filter_region_type_list[@]}
                do
                {
                    #20190619_LSQ10_19302_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_Nonrepetitive.vcf
                    input_blat_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_${filter_region_type}_blat.vcf
                    input_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_${filter_region_type}.vcf
                    output_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_${filter_region_type}.vcf
                    awk 'FILENAME==ARGV[1]{array[$1]++}FILENAME==ARGV[2] && array[$1]{print}' $input_blat_file $input_file >$output_file
                }
                done
            }
            fi
        }
        done
    }
    done
}
bcp_summary(){
    m_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target
}
test_run(){
    ref_genome_path_blat=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.2bit
    #gfServer start 127.0.0.1 49252 -repMatch=2253 -stepSize=5 -log=untrans.log $ref_genome_path_blat &
    
    bmc_file=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc/4-0-0Editing_sites/20141031_293FT_ADAR_scr_polyAplus/20141031_293FT_ADAR_scr_polyAplus.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_recal
    out_blat_filter=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc/4-0-0Editing_sites/20141031_293FT_ADAR_scr_polyAplus/20141031_293FT_ADAR_scr_polyAplus.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_blat_recal
    bam_file=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc/3-0-0Combine_bam/20141031_293FT_ADAR_scr_polyAplus_combine_readgroup_sort_recal_gatk4.bam
    #prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter
    m_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target
    #filter_region_type_list=(Repetitive_non-Alu Alu Nonrepetitive)
    filter_region_type_list=(Repetitive_non-Alu)
    method_list=(gatk)
    #all_treat_srr_list
    #all_control_srr_list
    srr_list=(SRR111896 SRR111897 SRR111898 SRR111899 SRR111900 SRR111895)
    for srr in ${srr_list[@]}
    do
    {
        for method in ${method_list[@]}
        do
        {
            if [ $method == "bmc" ];then
            {
                bam_file=${m_path}/hg38/bmc/3-0-0Combine_bam/${srr}_combine_readgroup_sort_recal_gatk4.bam
                for filter_region_type in  ${filter_region_type_list[@]}
                do
                {
                    #bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_${filter_region_type}_recal
                    bmc_file=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_recal
                    #out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_blat_${filter_region_type}_recal
                    out_blat_filter=${m_path}/hg38/bmc/4-0-0Editing_sites/${srr}/${srr}.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_recal
                    prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter $method
                }
                done
            }
            fi
            if [ $method == "gatk" ];then
            {
                bam_file=${m_path}/hg38/gatk/3-0-0sam_fine-tune_2pass/${srr}_dedupped_split_recal.bam
                for filter_region_type in  ${filter_region_type_list[@]}
                do
                {
                    #bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_${filter_region_type}.vcf
                    bmc_file=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer.vcf
                    #out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_${filter_region_type}_blat.vcf
                    out_blat_filter=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat.vcf
                    prepare_for_blat_2013_Natmethods_version $bmc_file $bam_file $out_blat_filter $method
                }
                done
            }
            fi
        }
        done
        #20141031_293FT_ADAR_scr_polyAplus.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_Repetitive_non-Alu_recal
    }
    done
    #gfServer stop 127.0.0.1 49252
}
2013_NatMethods_data_prepare_download(){
    #srr_list=(SRR085726 SRR036966 SRR039629 SRR039628)
    SRR_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/SRR
    ori_fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin
    all_srr=(ERR030890 ERR030882 SRR039630 SRR039632 SRR039631 SRR039633 SRR085473 SRR087416 SRR309137 SRR309133 SRR309135 SRR309138 SRR309134 SRR309136 SRR014262 SRR085474 SRR309139 SRR085471 SRR309141 SRR309143 SRR309140 SRR309142 SRR309144 SRR085725 SRR090440 SRR090441 SRR090442 SRR107727 SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900 SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306839 SRR306840 SRR306841 SRR306842 SRR306844 SRR309262 SRR085726 SRR036966 SRR039628 SRR039629)
    srr_list=(SRR085726 SRR036966 SRR039628_9add)
    all_srr_remain=(ERR030890 ERR030882 SRR039630 SRR039632 SRR039631 SRR039633 SRR085473 SRR087416 SRR309137 SRR309133 SRR309135 SRR309138 SRR309134 SRR309136 SRR014262 SRR085474 SRR309139 SRR085471 SRR309141 SRR309143 SRR309140 SRR309142 SRR309144 SRR085725 SRR090440 SRR090441 SRR090442 SRR107727 SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900 SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306839 SRR306840 SRR306841 SRR306842 SRR306844 SRR309262)
    single_end_srr_all=(SRR309139 SRR014262 SRR309136 SRR085474 SRR111899 SRR309142 SRR111900 SRR309143 SRR111902 SRR039630 SRR111906 SRR111905 SRR111903 SRR111895 SRR111898 SRR085471 SRR039632 SRR309135 SRR087416 SRR111904 SRR111897 SRR111901 SRR107727 SRR309140 SRR306839 SRR085725 SRR039633 SRR309141 SRR309134 SRR306844 SRR309138 SRR306841 SRR085473 SRR111907 SRR309133 SRR111896 SRR309144 SRR309137 SRR039631 ERR030890 SRR085726 SRR036966 SRR039629 SRR039628) #44
    single_end_srr_remain=(SRR309139 SRR014262 SRR309136 SRR085474 SRR111899 SRR309142 SRR111900 SRR309143 SRR111902 SRR039630 SRR111906 SRR111905 SRR111903 SRR111895 SRR111898 SRR085471 SRR039632 SRR309135 SRR087416 SRR111904 SRR111897 SRR111901 SRR107727 SRR309140 SRR306839 SRR085725 SRR039633 SRR309141 SRR309134 SRR306844 SRR309138 SRR306841 SRR085473 SRR111907 SRR309133 SRR111896 SRR309144 SRR309137 SRR039631) #39
    pair_end_srr_all=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306840 SRR306842 SRR309262 ERR030882) #16
    pair_end_srr_remain=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675 SRR306840 SRR306842 SRR309262) #15
    srr=SRR112600
    if [ 1 == 0 ];then
    srr=ERR030890
    ori_fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin
    wget -O ${ori_fq_path}/${srr}.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030890/ERR030890.fastq.gz
    srr=ERR030882
    wget -O ${ori_fq_path}/${srr}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030882/ERR030882_1.fastq.gz
    wget -O ${ori_fq_path}/${srr}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030882/ERR030882_2.fastq.gz
    fi

    n=0

    pair_end_srr_remain_19_7_25=(SRR090440 SRR090441 SRR090442 SRR111936 SRR111937 SRR112672 SRR112673 SRR112674 SRR112675 SRR306840 SRR306842 SRR309262) #12
    pair_end_srr_remain_19_7_26=(SRR112672 SRR112673 SRR112674 SRR112675 SRR309262) #5
    pair_end_srr_remain_19_7_26_2=(SRR112675) #5
    for srr in ${pair_end_srr_remain_19_7_26_2[@]}
    do
    if [ "$n" -eq "3" ];then
    {
        n=0
        wait
    }
    fi
    let n+=1
    {
        prefetch  -t http -O ${SRR_path} ${srr}
        fastq-dump --split-files --gzip -O ${ori_fq_path} ${SRR_path}/${srr}.sra 
    } &
    done
    exit

    for srr in ${single_end_srr_remain[@]}
    do
    {
        wget --user-agent="Mozilla/5.0" "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=${srr}&format=fastq" -cO ${ori_fq_path}/${srr}.fastq.gz 

    }
    done
    exit

    #zcat ${ori_fq_path}/SRR039628.fastq.gz ${ori_fq_path}/SRR039629.fastq.gz |gzip >${ori_fq_path}/SRR039628_9add.fastq.gz &
    for srr in ${srr_list[@]}
    do
    {
        java -jar ${user_bin}/trimmomatic-0.38.jar SE -threads  5 /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_trimed.fastq.gz LEADING:3 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36  
        fastqc -o /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/fastqc_out -d ~/tmp -q /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}.fastq.gz

    }
    done
    srr=SRR085726
    java -jar ${user_bin}/trimmomatic-0.38.jar SE -threads  5 /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_trimed.fastq.gz LEADING:3 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35 
    srr=SRR036966
    java -jar ${user_bin}/trimmomatic-0.38.jar SE -threads  5 /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_trimed.fastq.gz LEADING:3 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:32
}
patch_combine_some_data(){
    #ERR030890.ERR030882
    #SRR039628.SRR039629
    #SRR039630.SRR039631
    #SRR039632.SRR039633
    #SRR309137.SRR309138
    #SRR309133.SRR309134
    #SRR309135.SRR309136
    #SRR309139.SRR309140
    #SRR309141.SRR309142
    #SRR309143.SRR309144
    fore_part_array=(ERR030890 SRR039628 SRR039630 SRR039632 SRR309137 SRR309133 SRR309135 SRR309139 SRR309141 SRR309143)
    suff_part_array=(ERR030882 SRR039629 SRR039631 SRR039633 SRR309138 SRR309134 SRR309136 SRR309140 SRR309142 SRR309144)
    fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin
    fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed
    used_fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
    #nohup fastqc -o $out_path -d ~/tmp -q /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/*_shorted*.fastq.gz &
    for num in `seq 0 $(echo ${#fore_part_array[@]}-1|bc)`
    do
    {
        srr1=${fore_part_array[$num]}
        srr2=${suff_part_array[$num]}
        srr_combined=${srr1}add1
        if [ $num -eq 0 ];then
        {
            srr_combined=${srr2}add1
            zcat ${fq_path}/${srr1}_shorted.fastq.gz ${fq_path}/${srr2}_shorted_1.fastq.gz |gzip >${fq_path}/${srr_combined}_shorted_1.fastq.gz
            ln -s ${fq_path}/${srr2}_shorted_2.fastq.gz ${fq_path}/${srr_combined}_shorted_2.fastq.gz
            ln -s ${fq_path}/${srr_combined}_shorted_1.fastq.gz ${used_fq_path}/${srr_combined}_shorted_1.fastq.gz
            ln -s ${fq_path}/${srr_combined}_shorted_2.fastq.gz ${used_fq_path}/${srr_combined}_shorted_2.fastq.gz
        }
        else
        {
            zcat ${fq_path}/${srr1}_shorted.fastq.gz ${fq_path}/${srr2}_shorted.fastq.gz |gzip > ${fq_path}/${srr_combined}_shorted.fastq.gz
            
            ln -s ${fq_path}/${srr_combined}_shorted.fastq.gz ${used_fq_path}/${srr_combined}_shorted.fastq.gz
        }
        fi
    }
    done
}
merge_ERR030882_ERR030890_bam(){
    srr1=ERR030890
    srr2=ERR030882
    srr_combined=${srr2}add1
    samtools1 merge -f ${star_map2}/${srr}add1 ${star_map2}/${srr}_Aligned.out.sort.bam
}

test_whether_download_all(){
for srr in ${NatMethods_50_brain_single_end_srr_all[@]}
do

test -e /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/SRR/${srr}.sra && vdb-validate /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/SRR/${srr}.sra
test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}.fastq.gz || echo $srr
test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}_shorted.fastq.gz ||echo $srr no short
done


for srr in ${NatMethods_50_brain_pair_end_srr_all[@]}
do
test -e /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/SRR/${srr}.sra && vdb-validate /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/SRR/${srr}.sra
test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}_1.fastq.gz || echo $srr
test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}_shorted_1.fastq.gz ||echo $srr no short
done

}
2013_NatMethods_data_prepare_trim_try1(){
    readlen=30
    oh=0
    index_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38
    junction_index=${index_path}/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa
    if [ 1 == 0 ];then
    test -e $junction_index || {
        cd $index_path
        python /picb/rnomics1/xuew/Human/backup/hg38_genome_junction/junction_index_R100_O6/BuildIndex_hg38.py --genome hg38 --length $readlen --overhang $oh
        cat hg38_all.fa hg38_jctns_readlen${readlen}_overhang${oh}.fa > hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa
    }
    test -e ${index_path}/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa.sa || /picb/rnomics4/rotation/fuzhican/download/bwa-0.7.12/bwa index -a bwtsw ${index_path}/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa &
    test -e ${index_path}/backtrack/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa || ln -s ${index_path}/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa ${index_path}/backtrack
    test -e ${index_path}/backtrack/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa.sa || /picb/rnomics4/rotation/fuzhican/download/bwa-0.7.12/bwa index -a is ${index_path}/backtrack/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa   
    ori_fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin
    fi
    srr=SRR111897
    srr=SRR111898
    srr=SRR111899
    fastqc  -o /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/fastqc_out/2013_NatMethods_origin -d ~/tmp -q ${ori_fq_path}/${srr}.fastq.gz
        
    fq_file=""
    /picb/rnomics1/xuew/software/bwa-0.5.9/bwa aln -t 20 /picb/rnomics1/xuew/Human/backup/hg38_genome_junction/junction_index_R100_O6/hg38_genome_jctns_readlen100_overhang6.fa $2/unmapped.fq > $2/unmapped.sai

    /picb/rnomics1/xuew/software/bwa-0.5.9/bwa samse -n 4 -f $2/unmapped.sam /picb/rnomics1/xuew/Human/backup/hg38_genome_junction/junction_index_R100_O6/hg38_genome_jctns_readlen100_overhang6.fa $2/unmapped.sai $2/unmapped.fq

}
bwa_backtrack_index_junction(){
    readlen=$1
    oh=$2
    sh_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/sh
    index_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/genome_jctns_bwa_0_5_9_index
    test -d $index_path || mkdir -p $index_path
    cd $index_path
    python ${sh_path}/BuildIndex_hg38.py  --genome hg38 --length $readlen --overhang $oh
    cat /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa hg38_jctns_readlen${readlen}_overhang${oh}.fa > hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa
    /picb/rnomics1/xuew/software/bwa-0.5.9/bwa index -a bwtsw hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa
}
bwa_backtrack_mapping(){

    test_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/Non_AGTC_manual_check

    fq_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
    index_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/genome_jctns_bwa_0_5_9_index
    srr=SRR8832240
    readlen=76
    oh=1
    in_db_fasta=${index_path}/hg38_genome_jctns_readlen${readlen}_overhang${oh}.fa
    /picb/rnomics1/xuew/software/bwa-0.5.9/bwa aln -t 10 $in_db_fasta ${fq_path}/${srr}_1.fastq.gz > ${test_path}/${srr}_1_bwa_backtrack.sai &
    /picb/rnomics1/xuew/software/bwa-0.5.9/bwa aln -t 10 $in_db_fasta ${fq_path}/${srr}_2.fastq.gz > ${test_path}/${srr}_2_bwa_backtrack.sai
    wait
    /picb/rnomics1/xuew/software/bwa-0.5.9/bwa sampe -n 4 $in_db_fasta ${test_path}/${srr}_1_bwa_backtrack.sai ${test_path}/${srr}_2_bwa_backtrack.sai ${fq_path}/${srr}_1.fastq.gz  ${fq_path}/${srr}_2.fastq.gz |samtools1 view -b |samtools1 sort -o ${test_path}/${srr}_bwa_backtrack_sort.bam
    samtools1 index -@ 5 ${test_path}/${srr}_bwa_backtrack_sort.bam
}
trim_1(){
    #srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    srr_list=(SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675)
    n=0
    for srr in ${srr_list[@]}
    do
    if [ $n == 4 ];then
    n=0
    wait
    fi
    let n+=2
    {
        
        ln -s /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq

        java -jar ${user_bin}/trimmomatic-0.38.jar SE -threads  5 /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_trimed.fastq.gz LEADING:3 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36  

        ln -s /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_trimed.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
        
        fastqc -o /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/fastqc_out -d ~/tmp -q /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_trimed.fastq.gz 

    } &
    done
    wait
    run_bmc_gatk_flow_trimed
    wait
}
run_bmc_gatk_flow(){
    ref_genome=hg38 
    srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    srr_list=(SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675)
    n=0
    #trimed_or_not="_trimed"
    #sleep 3600
    for srr in ${srr_list[@]}
    do

    if [ $n -eq 2 ];then
    wait
    fi
    let n+=1
    {
        bash sh/bmc_HISAT2_BWA_flow_6_27.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/bmc/${srr}_`date +%Y_%m_%d`.log 2>&1 &
        bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/gatk/${srr}_`date +%Y_%m_%d`.log 2>&1 
        #bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh -q ${srr} -R ${ref_genome} -p >log/${ref_genome}/gatk/${srr}_`date +%Y_%m_%d`.log 2>&1 &
    } &
    done
    wait
}

run_bmc_gatk_flow_trimed(){
    ref_genome=hg38
    need_combine=(ERR030890 SRR039628 SRR039630 SRR039632 SRR309137 SRR309133 SRR309135 SRR309139 SRR309141 SRR309143 ERR030882 SRR039629 SRR039631 SRR039633 SRR309138 SRR309134 SRR309136 SRR309140 SRR309142 SRR309144)
    #srr_list=(SRR111895 SRR111896 SRR111897 SRR111898 SRR111899 SRR111900)
    #srr_list=(SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675)
    n=0
    trimed_or_not="_shorted"
    #patch_combine_some_data
    patch_combine_some_data_pair_7_31=(ERR030882)
    patch_combine_some_data_single_7_31=(SRR039628 SRR039630 SRR039632 SRR309137 SRR309133 SRR309135 SRR309139 SRR309141 SRR309143)
    patch_single_end_ClipRead_srr_list=(SRR036966 SRR085473 SRR085725 SRR085726 SRR107727 SRR111896 SRR111897 SRR111901 SRR111904 SRR111907 SRR306839 SRR306841 SRR306844)
    #for srr in ${NatMethods_50_brain_pair_end_srr_all[@]}
    for srr in ${patch_single_end_ClipRead_srr_list[@]}
    #for srr in ${patch_combine_some_data_pair_7_31[@]}
    #for srr in ${patch_combine_some_data_single_7_31[@]}
    do
    if [ $n -eq 3 ];then
    n=0
    wait
    fi
    echo ${need_combine[@]}|grep -q $srr && continue
    let n+=1
    {
        #bash sh/bmc_HISAT2_BWA_flow_6_27.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/bmc/${srr}${trimed_or_not}_`date +%Y_%m_%d`.log 2>&1 &
        #bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/gatk/${srr}${trimed_or_not}_`date +%Y_%m_%d`.log 2>&1 
        #srr=${srr}add1
        #bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -p >log/${ref_genome}/gatk/${srr}${trimed_or_not}_`date +%Y_%m_%d`.log 2>&1 
        #bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh_tmp_begin_at_STAR -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/gatk/${srr}${trimed_or_not}_`date +%Y_%m_%d`.log 2>&1 
        #bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/gatk/${srr}${trimed_or_not}_`date +%Y_%m_%d`.log 2>&1 
        bash sh/GATK_RNA_seq_STAR_flow_19_6_27_patch_8_6.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/gatk/${srr}${trimed_or_not}_patch_8_6_`date +%Y_%m_%d`.log 2>&1 
    } &
    
    done
    wait
}
tmp_patch_run_bmc_gatk_flow_trimed(){
NatMethods_50_brain_single_end_srr_not_execute_2013_filter_7_26=(SRR309139 SRR014262 SRR309136 SRR085474 SRR111899 SRR309142 SRR111900 SRR309143)
    trimed_or_not="_shorted"
    sleep 1800
for srr in ${NatMethods_50_brain_single_end_srr_not_execute_2013_filter_7_26[@]}
do
{
 #bash sh/GATK_RNA_seq_STAR_flow_19_6_27_2013_filter.sh -q ${srr}${trimed_or_not} -R ${ref_genome} -s >log/${ref_genome}/gatk/${srr}${trimed_or_not}_`date +%Y_%m_%d`.log 2>&1 
bash sh/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core gatk ${srr}${trimed_or_not}  "/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target"
}
done


}
trim_to_specified_length(){
    
    NatMethods_50_brain_srr_readlen_list=( 50 75 45 45 45 45 45 45 35 35 36 76 76 74 74 74 74 32 32 35 76 76 36 76 76 73 73 36 76 76 76 100 100 100 100 100 100 100 100 100 100 100 100 100 100 50 50 50 50 50 50 50 50 50 76 76 76 76 76 50) 
    srr_list=(SRR111901 SRR111902 SRR111903 SRR111904 SRR111905 SRR111906 SRR111907 SRR111935 SRR111936 SRR111937 SRR112600 SRR112601 SRR112672 SRR112673 SRR112674 SRR112675)
    n=0
    declare -A brain_srr_mapping_to_readlength
		for num in `seq 0 $(echo ${#NatMethods_50_brain_srr_list[@]}-1|bc)`
		do
		if [ $n == 2 ];then
		n=0
		wait
		fi
		let n+=1
		{
        read_length=${NatMethods_50_brain_srr_readlen_list[$num]}
		srr=${NatMethods_50_brain_srr_list[$num]}
		single_or_not=`echo ${NatMethods_50_brain_single_end_srr_all[@]}|grep -q $srr &&echo "y"||echo "n"`
		paired_or_not=`echo ${NatMethods_50_brain_pair_end_srr_all[@]}|grep -q $srr &&echo "y"||echo "n"`

		echo "begin" $srr
		brain_srr_mapping_to_readlength[${srr}]=$read_length
		if [ $single_or_not == "y" ];then
		{
        test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}.fastq.gz|| ln -s /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
				cutadapt -j 15 -l $read_length -o /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_shorted.fastq.gz  /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}.fastq.gz
        test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}_shorted.fastq.gz ||ln -s /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_shorted.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq

		}
		elif [ $paired_or_not == "y" ];then
		{
        pair_end_srr_remain_19_7_26=(SRR112672 SRR112673 SRR112674 SRR112675 SRR309262)
        pair_end_srr_remain_19_7_27=(SRR112675)
        #echo ${pair_end_srr_remain_19_7_27[@]}|grep -q "$srr" || continue
        echo "begin" $srr
        #exit
        test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}_1.fastq.gz|| ln -s /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}_*.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
				cutadapt -j 15 -l $read_length -o /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_shorted_1.fastq.gz  -p  /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_shorted_2.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}_1.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/origin/${srr}_2.fastq.gz

        test -e /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/${srr}_shorted_1.fastq.gz ||ln -s /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain/trimed/${srr}_shorted_*.fastq.gz /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
		}
		fi
		}
		done
		declare -p brain_srr_mapping_to_readlength



}
tmp_test_short_readlength(){
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
    fq_list=($(ls *_shorted*.fastq.gz))
    for fq in ${fq_list[@]}
    do
    {
        rl=`expr length $(zcat $fq|head|sed -n '2p' )`
        echo $fq $rl
    }
    done
}
splitPairedEndReads(){
    # @SRR112600.1.1 HWI-EAS244:1:1:0:504 length=76
    # NATCGTCATAGTCCTGATTATCGGGGGGGTGGGTGACTTTATAATGTCTCGTTTGTTTGGAGTGACGGGGTGATCT
    # +SRR112600.1.1 HWI-EAS244:1:1:0:504 length=76
    # !CB1@#######################################################################
    # @SRR112600.1.2 HWI-EAS244:1:1:0:504 length=76
    # GCAATGTGCAAAGGAATTTTGGTCACTTCCATCAAACATTAAATAGGGAAAAATACAGAGTGGTTGCACGCAGAGG
    # +SRR112600.1.2 HWI-EAS244:1:1:0:504 length=76
    # **46;9B<:<C9=###############################################################
    # @SRR112600.2.1 HWI-EAS244:1:1:0:93 length=76
    # NACCCTGCAGGGTCCAGCTGCAGGAGGGCCCTTCTCCCAGGAAAGCTTTGGGAGGAAGGGAGAAGTAGGGGAGGGC
    file_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq/2013_NatMethods_HumanBrain

    input_file=${file_path}/SRR112600.fastq.gz
    read1_file=${file_path}/SRR112600_1.fastq
    read2_file=${file_path}/SRR112600_2.fastq
    zcat $input_file|awk 'NR%8<=4{print $0 >"'${read1_file}'"}NR%8>4{print $0 >"'${read2_file}'"}' 
    pigz -p 3 $read1_file &
    pigz -p 3 $read2_file
    wait
}
bamSurgeon(){
    exit
    #cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/bamSurgeon/gatk
    ref_genome=hg38
    ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_hg38/hg38_all.fa
    tmp_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/bamSurgeon/gatk
    tmp_work_path=${tmp_path}
    var_addsnv=${tmp_work_path}/BE3_unique_chr21.addsnv
    srr=TREAT-Cas9_rep1
    bam_file=${tmp_work_path}/${srr}_dedupped_split_recal_chr21.bam
    out_bam=${tmp_work_path}/${srr}_dedupped_split_recal_chr21_out.bam
    star_map2=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/gatk/2-0-0STAR_map_2pass
    star_index2=${star_map2}/star_index_2pass_${srr}
    ##addsnv.py -v $var_addsnv -f ${bam_file} -r ${ref_genome_path} -o ${out_bam} -p 10  --picardjar /picb/rnomics4/rotation/fuzhican/bin/picard.jar --tmpdir ${tmp_work_path}/${srr}_addsnv_tmp --aligner STAR --alignopts STARrefdir:${star_index2},outFilterMultimapNmax:1,runThreadN:20,outFileNamePrefix:${star_map2}/${srr}_


}
tmp_test_avimimic_uniq(){
    avimimic_path=${m_path}/hg38/gatk/6-0-0Summary/reshape_avimimic
    cd $avimimic_path
    avifiles=($(ls $avimimic_path))
    for avifile in ${avifiles[@]}
    do
    {
        total_lines=`cat $avifile|wc -l `
        uniq_lines=`cut -f1 $avifile|sort|uniq|wc -l`
        test $total_lines -eq `expr $uniq_lines - 1`||echo $avifile ,$total_lines , $uniq_lines
    }
    done
    #result all ok!
}

git_submit(){
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/sh
    git add --ignore-removal .
    git checkout 19_6_30_branch
    git commit -m "`date +%Y_%m_%d`"
    git push -u origin 19_6_30_branch
}
git_submit_hour(){
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/sh

    git add --ignore-removal .
    git commit -m "`date +%Y_%m_%d_%H_%M`"
    git push -u origin 19_6_30_branch
}
crontab_backed(){
    20 9 * * * bash /picb/rnomics4/rotation/fuzhican/project/Comparsion_of_alingers/sh/flow_2019_07_05.sh git_submit
    20 9 * * * bash /picb/rnomics4/rotation/fuzhican/project/RNA_editing_Time_flow_fine_tune/sh/time_flow_19_7_11.sh git_submit
    20 9 * * * bash /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/sh/support_flow_19_7_8.sh git_submit
    20 * * * * /picb/rnomics4/rotation/fuzhican/bin/chmod_add.sh /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk/5-0-0vcf_filter
    20 * * * * /picb/rnomics4/rotation/fuzhican/bin/chmod_add.sh /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/bmc/4-0-0Editing_sites
    20 * * * * /picb/rnomics4/rotation/fuzhican/bin/chmod_add.sh /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/SameBam_TwoCaller/STAR_bam/hg38/bmc/4-0-0Editing_sites
    20 * * * * /picb/rnomics4/rotation/fuzhican/bin/chmod_add.sh /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/SameBam_TwoCaller/STAR_bam/hg38/gatk/5-0-0vcf_filter

}
tmp_run_SRR6191044_SRR619105(){
    run_2013_NatMethods_filter_SNP_begin_at_SNP
    tmp_run_2013_NatMethods_filter_regions_in_bed_core_flexible
}
tmp_check_multi_script(){
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    inter_name="_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_Nonrepetitive"
    all_srr=(SRR8570462 SRR8832240 SRR8570460 TREAT-Cas9_rep1 NT-Cas9_rep1 20190619_LSQ10_19302 20141031_293FT_ADAR_scr_polyAplus)
    test_path=${m_path}/test_ground/test_script_multi_sample_base_change
    combine_file_part_name=${test_path}/combine_of_2019_Nature_AND_lab_made_seven_control_Nonrepetitive
    test -d $test_path|| mkdir -p $test_path
    suffix=".vcf"
    test -e ${combine_file_part_name}${suffix}||rm ${combine_file_part_name}${suffix}
    for srr in ${all_srr[@]}
    do
    {
        tmp_work_path=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}
        
        input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
        total_num=`cat $input_result|cut -f1,2|wc -l`
        unique_num=`cat $input_result|cut -f1,2|sort|uniq|wc -l`
        #cat $input_result >> ${combine_file_part_name}${suffix}
    }
    done
    for count_num in `seq 1 ${#all_srr[@]}`
    do
    {
        awk 'BEGIN{ORS=""}{a[$1":"$2]++;ab[$1":"$2]=ab[$1":"$2]$0"\n"}END{for (i in a){if (a[i]>="'${count_num}'"){print ab[i]}}}' ${combine_file_part_name}${suffix} >${combine_file_part_name}_count${count_num}${suffix}
        #echo ${count_num} 
        #me_bcp.sh gatk ${combine_file_part_name}_count${count_num}${suffix}|awk '{a[$1]=$2}END{print a["A_G"]+a["T_C"]}'
        echo "${count_num},uniq"
        tmp_file=${combine_file_part_name}_count${count_num}_tmp${suffix}
        cat ${combine_file_part_name}_count${count_num}${suffix}|cut -f1-5|sort|uniq >$tmp_file
        me_bcp.sh gatk $tmp_file|awk '{a[$1]=$2}END{print a["A_G"]+a["T_C"]}'
        rm $tmp_file
    }
    done
}
tmp_test_known_sites_overlap(){
    m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    
    filter_region_type_list='Repetitive_non-Alu Alu Nonrepetitive' 
    all_srr=(SRR090440 SRR090441 SRR090442 SRR111935 SRR111936 SRR111937)
    test_path=${m_path}/test_ground/known_sites_overlap
    
    test -d $test_path|| mkdir -p $test_path
    suffix=".vcf"
    
    for filter_region_type in ${filter_region_type_list[@]}
    do
    {
        echo $filter_region_type
        test $filter_region_type == "Alu" &&{
            inter_name="_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_${filter_region_type}"
        }||{
            inter_name="_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_${filter_region_type}"
        }
        combine_file_part_name=${test_path}/combine_of_2013_Natmethods_Six_paireddata_${filter_region_type}
        test -e ${combine_file_part_name}${suffix}||rm ${combine_file_part_name}${suffix}
        
        for srr in ${all_srr[@]}
        do
        {
            srr=${srr}"_shorted"
            tmp_work_path=${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}
            
            input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
            total_num=`cat $input_result|cut -f1,2|wc -l`
            unique_num=`cat $input_result|cut -f1,2|sort|uniq|wc -l`
            #cat $input_result >> ${combine_file_part_name}${suffix}
        }
        done
        for count_num in `seq 1 ${#all_srr[@]}`
        do
        {
            awk 'BEGIN{ORS=""}{a[$1":"$2]++;ab[$1":"$2]=ab[$1":"$2]$0"\n"}END{for (i in a){if (a[i]>="'${count_num}'"){print ab[i]}}}' ${combine_file_part_name}${suffix} >${combine_file_part_name}_count${count_num}${suffix}
            #echo ${count_num} 
            #me_bcp.sh gatk ${combine_file_part_name}_count${count_num}${suffix}|awk '{a[$1]=$2}END{print a["A_G"]+a["T_C"]}'
            echo "${count_num}"
            tmp_file=${combine_file_part_name}_count${count_num}_tmp${suffix}
            cat ${combine_file_part_name}_count${count_num}${suffix}|cut -f1-5|sort|uniq >$tmp_file
            overlapped_file_1=${combine_file_part_name}_count${count_num}_overlapped${suffix}
            known_sites_overlap $tmp_file $overlapped_file_1
            me_bcp.sh gatk $overlapped_file_1 #|awk '{a[$1]=$2}END{print a["A_G"]+a["T_C"]}'
            rm $tmp_file &
        }
        done
    }
    done
}
known_sites_overlap(){
    known_sites_file=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/editing_known_v6_hg38_add1.txt
    test_file=$1
    ovelaped_file=$2
    method=$3
    test "$method" == "" &&{
        method=gatk
    }
    test $method == "gatk" &&{
    overlapped_num=`awk 'FILENAME==ARGV[1]{a[$1":"$2]++}FILENAME==ARGV[2]&&a[$1":"$2]' $known_sites_file $test_file |tee $ovelaped_file|wc -l`
    }
    test $method == "bmc" &&{
    overlapped_num=`awk 'FILENAME==ARGV[1]{a[$1":"$2]++}FILENAME==ARGV[2]&&a[$1]' $known_sites_file $test_file |tee $ovelaped_file|wc -l`
    }
    total_num=`cat $test_file|wc -l`
    #echo "overlapped ratio: `echo $overlapped_num/$total_num|bc`"
    test $total_num -eq 0 && echo 0 || echo `echo "scale=2;100*$overlapped_num/$total_num"|bc`%

}
Novel_sites_percentage(){
    known_sites_file=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/editing_known_v6_hg38_add1.txt
    test_file=$1
    Novel_file=$2
    method=$3
    test "$method" == "" &&{
        method=gatk
    }
    test $method == "gatk" &&{
    Novel_num=`awk 'FILENAME==ARGV[1]{a[$1":"$2]++}FILENAME==ARGV[2]&&! a[$1":"$2]' $known_sites_file $test_file |tee $Novel_file|wc -l`
    }
    test $method == "bmc" &&{
    Novel_num=`awk 'FILENAME==ARGV[1]{a[$1":"$2]++}FILENAME==ARGV[2]&&! a[$1]' $known_sites_file $test_file |tee $Novel_file|wc -l`
    }
    total_num=`cat $test_file|wc -l`
    #echo "overlapped ratio: `echo $overlapped_num/$total_num|bc`"
    echo `echo "scale=2;100*$Novel_num/$total_num"|bc`%

}
order_fastq_name(){
    a=(SRR8832240 SRR8570460 SRR8570462 SRR8570464)
    for srr in ${a[@]}
    do
    for count in `seq 1 2`
    do
    ln -s `pwd`/${srr}_${count}.fastq.gz ${srr}_R${count}.fastq.gz
    done
    done
}
check_riboMinus_data(){
    rRNA_index=/picb/rnomics4/rotation/fuzhican/download/human/hg38/sequence/RNA_45S5/RNA45S5
    p_fq_prefix_list=(20190619_LSQ10_19302 20190619_LSQ9_19301 NT-Cas9_rep1 TREAT-Cas9_rep1 SRR8832240 SRR8570460 SRR8570462 SRR8570464)

    fq_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
    out_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/check_riboMinus_data
    
    for p_fq_prefix in ${p_fq_prefix_list[@]}
    do
    out_sam_1=${out_path}/${p_fq_prefix}_rRNA_mapped.sam
    log_file=${out_path}/${p_fq_prefix}_rRNA_mapped.log
    bowtie2 -p 10 --score-min L,-16,0 --rfg 0,7 --rdg 0,7 --mp 7,7 -S $out_sam_1 -x $rRNA_index -1 ${fq_path}/${p_fq_prefix}_R1.fastq.gz -2 ${fq_path}/${p_fq_prefix}_R2.fastq.gz > ${log_file} 2>&1 
    #bowtie2 -p 5 --score-min L,-16,0 --rfg 0,7 --rdg 0,7 --mp 7,7 -S /picb/rnomics4/rotation/fuzhican/project/Comparsion_of_alingers/test_ground/check_riboMinus_data/${s_fq_prefix}_rRNA.sam -x /picb/rnomics4/rotation/fuzhican/download/human/hg38/sequence/RNA_45S5/RNA45S5 -U ${fq_path}/${s_fq_prefix}.fastq.gz > /picb/rnomics4/rotation/fuzhican/project/Comparsion_of_alingers/test_ground/check_riboMinus_data/${s_fq_prefix}_rRNA.log 2>&1 
    done
    wait
}
filter_riboMinus_data_rRNA_mapped_reads(){
    srr=20190619_LSQ9_19301
    fq_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/fastq
    tmp_work_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/check_riboMinus_data
    samtools1 view -@ 3 -f4 ${tmp_work_path}/${srr}_rRNA_mapped.bam|cut -f1 >${tmp_work_path}/${srr}_rRNA_mapped.readid
    seqtk subseq ${fq_path}/${srr}_R1.fastq.gz  ${tmp_work_path}/${srr}_rRNA_mapped.readid |gzip >${fq_path}/${srr}_derRNA_R1.fastq.gz &
    seqtk subseq ${fq_path}/${srr}_R2.fastq.gz  ${tmp_work_path}/${srr}_rRNA_mapped.readid |gzip >${fq_path}/${srr}_derRNA_R2.fastq.gz 
}
rm_sam(){
    a=($(ls *sam))
    for file1 in ${a[@]}
    do
    related_bam=`echo $file1|cut -d. -f1`.bam
    test -e  $related_bam && rm $file1 || samtools1 view -bSh@ 3 $file1|samtools1 sort -o $related_bam
    test -e $related_bam && rm $file1
    done
}
rm_bmc_internal_file(){
    need_bam=(_unique_mismatch2.sam _bwa_unique_mis6_mapq0.sam )
    bam_sam_file_list=( _hisat2_unmap.readid _unmapped_1.fastq.gz _unmapped_2.fastq.gz  _HISAT2_unmapped.bam _HISAT2_mapped.sam _HISAT2_unmapped.fastq.gz _bwa_mapped.sam _accepted_hits.nbam  _unmapped.sort.bam _bwa.header  _accepted_hits.sort.bam _combine.bam _combine_readgroup_sort.bam _accepted_hits.nsam _combine_readgroup.bam _combine_readgroup_sort.bam.bai _combine.bam.bai _unmapped.nbam _combine_readgroup_sort_recal_gatk4.grv)
    m1_path=/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/hg19/bmc
    path_list=(1-1-0HISAT_map 2-0-0BWA_map 3-0-0Combine_bam)
    #bam_sam_file_list=(_Aligned.out.bam   _Aligned.out.bam _rg_added_sorted.bam  _dedupped.bam _dedupped.bai _dedupped_split.bam _dedupped_split.bai _dedupped_split_recal.grv _output.metrics)
    for path1 in ${path_list[@]}
    do
        cd ${m1_path}/$path1
        for bam_file in ${need_bam[@]}
        do
            rm_files=($(ls *$bam_file))
            for rm_file in ${rm_files[@]}
            do
            {
                related_bam=`echo $rm_file|cut -d. -f1`.bam
                test -e $rm_file && {
                     test -e $related_bam && rm $rm_file || {
                     samtools1 view -bSh@ 3 $rm_file |samtools1 sort -o $related_bam
                     } 
                     }
                test -e $rm_file && {
                     test -e $related_bam && rm $rm_file || {
                     samtools1 view -bSh@ 3 $rm_file |samtools1 sort -o $related_bam
                     } 
                     }
            }
            done
        done

        for bam_file in ${bam_sam_file_list[@]}
        do
            rm_files=($(ls *$bam_file))
            for rm_file in ${rm_files[@]}
            do
            test -e $rm_file && rm  $rm_file
            done
        done
    done

    a=($(ls **/*.BQ20o6ES95v0_NoRecal))
    for i in ${a[@]}
    do
    pigz -p 3 $i
    done
}
change_bulk_samTobam(){
    a=($(ls **/*sam))
    for sam_file in ${a[@]}
    do
    related_bam=`echo $sam_file|cut -d. -f1`.bam
    samtools1 view -bSh@ 3 $sam_file |samtools1 sort -o $related_bam
    test -e $related_bam && rm $sam_file
    done
}
tmp_ribo_unqiue_too_low(){
    star_index1=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/hg38/star_index_hg38
    star_map1=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk/1-0-0STAR_map_1pass
    srr=20190619_LSQ10_19302_derRNA
    STAR=/picb/rnomics4/rotation/fuzhican/bin/STAR
    test -d ${star_map1}/${srr}_STARtmp_unique && rm -r ${star_map1}/${srr}_STARtmp_unique
    $STAR --genomeDir $star_index1 --readFilesIn ${dep_path}/fastq/${srr}_R1.fastq.gz ${dep_path}/fastq/${srr}_R2.fastq.gz --outTmpDir ${star_map1}/${srr}_STARtmp_unique --readFilesCommand zcat --runThreadN 22 --outFileNamePrefix ${star_map1}/${srr}_unique_ --outSAMtype BAM Unsorted --outSAMunmapped Within --outFilterMultimapNmax 1
    
    test -d ${star_map1}/${srr}_relax_length_both_STARtmp_unique && rm -r ${star_map1}/${srr}_relax_length_both_STARtmp_unique
    $STAR --genomeDir $star_index1 --readFilesIn ${dep_path}/fastq/${srr}_R1.fastq.gz ${dep_path}/fastq/${srr}_R2.fastq.gz --outTmpDir ${star_map1}/${srr}_relax_length_both_STARtmp_unique --readFilesCommand zcat --runThreadN 22 --outFileNamePrefix ${star_map1}/${srr}_relax_length_both_unique_ --outSAMtype BAM Unsorted --outSAMunmapped Within --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMultimapNmax 1

}
tmp_read_average_mapped_length(){

    #match($9,/transcript_id \"([^;])+\"/);geneid=substr($9,RSTART+15,RLENGTH-16)
    cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/hg38/gatk/2-0-0STAR_map_2pass
    samtools1 view TREAT-Cas9_rep1_Aligned.out.bam|awk -F"\t" '{readid=$1;if ($9<0){span=-1*$9}else{span=$9};match($6,/[0-9]+N/);nnum=substr($6,RSTART,RLENGTH-1);if (nnum>0){span=span-nnum};if(span<readid_span_a[readid] || !readid_span_a[readid]){readid_span_a[readid]=span}}END{for (i in readid_span_a){print i,readid_span_a[i]}}' >TREAT-Cas9_rep1_read_span
}
tmp_dl_qc(){
    test -d ${dep_path}/fastq/GEO_HEK293/|| mkdir -p ${dep_path}/fastq/GEO_HEK293/
    wget -cO ${dep_path}/fastq/GEO_HEK293/SRR6191044.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR619/004/SRR6191044/SRR6191044.fastq.gz 2>/dev/null
    wget -cO ${dep_path}/fastq/GEO_HEK293/SRR6191045.fastq.gz  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR619/005/SRR6191045/SRR6191045.fastq.gz 2>/dev/null
    wget -cO ${dep_path}/fastq/GEO_HEK293/SRR6191046.fastq.gz  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR619/006/SRR6191046/SRR6191046.fastq.gz 2>/dev/null 
    fastqc -o /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/fastqc_out -d ~/tmp -q ${dep_path}/fastq/GEO_HEK293/*

}
#tmp_run_GEO_HEK293(){
tmp_run_20190619_LSQ10_19302_derRNA(){
    
    srr_list=(SRR6191044 SRR6191045 SRR6191046)
    srr_list=(20190619_LSQ10_19302_derRNA 20190619_LSQ9_19301_derRNA)
    ref_genome=hg38
    for srr in ${srr_list[@]}
    do
    {
        test -e ${dep_path}/fastq/${srr}.fastq.gz || ln -s ${dep_path}/fastq/GEO_HEK293/${srr}.fastq.gz ${dep_path}/fastq/${srr}.fastq.gz
        #bash sh/bmc_HISAT2_BWA_flow_6_27.sh -q $srr -R ${ref_genome}  >log/${ref_genome}/bmc/${srr}_run_main_flow_`date +%Y_%m_%d`.log 2>&1 &
        bash sh/GATK_RNA_seq_STAR_flow_19_6_27.sh -q $srr -R ${ref_genome}  >log/${ref_genome}/gatk/${srr}_run_main_flow_`date +%Y_%m_%d`.log 2>&1 &

    }
    done
    wait
    
}
rerun_somedata(){
    srr_list=(20190619_LSQ9_19301_derRNA 20190619_LSQ10_19302_derRNA)
    sh_path=${ABE_m_path}/sh
    for srr in ${srr_list[@]}
    do
    bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27.sh -q $srr -R hg38 -p -m ${ABE_m_path} -l 150 >${ABE_m_path}/log/hg38/gatk/${srr}_sjdbGTFfile_mian_flow`date +%Y_%m_%d`.log 2>&1 &
    done
    wait
    srr_list=(SRR6191044 SRR6191045 SRR6191046)
    for srr in ${srr_list[@]}
    do
    bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27.sh -q $srr -R hg38 -s -m ${ABE_m_path} -l 50 >${ABE_m_path}/log/hg38/gatk/${srr}_sjdbGTFfile_mian_flow`date +%Y_%m_%d`.log 2>&1 
    done
}
run_SameBam_TwoCaller(){
    srr_list=(TREAT-Cas9_rep1 NT-Cas9_rep1  20141031_293FT_ADAR_scr_polyAplus SRR8570460 SRR8570462 SRR8832240 20190619_LSQ9_19301_derRNA 20190619_LSQ10_19302_derRNA)
    sh_path=${ABE_m_path}/sh
    
    ref_genome=hg38
    n=0
    #srr_list=(Human_Forebrain_11Week_PostConception_Female_rep0)
    #srr_list=(Human_Forebrain_11Week_PostConception_Male_rep0)
    #srr_list=(Human_Forebrain_11Week_PostConception_Male_rep0 Human_Forebrain_11Week_PostConception_Male_rep1 Human_Forebrain_11Week_PostConception_Male_rep2 Human_Forebrain_12Week_PostConception_Male_rep0 Human_Forebrain_13Week_PostConception_Female_rep0)
    for srr in ${srr_list[@]}
    do
    if [ $n == 3 ];then
    n=0
    wait
    fi
    let n+=1
    {
        echo "[`date`]BEGIN $srr"
        #m0_path=$Time_m_path
        #m_path="${m0_path}/SameBam_TwoCaller"
        bash ${sh_path}/SameBam_TwoCaller.sh -q $srr -R hg38 -s -m ${ABE_m_path}
    }&
    done
    wait
}
tmp_run_run_SameBam_TwoCaller(){
    rerun_somedata &
    run_SameBam_TwoCaller
    wait
    
}
ABE_CBE_filter(){
    all_control_srr_list=(TREAT-Cas9_rep1 NT-Cas9_rep1 SRR8570460 SRR8570462 SRR8832240)
    all_treat_srr_list=(SRR8570461 SRR8570463 SRR8570465 TREAT-BE3_rep1 NT-BE3_rep1)
    srr_list=(${all_control_srr_list[@]} ${all_treat_srr_list[@]})
    #srr_list=(TREAT-Cas9_rep1 NT-Cas9_rep1 TREAT-BE3_rep1 NT-BE3_rep1)
    #srr_list=(SRR8570461 SRR8570463 SRR8570465 SRR8570460 SRR8570462 SRR8832240)
    srr_list=(NT-BE3_rep2 NT-nCas9_rep2 EMX1-BE3_rep2 EMX1-nCas9_rep2 TREAT-BE3_rep1 NT-BE3_rep1 TREAT-Cas9_rep1 NT-Cas9_rep1)
    srr_list=(NT-BE3_rep2 NT-nCas9_rep2 EMX1-BE3_rep2 EMX1-nCas9_rep2)
    srr_list=(SRR8570463 SRR8570465)
    srr_list=(TREAT-BE3_rep1 NT-BE3_rep1)

    srr_list=(NT-BE3_rep2 NT-nCas9_rep2)
    srr_list=(NT-BE3_rep2)
    srr_list=(NT-BE3_rep2 NT-nCas9_rep2 EMX1-BE3_rep2 EMX1-nCas9_rep2)
    srr_list=(SRR111895_shorted SRR111896_shorted SRR111897_shorted SRR111898_shorted SRR111899_shorted SRR111900_shorted SRR111901_shorted SRR111902_shorted SRR111903_shorted)
    srr_list=(NT-BE3_rep2)
    srr_list=(SRR111895_shorted SRR111896_shorted SRR111897_shorted SRR111898_shorted SRR111899_shorted SRR111900_shorted SRR111901_shorted SRR111902_shorted SRR111903_shorted SRR107727_shorted)
    #srr_list=(NT-BE3_rep1)
    #srr_list=(SRR111903_shorted)
    #srr_list=(NT-nCas9_rep2_dedupped NT-Cas9_rep1_dedupped)
    
    srr_list=(SRR107727_shorted)
    srr_list=(NT-Cas9_rep1)
    ref_genome=hg38
    n=0
    for srr in ${srr_list[@]}
    do
    if [ $n == 2 ];then
    n=0
    wait
    fi
    let n+=1
    {
        echo "[`date`]BEGIN $srr"
        m0_path=$ABE_m_path
        m_path="${m0_path}"
        inter_name="_HaplotypeCaller_Variants_2pass_SNP"
        
        gatk_inter_name=$inter_name
        bam_m_path="${m0_path}/${ref_genome}/gatk"
        sam_fine_tune=${bam_m_path}/3-0-0sam_fine-tune_2pass
        from_bam=${sam_fine_tune}/${srr}_dedupped_split_recal.bam
        ####awk -F":" '$7>=10' ${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}${inter_name}.vcf >${m_path}/hg38/gatk/5-0-0vcf_filter/${srr}/${srr}${inter_name}_hits10.vcf
        #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m_path} ${inter_name}

        #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr "${m_path}" "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam True 

        inter_name=".BQ20o6ES95v2.allvariants"
        #inter_name=".BQ20o6ES95v2.allvariants_deFPS"
        bmc_inter_name=$inter_name
        bam_m_path="${m0_path}/${ref_genome}/bmc"
        combine_bam=${bam_m_path}/3-0-0Combine_bam
        from_bam=${combine_bam}/${srr}_combine_readgroup_sort_recal_gatk4.bam
        from_bam=${combine_bam}/${srr}_combine_readgroup_sort_dedupped_recal_gatk4.bam

        #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible bmc $srr ${m_path} ${inter_name}
        
        #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible bmc $srr ${m_path} "${bmc_inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS_HPB3_ER5" "" $from_bam True 
        #run_2013_NatMethods_filter_regions_in_bed_core_flexible bmc $srr ${m_path} "${bmc_inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS_HPB3_ER5" "" $from_bam True 
        bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible bmc $srr ${m_path} "${bmc_inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam True 
        
    }&
    done
    wait
}
tmp_ABE_CBE_filter_bmc_patch(){
srr=NT-BE3_rep2
srr=EMX1-BE3_rep2
#srr=EMX1-nCas9_rep2
inter_name=".BQ20o6ES95v2.allvariants"
from_bam=/data/rnomics5/wangying/Nature_2019/02_HISAT2_BWA/05_baseRecalibrator/${srr}.recal.bam
#bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible bmc $srr ${ABE_m_path} ${inter_name}
run_2013_NatMethods_filter_regions_in_bed_core_flexible bmc $srr ${ABE_m_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS_HPB3_ER5" "" $from_bam True
#run_2013_NatMethods_filter_regions_in_bed_core_flexible bmc $srr ${ABE_m_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam True
}
ABE_CBE_filter_pacth(){
    srr_list=(SRR8570461 SRR8570463 SRR8570465 SRR8570460 SRR8570462 SRR8832240)
    srr_list=(SRR8570465 SRR8570460 SRR8570462 SRR8832240)
    n=0
    for srr in ${srr_list[@]}
    do
    echo "[`date`]BEGIN $srr"
    if [ $n -eq 2 ];then
    n=0
    wait
    fi
    let n+=1
    bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27_patch_19_9_3.sh -q $srr -R hg38 -s -m ${ABE_m_path} -l 101 >log/hg38/gatk/${srr}_ABE_CBE_filter_pacth_`date +%Y_%m_%d`.log 2>&1 &
    done
}
single_run_ABE_CBE_filter_pacth(){
    ref_genome=hg38
    n=0
    srr=SRR8570462
    Time_m_path=/data/rnomics6/fuzhican/project/RNA_editing_Time_flow_fine_tune
    ABE_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
    sh_path=${ABE_m_path}/sh
    nohup memi bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27_patch_19_9_3.sh -q $srr -R hg38 -s -m ${ABE_m_path} -l 101 >log/hg38/gatk/${srr}_ABE_CBE_filter_pacth_`date +%Y_%m_%d`.log 2>&1 &
}
ABE_CBE_rep2_link_bmc(){
    ###### Tue Sep 3 19:27:12 CST 2019
    origin_path="/data/rnomics5/wangying/Nature_2019/01_star_uniq/08_HaplotypeCaller"
    origin_path="/data/rnomics5/wangying/Nature_2019/02_HISAT2_BWA/06_mpileup"
    target_path="${ABE_m_path}/hg38/bmc/4-0-0Editing_sites"
    srr_list=(NT-BE3_rep2 NT-nCas9_rep2 EMX1-BE3_rep2 EMX1-nCas9_rep2)
    bwa_index_path=${dep_path}/bwa_mem_index_${ref_genome}
    ref_genome_path=${bwa_index_path}/${ref_genome}_all.fa
    gatk4=/picb/rnomics4/rotation/fuzhican/bin/gatk4
    suffix=".BQ20o6ES95v0"

    for srr in ${srr_list[@]}
    do
    #test -d ${target_path}/${srr}||mkdir -p ${target_path}/${srr}
    ln -s ${origin_path}/${srr}${suffix} ${target_path}/${srr}/${srr}${suffix}_recal
    done
}
ABE_CBE_rep2_link_gatk(){
    ###### Tue Sep 3 19:27:12 CST 2019
    origin_path="/data/rnomics5/wangying/Nature_2019/01_star_uniq/08_HaplotypeCaller"
    #origin_path="/data/rnomics5/wangying/Nature_2019/02_HISAT2_BWA/06_mpileup"
    target_path="${ABE_m_path}/hg38/bmc/4-0-0Editing_sites"
    srr_list=(NT-BE3_rep2 NT-nCas9_rep2 EMX1-BE3_rep2 EMX1-nCas9_rep2)
    bwa_index_path=${dep_path}/bwa_mem_index_${ref_genome}
    ref_genome_path=${bwa_index_path}/${ref_genome}_all.fa
    gatk4=/picb/rnomics4/rotation/fuzhican/bin/gatk4
    suffix=".BQ20o6ES95v0"
    suffix=".vcf"
    for srr in ${srr_list[@]}
    do
    #test -d ${target_path}/${srr}||mkdir -p ${target_path}/${srr}
    #ln -s ${origin_path}/${srr}${suffix} ${target_path}/${srr}/
    {
    var_call=${ABE_m_path}/hg38/gatk/4-0-0var_calling
    vcf_filter_path=${ABE_m_path}/hg38/gatk/5-0-0vcf_filter/${srr}
    test -d ${vcf_filter_path}||mkdir -p ${vcf_filter_path}
    tag="_2pass"
    ln -s ${origin_path}/${srr}${suffix} ${var_call}/${srr}_HaplotypeCaller_Variants${tag}.vcf
    ln -s ${origin_path}/${srr}${suffix}.idx ${var_call}/${srr}_HaplotypeCaller_Variants${tag}.vcf.idx
    $gatk4 SelectVariants -select-type SNP -R ${ref_genome_path} -V ${var_call}/${srr}_HaplotypeCaller_Variants${tag}.vcf -O ${vcf_filter_path}/${srr}_HaplotypeCaller_Variants${tag}_SNP.vcf
    }
    done
}
ABE_CBE_rep2_run(){
    ###### Tue Sep 3 19:27:15 CST 2019
    srr_list=(NT-BE3_rep2 NT-nCas9_rep2 EMX1-BE3_rep2 EMX1-nCas9_rep2)
    srr_list=(EMX1-BE3_rep2 EMX1-nCas9_rep2)
    #srr_list=(NT-BE3_rep2)
    ref_genome=hg38
    n=0
    for srr in ${srr_list[@]}
    do
    echo "[`date`]BEGIN $srr"
    if [ $n -eq 2 ];then
    n=0
    wait
    fi
    let n+=1
    m0_path=$ABE_m_path
    m_path="${m0_path}"
    inter_name="_HaplotypeCaller_Variants_2pass_SNP"
    
    gatk_inter_name=$inter_name
    bam_m_path="${m0_path}/${ref_genome}/gatk"
    {
    from_bam=/data/rnomics5/wangying/Nature_2019/01_star_uniq/07_BaseRecalibrator/${srr}_dedupped_splitN_recal.bam

    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m_path} ${inter_name}

    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr "${m_path}" "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam True
    bash sh/bmc_HISAT2_BWA_flow_6_27_patch_19_9_3.sh -q $srr -R ${ref_genome} -s >log/${ref_genome}/bmc/${srr}_main_flow`date +%Y_%m_%d`.log 2>&1 
    } &
    done
    wait
}
2013_NatMethods_data_filter_bmc(){
    srr_list=(SRR111895_shorted SRR111896_shorted SRR111897_shorted SRR111898_shorted SRR111899_shorted SRR111900_shorted SRR111901_shorted SRR111902_shorted SRR111903_shorted)
    n=0
    for srr in ${srr_list[@]}
    do
    echo "[`date`]BEGIN $srr"
    if [ $n -eq 2 ];then
    n=0
    wait
    fi
    let n+=1
    bash sh/bmc_HISAT2_BWA_flow_6_27.sh -q $srr -R ${ref_genome} -s >log/${ref_genome}/bmc/${srr}_main_flow_`date +%Y_%m_%d`.log 2>&1 &
    done
}
2013_NatMethods_data_filter_pacth(){
    srr_list=(SRR8570461 SRR8570463 SRR8570465 SRR8570460 SRR8570462 SRR8832240)
    srr_list=(SRR8570465 SRR8570460 SRR8570462 SRR8832240)
    srr_list=(SRR111895_shorted SRR111896_shorted SRR111897_shorted SRR111898_shorted SRR111899_shorted SRR111900_shorted SRR111901_shorted SRR111902_shorted SRR111903_shorted)
    srr_list=(SRR107727_shorted)
    n=0
    for srr in ${srr_list[@]}
    do
    echo "[`date`]BEGIN $srr"
    if [ $n -eq 2 ];then
    n=0
    wait
    fi
    let n+=1
    bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27_end_before_calling.sh -q $srr -R hg38 -s -m ${ABE_m_path} -l 100 >log/hg38/gatk/${srr}_2013_NatMethods_data_filter_pacth_`date +%Y_%m_%d`.log 2>&1 &
    done
    wait
}
2013_NatMethods_data_filter_bmc_run(){
    n=0
    srr_list=(SRR107727_shorted SRR111895_shorted SRR111896_shorted SRR111897_shorted)
    srr_list=(SRR111898_shorted SRR111899_shorted SRR111900_shorted SRR111901_shorted SRR111902_shorted SRR111903_shorted)
    srr_list=(SRR111901_shorted SRR111902_shorted SRR111903_shorted)
    for srr in ${srr_list[@]}
    do
    echo "[`date`]BEGIN $srr"
    if [ $n -eq 3 ];then
    n=0
    wait
    fi
    let n+=1
    {
        bash ${sh_path}/bmc_HISAT2_BWA_flow_6_27.sh -q $srr -R ${ref_genome} -s >log/${ref_genome}/bmc/${srr}_main_flow_`date +%Y_%m_%d`.log 2>&1 
        echo "bash ${sh_path}/bmc_HISAT2_BWA_flow_6_27.sh -q $srr -R ${ref_genome} -s >log/${ref_genome}/bmc/${srr}_main_flow_`date +%Y_%m_%d`.log 2>&1" 
    }&
    done
    wait

}
dedupped_NT_Cas9_rep12(){
    srr_list=(NT-nCas9_rep2 NT-Cas9_rep1)

    for srr in ${srr_list[@]}
    do
        bash ${sh_path}/bmc_HISAT2_BWA_flow_6_27.sh -q $srr -R ${ref_genome} -m ${ABE_m_path} -p >log/${ref_genome}/bmc/${srr}_dedup_main_flow_`date +%Y_%m_%d`.log 2>&1 &
    done
    wait
}
gifts_test_gff3toGenePred(){
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -useName -rnaNameAttr=gene_name -geneNameAttr=transcript_id gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_useName_rnaNameAttr_geneNameAttr
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -useName -geneNameAttr=transcript_id gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_useName_geneNameAttr &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -useName -rnaNameAttr=gene_name  gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_useName_rnaNameAttr &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -useName gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_useName &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred  -geneNameAttr=transcript_id gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_geneNameAttr &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -rnaNameAttr=gene_name  gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_rnaNameAttr &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred  gencode.v31.annotation.gff3 gencode.v31.annotation.genePred &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -refseqHacks gencode.v31.annotation.gff3 gencode.v31.annotation.genePred_refseqHacks &

    in_gff3=gencode.v31.annotation.gff3
    out_genePred=gencode.v31.annotation.genePred
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred $in_gff3 ${out_genePred}_tmp1 &
    /picb/rnomics4/rotation/fuzhican/download/ucsc_utils/gff3ToGenePred -rnaNameAttr=gene_name $in_gff3 ${out_genePred}_rnaNameAttr 
    wait
    cut -f2- /picb/rnomics1/database/Human/hg38/annotation/ref_all.txt|/picb/rnomics4/rotation/fuzhican/download/ucsc_utils/genePredToGtf file stdin ref_all.gtf
}
split_strand_SameBam_TwoCaller(){
    srr=NT-Cas9_rep1

    bash ${sh_path}/SameBam_TwoCaller.sh -q $srr -R hg38 -s -m ${ABE_m_path} &
    #:q:R:b:m:spu 
    #bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27_patch_8_6.sh -q $srr -R hg38 -s -m ${ABE_m_path} -b "" &
    srr=NT-Cas9_rep1_R1
    bash ${sh_path}/SameBam_TwoCaller.sh -q $srr -R hg38 -s -m ${ABE_m_path} &
    #:q:R:b:m:spu 
    #bash ${sh_path}/GATK_RNA_seq_STAR_flow_19_6_27_patch_8_6.sh -q $srr -R hg38 -s -m ${ABE_m_path} -b "" &
    wait
}
split_strand_SameBam_TwoCaller_patch(){
    inter_name=".BQ20o6ES95v2.allvariants"
    srr=NT-Cas9_rep1
    m2_path=${ABE_m_path}/SameBam_TwoCaller/STAR_bam
    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m_path} ${inter_name}
    bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible bmc $srr ${m2_path} ${inter_name}
    bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible bmc $srr ${m2_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam True 
}
SameCaller_TwoBam_NT_Cas9_rep1(){
    srr=NT-Cas9_rep1

    bash ${sh_path}/SameBam_TwoCaller.sh -q $srr -R hg38 -s -m ${ABE_m_path}        
}
SameCaller_TwoBam_NT_Cas9_rep1_R1(){
    srr=NT-Cas9_rep1_R1

    bash ${sh_path}/SameBam_TwoCaller.sh -q $srr -R hg38 -s -m ${ABE_m_path}        
}
test_overlap_between_dbSNP_All_and_known_sites(){
    known_sites_file=${ABE_m_path}/dep_files/hg38/editing_known_v6_hg38_add1.txt
    dbSNP_file_path=${ABE_m_path}/dep_files/hg38/SNP/dbSNP_b151/split_chr
    chr_list=(`ls $dbSNP_file_path|cut -d. -f1`)
    for chr in ${chr_list[@]}
    do
    chr_known_sites_num=`cat $known_sites_file|awk '$1=="'$chr'"'|wc -l`
    chr_known_sites_in_dbSNP_num=`awk 'FILENAME==ARGV[1]{array_1[$1":"$2]}FILENAME==ARGV[2]&&array_1["chr"$1":"$2]' $known_sites_file <(zcat $dbSNP_file_path/${chr}.gz)|wc -l`
    echo -e "$chr\t$chr_known_sites_num\t$chr_known_sites_in_dbSNP_num"
    done


}
test_overlap_between_dbSNP_common_and_known_sites(){
    known_sites_file=${ABE_m_path}/dep_files/hg38/editing_known_v6_hg38_add1.txt
    dbSNP_file_file=${ABE_m_path}/dep_files/hg38/SNP/dbSNP_b151/NCBI_dbSNP_b151_common_hg38_coordinate.vcf
    known_sites_in_dbSNP_num=`awk 'FILENAME==ARGV[1]{array_1[$1":"$2]}FILENAME==ARGV[2]&&array_1[$1":"$2]' $known_sites_file $dbSNP_file_file|wc -l`
    known_sites_num=`cat $known_sites_file|wc -l`
    echo -e "known_sites_num\tknown_sites_in_dbSNP_num\tratio\n$known_sites_num\t$known_sites_in_dbSNP_num\t`echo $known_sites_in_dbSNP_num/$known_sites_num|bc`"


}

"$@"







