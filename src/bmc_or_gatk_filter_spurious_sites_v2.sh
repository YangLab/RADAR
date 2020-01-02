#######################################################################
# File Name: bmc_or_gatk_filter_spurious_sites_v2.sh
# Author: Zhi-Can Fu
# mail: yyqin90@gmail.com
# Created time: 2015-09-03 00:08:03
# Last modified: 2015-10-27 21:56:45
#!/bin/bash
echo `date` dir:`pwd` command:$0 $* >>$0.arg


#set -eo pipefail
######################################usage############################################################
#set -eo pipefail
######################################usage############################################################

#method="bmc"; inter_name=".BQ20o6ES95v2.allvariants"
#or 
#method="gatk"; inter_name="_HaplotypeCaller_Variants_BQ20_MAPQ0_2pass_SNP"

##filter SNP: 
#bash  bmc_or_gatk_filter_spurious_sites.sh run_2013_NatMethods_filter_SNP_flexible ${method} $outname ${ABE_outdir} ${inter_name}
##filter_spurious_sites:
#bash bmc_or_gatk_filter_spurious_sites.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible ${method} $outname ${ABE_outdir} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS"

######################################usage############################################################
MYDIR=`dirname $0`
dir_tools=${MYDIR}/../tools
dir_scripts=${MYDIR}/../src
chmod +x ${dir_scripts}/annotate_variation.pl
chmod +x ${dir_scripts}/convert2annovar.pl
chmod +x ${dir_scripts}/pslScore.pl
chmod +x ${dir_tools}/random
chmod +x ${dir_tools}/gfServer
chmod +x ${dir_tools}/gfClient
chmod +x ${dir_tools}/netstat



run_2013_NatMethods_filter_SNP_flexible(){

        GETOPT_ARGS=`getopt  -o m:n:o:i:g:  -al call_var_method:,outname:,outdir:,inter_name:,genome_build_version:,SNP_dbSNP_divided_by_chromosome:,SNP_1000Genome_divided_by_chromosome:,SNP_EVS_divided_by_chromosome: -- "$@"`
        echo ${GETOPT_ARGS}
        eval set -- "$GETOPT_ARGS"
        while [ -n "$1" ]
        do
                case "$1" in
                        -m | --call_var_method | -call_var_method) method=$2; shift 2;;
                        -n | --outname | -outname) outname=$2; shift 2;;
                        -o | --outdir | -outdir) outdir=$2; shift 2;;
                        -i | --inter_name | -inter_name) tmp_inter_name=$2; shift 2;;
                        -g | --genome_build_version | -genome_build_version) ref_genome=$2; shift 2;;
                        --SNP_dbSNP_divided_by_chromosome | -SNP_dbSNP_divided_by_chromosome) SNP_dbSNP_divided_by_chromosome=$2; shift 2;;
                        --SNP_1000Genome_divided_by_chromosome | -SNP_1000Genome_divided_by_chromosome) SNP_1000Genome_divided_by_chromosome=$2; shift 2;;
                        --SNP_EVS_divided_by_chromosome | -SNP_EVS_divided_by_chromosome) SNP_EVS_divided_by_chromosome=$2; shift 2;;
                        --) break ;;
                        *) echo $1,$2,$show_usage; break ;;
                esac
        done
<<COMMENT
    while getopts :m:n:o:i:g: ARGS  
    do  
    case $ARGS in   
        m)  
            method=$OPTARG
            ;;  
        n)  
            outname=$OPTARG  
            ;;
        o)  
            outdir=$OPTARG
            ;;  
        i)  
            tmp_inter_name=$OPTARG
            ;;  
        g)  
            ref_genome=$OPTARG
            ;;  
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
COMMENT
    
    #outdir=$outdir
    
    test $method == "gatk" && {
        suffix=".vcf"
        tmp_work_path=${outdir}/${outname}_tmp/${method}/5-0-0vcf_filter/${outname}           
        inter_name="_HaplotypeCaller_Variants_2pass_SNP"                
    }
    test $method == "bmc" && {
        suffix="_recal"
        tmp_work_path=${outdir}/${outname}_tmp/${method}/4-0-0Editing_sites/${outname}           
        inter_name=".BQ20o6ES95v2.allvariants_deFPS"
        
    }
    tmp_a="aa${tmp_inter_name}aa"
    if [ "$tmp_a" != "aaaa" ];then
    inter_name=$tmp_inter_name
    fi
    input_file=${tmp_work_path}/${outname}${inter_name}${suffix}
    #ln -s ${tmp_work_path}/${outname}.BQ20o6ES95v2.allvariants_recal_deFPS $input_file
    output_file=${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151${suffix}
    test -e $input_file ||{
        echo "$input_file do not exists,skip it!"
#        continue
    }
    test -e ${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix} && tmp_a=`head ${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix}`

    test "$tmp_a" != ""  && {
        echo "${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix} exists,skip it!"
#        continue
    }


    2013_NatMethods_filter_SNP $method $input_file $output_file ${SNP_dbSNP_divided_by_chromosome}
	#${dep_path}/${ref_genome}/SNP/dbSNP_b151/split_chr

    input_file=${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151${suffix}
    output_file=${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151_1000genomes${suffix}
    2013_NatMethods_filter_SNP $method $input_file $output_file  ${SNP_1000Genome_divided_by_chromosome}
	#${dep_path}/${ref_genome}/SNP/1000genomes/split_chr

    input_file=${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151_1000genomes${suffix}
    output_file=${tmp_work_path}/${outname}${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS${suffix}
    2013_NatMethods_filter_SNP $method $input_file $output_file ${SNP_EVS_divided_by_chromosome}
	#${dep_path}/${ref_genome}/SNP/EVS/split_chr
}

2013_NatMethods_filter_SNP(){
    #download cat 20170504_GRCh38_positions_manifest.txt|grep '_sites'|cut -f1|xargs -I{} wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/\{\}
    #link
    ref_genome=${ref_genome}  
    method=$1 
    input_file=$2
    output_file=$3
    SNP_path=$4 #/data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/dep_files/${ref_genome}/SNP/dbSNP_b151/split_chr
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

run_2013_NatMethods_filter_regions_in_bed_core_flexible(){ 
    #$method $outname ${Time_outdir}/SameBam_TwoCaller/STAR_bam "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS"
    ###### Tue Oct 8 19:47:21 CST 2019
    #UPDATE:
    #1.strandness Alu

	GETOPT_ARGS=`getopt  -o m:n:o:i:t:b:g:v:  -al call_var_method:,outname:,outdir:,inter_name:,bam_file:,genome_build_version:,genome_fasta:,genome_index_blat:,annotation_Alu:,annotation_Repetitive_non_Alu:,annotation_All_repetitive:,annotation_RepeatMasker_simple_repeats:,annotation_intronic_4site:,annotation_gene_transcribed_strands:,stranded: -- "$@"`
	echo ${GETOPT_ARGS}
        eval set -- "$GETOPT_ARGS"
        while [ -n "$1" ]
        do
                case "$1" in
			-v) need_nonAlu_v3=$2; shift 2;;
                        -m | --call_var_method | -call_var_method) method=$2; shift 2;;
                        -n | --outname | -outname) outname=$2; shift 2;;
                        -o | --outdir | -outdir) outdir=$2; shift 2;;
                        -i | --inter_name | -inter_name) tmp_inter_name=$2; shift 2;;
                        -t) tmp_work_path_input=$2; shift 2;;
                        -b | --bam_file | -bam_file) tmp_bam_file=$2; shift 2;;
			--stranded | -stranded) stranded=$2; shift 2;;
                        -g | --genome_build_version | -genome_build_version) ref_genome=$2; shift 2;;
                        --genome_index_blat | -genome_index_blat) genome_index_blat=$2; shift 2;;
                        --genome_fasta | -genome_fasta) genome_fasta=$2; shift 2;;
                        --annotation_Alu | -annotation_Alu) annotation_Alu=$2; shift 2;;
                        --annotation_Repetitive_non_Alu |annotation_Repetitive_non_Alu) annotation_Repetitive_non_Alu=$2; shift 2;;
                        --annotation_All_repetitive | -annotation_All_repetitive) annotation_All_repetitive=$2; shift 2;;
                        --annotation_RepeatMasker_simple_repeats | -annotation_RepeatMasker_simple_repeats) annotation_RepeatMasker_simple_repeats=$2; shift 2;;
                        --annotation_intronic_4site | -annotation_intronic_4site) annotation_intronic_4site=$2; shift 2;;
                        --annotation_gene_transcribed_strands | -annotation_gene_transcribed_strands) annotation_gene_transcribed_strands=$2; shift 2;;
                        --) break ;;
                        *) echo $1,$2,$show_usage; break ;;
                esac
        done
<<COMMENT
    while getopts :v:b:t:i:m:o:n:g: ARGS  
    do  
    case $ARGS in   
        v)  
            need_nonAlu_v3=$OPTARG
            ;;
        b)  
            tmp_bam_file=$OPTARG
            ;; 
        t)  
            tmp_work_path_input=$OPTARG
            ;;  
        i)  
            tmp_inter_name=$OPTARG  
            ;;
        m)  
            method=$OPTARG
            ;;  
        o)  
            outdir=$OPTARG
            ;;  
        n)  
            outname=$OPTARG
            ;; 
        g)  
            ref_genome=$OPTARG
            ;;  
        
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
COMMENT

    #outdir=$outdir


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
    tmp_a="aa${tmp_inter_name}aa"
    if [ "$tmp_a" != "aaaa" ];then
    inter_name=$tmp_inter_name
    fi

    tmp_work_path=${outdir}/${outname}_tmp/${method}/${result_path}/${outname}
    tmp_a="aa${tmp_work_path_input}aa"
    if [ "$tmp_a" != "aaaa" ];then
    tmp_work_path=$tmp_work_path_input
    fi
    tmp_a="aa${need_nonAlu_v3}aa"
    if [ "$tmp_a" == "aaaa" ];then
    need_nonAlu_v3=True
    fi
    ##############################variants>=2####################################
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    test $method == "gatk" && {
        out_final_part_name=${inter_name}"_v2"
        output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
        awk '$0!~/^#/{split($10,tmp1_a,":");dp=tmp1_a[3];split(tmp1_a[2],tmp1_c,",");v=tmp1_c[2];if (v>=2){print}}' $input_result >$output_result
    }
    test $method == "bmc" && {
        out_final_part_name=${inter_name}
    }

    ################################split all to Alu and non-Alu#################################
    in_bed=${annotation_Alu}
	#${dep_path}/${ref_genome}/Alu.bed
    inter_name=$out_final_part_name
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    out_final_part_name_1=${inter_name}"_non-Alu"
    out_final_part_name_2=${inter_name}"_Alu"
    output_result_1="${tmp_work_path}/${outname}${out_final_part_name_1}${suffix}"
    output_result_2="${tmp_work_path}/${outname}${out_final_part_name_2}${suffix}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result_1"

    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result_1
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result_2 reverse

    ################################strandness Alu################################

    inter_name=$out_final_part_name_2 
    out_final_part_name=${inter_name}"_annotationStrand" 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    if [[  "$stranded" == true   ]];then
    	echo "distinguish_plus_minus_main_annotation_version $method $input_result $output_result"
    	distinguish_plus_minus_main_annotation_version $method $input_result $output_result ${annotation_gene_transcribed_strands}
    else
	cp ${input_result} ${output_result}
    fi

    Alu_final_result="${output_result}"

    ##############################non Alu variants>=3####################################
    if [ $need_nonAlu_v3 == "True" ];then
    inter_name=$out_final_part_name_1 
    out_final_part_name=${inter_name}"_v3" 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    echo "2013_NatMethods_filter_variants3 $method $input_result $output_result"
    2013_NatMethods_filter_variants3 $method $input_result $output_result
    else
    inter_name=$out_final_part_name_1 
    out_final_part_name=${inter_name}
    fi



    ################################ filter SimpleRepeat ################################

    in_bed=${annotation_RepeatMasker_simple_repeats}
	#"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_${ref_genome}.bed"
    inter_name=$out_final_part_name
    out_final_part_name=${inter_name}"_deSimpleRepeat" 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result

    ################################ filter intronic4bp ################################

    in_bed=${annotation_intronic_4site}
	#${dep_path}/${ref_genome}/${ref_genome}_intronic_4site.bed #${dep_path}/${ref_genome}/Alu_str_${ref_genome}.bed #"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_${ref_genome}.bed"
    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_intronic4bp" 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    echo "2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result"
    2013_NatMethods_filter_regions_in_bed $method $in_bed $input_result $output_result

    ################################ filter Homopolymer ################################

    if [ $method == "bmc" ];then
    {
        bam_file=${outdir}/${outname}_tmp/bmc/3-0-0Combine_bam/${outname}_combine_readgroup_sort_dedupped_recal_gatk4.bam
    }
    elif [ $method == "gatk" ];then
    {
        bam_file=${outdir}/${outname}_tmp/gatk/3-0-0sam_fine-tune_2pass/${outname}_dedupped_split_recal.bam
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
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    echo "2013_Natmethods_deHomopolymer $input_result $bam_file $output_result $method"
    2013_Natmethods_deHomopolymer $input_result $bam_file $output_result $method ${genome_fasta}

    ################################ filter blat ################################

    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_blat" 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    echo "prepare_for_blat_2013_Natmethods_version $input_result $bam_file $output_result $method"
    prepare_for_blat_2013_Natmethods_version $input_result $bam_file $output_result $method ${genome_fasta} ${genome_index_blat} ${ref_genome}

    ################################strandness Non_Alu################################

    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_annotationStrand" 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${outname}${out_final_part_name}${suffix}"
    if [[  "$stranded" == true   ]];then
    	echo "distinguish_plus_minus_main_annotation_version $method $input_result $output_result"
    	distinguish_plus_minus_main_annotation_version $method $input_result $output_result ${annotation_gene_transcribed_strands}
    else
        cp ${input_result} ${output_result}
    fi    
    ################################split_non-alu_to_nonAlu_nonRepeat#################################

    #_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.
    inter_name=$out_final_part_name 
    input_result="${tmp_work_path}/${outname}${inter_name}${suffix}"
    #echo $input_result
    echo "core_of_split_non_alu_to_nonAlu_nonRepeat $method $input_result #output name is auto determined by input "
    core_of_split_non_alu_to_nonAlu_nonRepeat $method $input_result ${annotation_Alu} ${annotation_Repetitive_non_Alu} ${annotation_All_repetitive} 

    ################################link_and_combine_the_final_output################################

    Nonrepetitive_final_result="${tmp_work_path}/${outname}${inter_name}_Nonrepetitive${suffix}"
    Repetitive_non_Alu_final_result="${tmp_work_path}/${outname}${inter_name}_Repetitive_non-Alu${suffix}"
    
    ln_target=`dirname $Nonrepetitive_final_result`/${outname}_recommended_final_result_Nonrepetitive.vcf
    test -h $ln_target && rm $ln_target
    ln -s $(basename $Nonrepetitive_final_result) $ln_target

    ln_target=`dirname $Repetitive_non_Alu_final_result`/${outname}_recommended_final_result_Repetitive_non-Alu.vcf
    test -h $ln_target && rm $ln_target
    ln -s $(basename $Repetitive_non_Alu_final_result) $ln_target

    ln_target=`dirname $Alu_final_result`/${outname}_recommended_final_result_Alu.vcf
    All_combined_result="${tmp_work_path}/${outname}_recommended_final_result_All.vcf"
    test -h $ln_target && rm $ln_target
    ln -s $(basename $Alu_final_result) $ln_target
    All_combined_result="${tmp_work_path}/${outname}_recommended_final_result_All.vcf"
    cat $Alu_final_result $Repetitive_non_Alu_final_result $Nonrepetitive_final_result >$All_combined_result
    rm ${tmp_work_path}/*log


    ################################reshape to avimimic################################
    #python solution
    #be_off_target_huge.py function:reshape_avimimic_filter_region_type
}

core_of_split_non_alu_to_nonAlu_nonRepeat(){
    ref_genome=${ref_genome}
    method=$1
    input_file=$2
    annotation_Alu=${3}
    annotation_Repetitive_non_Alu=${4}
    annotation_All_repetitive=${5}
    deSNPtag=de_commonSNP

    filter_region_type_list=(Repetitive_non-Alu Alu)
    if [ $method == "bmc" ];then
    {
        bmc_result_file=$input_file
        editing_sites=`dirname $bmc_result_file`
        outname=`basename $editing_sites`
        tag=_recal
        avi_input=${editing_sites}/${outname}.BQ20o6ES95v2.allvariants_${deSNPtag}_HPB3_ER5.avinput${tag}
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $bmc_result_file >${avi_input}
<<COMMENT
        for filter_region_type in ${filter_region_type_list[@]}
        do
        {
            out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}${tag}
            out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
            ${dir_scripts}/annotate_variation.pl ${avi_input} ${dep_path}/${ref_genome} --buildver ${ref_genome}  -bedfile ${filter_region_type}.bed -dbtype bed -regionanno -out ${out_1baseBed}
            awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >${out_resultFile}
            rm ${out_1baseBed}.${ref_genome}_bed 
        }
        done
COMMENT
	## Repetitive_non-Alu
        filter_region_type=Repetitive_non-Alu
        annotation_file_path=`dirname $annotation_Repetitive_non_Alu`
        annotation_file_name=`basename $annotation_Repetitive_non_Alu`
        out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}${tag}
        out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
        ${dir_scripts}/annotate_variation.pl ${avi_input} ${annotation_file_path} --buildver ${ref_genome}  -bedfile ${annotation_file_name} -dbtype bed -regionanno -out ${out_1baseBed}
        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >$out_resultFile
        rm ${out_1baseBed}.${ref_genome}_bed $avi_input	

	## Alu
        filter_region_type=Alu
	annotation_file_path=`dirname $annotation_Alu`
	annotation_file_name=`basename $annotation_Alu`
        out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}${tag}
        out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
        ${dir_scripts}/annotate_variation.pl ${avi_input} ${annotation_file_path} --buildver ${ref_genome}  -bedfile ${annotation_file_name} -dbtype bed -regionanno -out ${out_1baseBed}
        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >$out_resultFile
        rm ${out_1baseBed}.${ref_genome}_bed $avi_input

	## non-repetitive
        filter_region_type=All_repetitive
        annotation_file_path=`dirname $annotation_All_repetitive`
        annotation_file_name=`basename $annotation_All_repetitive`
        out_resultFile=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_Nonrepetitive${tag}
        out_1baseBed=$(echo $bmc_result_file|awk -F "${tag}" '{print $1}')_${filter_region_type}_contain${tag}
        ${dir_scripts}/annotate_variation.pl ${avi_input} ${annotation_file_path} --buildver ${ref_genome}  -bedfile ${annotation_file_name} -dbtype bed -regionanno -out ${out_1baseBed}
        awk 'FILENAME==ARGV[1]{a[$3":"$4]++}FILENAME==ARGV[2]{if (! a[$1]){print $0}}' ${out_1baseBed}.${ref_genome}_bed $bmc_result_file >$out_resultFile
        rm ${out_1baseBed}.${ref_genome}_bed $avi_input
    }
    elif [ $method == "gatk" ];then
    {
        #convert2annovar=/picb/rnomics4/rotation/fuzhican/bin/convert2annovar.pl
        #annotate_variation=/picb/rnomics4/rotation/fuzhican/bin/annotate_variation.pl
        convert2annovar=${dir_scripts}/convert2annovar.pl
        annotate_variation=${dir_scripts}/annotate_variation.pl
        vcf_name=$input_file
        avinput_name=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}.avinput`
        $convert2annovar -format vcf4 ${vcf_name} > ${avinput_name}
<<COMMENT
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
COMMENT
	## Repetitive_non-Alu
        filter_region_type=Repetitive_non-Alu
        annotation_file_path=`dirname $annotation_Repetitive_non_Alu`
        annotation_file_name=`basename $annotation_Repetitive_non_Alu` 
        alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
        alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}.vcf`
        $annotate_variation ${avinput_name} ${annotation_file_path} --buildver ${ref_genome}  -bedfile ${annotation_file_name} -dbtype bed -regionanno -out  ${alu_contain_bed}
        awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if ( a[$1"_"$2]){print $0}}'  ${alu_contain_bed}.${ref_genome}_bed  ${vcf_name} > ${alu_filtered_vcf}
        rm  ${alu_contain_bed}.${ref_genome}_bed

	## Alu
        filter_region_type=Alu
        annotation_file_path=`dirname $annotation_Alu`
        annotation_file_name=`basename $annotation_Alu`
        alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
        alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}.vcf`
        $annotate_variation ${avinput_name} ${annotation_file_path} --buildver ${ref_genome}  -bedfile ${annotation_file_name} -dbtype bed -regionanno -out  ${alu_contain_bed}
        awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if ( a[$1"_"$2]){print $0}}'  ${alu_contain_bed}.${ref_genome}_bed  ${vcf_name} > ${alu_filtered_vcf}
        rm  ${alu_contain_bed}.${ref_genome}_bed

	## non-repetitive
        filter_region_type=All_repetitive
        annotation_file_path=`dirname $annotation_All_repetitive`
        annotation_file_name=`basename $annotation_All_repetitive`
        alu_contain_bed=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_${filter_region_type}_contain.bed`
        alu_filtered_vcf=`echo $vcf_name|cut -d. -f1 |xargs -I{} echo {}_Nonrepetitive.vcf`
        $annotate_variation ${avinput_name} ${annotation_file_path} --buildver ${ref_genome}  -bedfile ${annotation_file_name} -dbtype bed -regionanno -out  ${alu_contain_bed}
        awk 'FILENAME==ARGV[1]{a[$3"_"$4]++}FILENAME==ARGV[2]{if (! a[$1"_"$2]){print $0}}' ${alu_contain_bed}.${ref_genome}_bed  ${vcf_name} > ${alu_filtered_vcf}
        rm  ${alu_contain_bed}.${ref_genome}_bed ${avinput_name}
    }
    fi
}

2013_NatMethods_filter_regions_in_bed(){
    ###this function filtering sites localed in regions specified within in_bed.
    #in detail, this function can be use to filter Alu, Simple Repeat, intronic.
    ref_genome=${ref_genome}
    method=$1 #bmc or gatk
    in_bed=$2 #"${dep_path}/${ref_genome}/UCSC_RepeatMask_SimpleRepeats_${ref_genome}.bed"
    input_result=$3 #absolute path of input result file
    output_result=$4 #absolute path of output result file
    filter_or_not=$5
    tmp_work_path=`dirname $input_result`
    outname=`basename $tmp_work_path` #####only suit for result_file deposited at director named by outname
#    inter_name=`echo $input_result|awk -F "$outname" '{print $3}'`
    inter_name=`echo $input_result|awk -F "/$outname/$outname" '{print $NF}'`
    out_aviput="${tmp_work_path}/${outname}${inter_name}_recal" #temp internal file
    out_1baseBed="${tmp_work_path}/${outname}${inter_name}_1baseBed_recal" #temp internal file
    test $method == "bmc" && ( test -s $out_aviput || awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $input_result >${out_aviput} )
    test $method == "gatk" && ( test -s $out_aviput ||${dir_scripts}/convert2annovar.pl -format vcf4 $input_result > ${out_aviput})
    ${dir_scripts}/annotate_variation.pl ${out_aviput} `dirname $in_bed` --buildver ${ref_genome}  -bedfile `basename $in_bed` -dbtype bed -regionanno -out ${out_1baseBed}
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
    echo $method $outname
    awk_base_change_percentage $method ${output_result}
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

2013_Natmethods_deHomopolymer(){
    input=$1
    ref_genome_path=$5
	#/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_${ref_genome}/${ref_genome}_all.fa
    out_Homopolymer_filter=$3
    bam_file=$2
    method=$4
    test -e $out_Homopolymer_filter && rm $out_Homopolymer_filter
    #faToTwoBit $ref_genome_path  `echo $ref_genome_path|cut -d. -f1`.2bit
    n=0
    while IFS="\t" read -r line
    do
    if [ "$n" -eq "10" ];then
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
        fa=$(samtools faidx $ref_genome_path ${chr1}:${pos_1}-${pos_2}|awk 'END{print $0}')
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
    #echo $method $outname
    #awk_base_change_percentage $method ${output_result}

}
prepare_for_blat_2013_Natmethods_version(){
    ##Preparing
    ref_genome_path_blat=$6
    ref_genome=$7

	#/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_${ref_genome}/${ref_genome}_all.2bit
    ps -fu $(id -un)|grep  -ve "grep" |grep  -e "gfServer" |grep -qe "${ref_genome}" &&{
        Port=`ps -fu $(id -un)|grep  -ve "grep" -e "awk"|grep  -e "gfServer" |awk -F"start" '{print $NF}'|cut -d" " -f 3`
        echo "gfServer(${ref_genome}) already start! The Port used for gfServer is ${Port}"
    }|| {
        Port=$(${dir_tools}/random /49152..65000/)

        ${dir_tools}/netstat -anp 2>/dev/null|grep -q $Port
        retval=$?

        while [ $retval -eq 0 ];do
        Port=$(${dir_tools}/random /49152..65000/)
        ${dir_tools}/netstat -anp 2>/dev/null|grep -q $Port
        retval=$?
        done

        echo "The Port used for gfServer(${ref_genome}) is $Port!"
        test -e ./gfServer_untrans.log && rm ./gfServer_untrans.log
        ${dir_tools}/gfServer start 127.0.0.1 $Port -repMatch=2253 -stepSize=5 -log=./gfServer_untrans.log $ref_genome_path_blat &
        sleep 30
        cat ./gfServer_untrans.log|grep -q "Server ready for queries"
        retval=$?
        while [ $retval -ne 0 ];do
        sleep 60
        cat ./gfServer_untrans.log|grep -q "Server ready for queries"
        retval=$?
        done
        
    }   ###### Fri Nov 22 16:00:41 CST 2019 fzc; automatically set $Port!
    
    #blat:Standalone BLAT v. 36x2
    #blat  -makeOoc=11.ooc /picb/rnomics1/database/Human/${ref_genome}/genome/${ref_genome}_all.fa  /picb/rnomics1/database/Human/${ref_genome}/genome/${ref_genome}_all.fa out.psl
    input=$1
    input_basename=`basename $input`
    #ooc=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_${ref_genome}/11.ooc
    ref_genome_path=$5
	#/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_${ref_genome}/${ref_genome}_all.fa

    
    #tmp_path=~/tmp/blat_filter
    tmp_path=`dirname $input | awk '{print $0"/tmp/blat_filter"}' `
    
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
        samtools mpileup -Q 20 -O -f $ref_genome_path --output-QNAME -r $chr -l ${tmp_path}/${input_basename}/split_chr/$chr $bam_file > ${tmp_path}/${input_basename}/mpilup_output/${chr}_mpileup_output &
    }
    done

echo "samtools end"


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
        output_fa=${tmp_path}/${input_basename}/tmp_${pos_ori}_$(${dir_tools}/random).fa
        out_client_psl=${tmp_path}/${input_basename}/tmp_out_client_${pos_ori}_$(${dir_tools}/random).psl
        out_with_score_psl=${tmp_path}/${input_basename}/tmp_out_with_score_${pos_ori}_$(${dir_tools}/random).psl
        #echo "pileup_extract_readseq  $bam_file ${pos_ori}-${pos_0}"
        #echo "`date`_token_1"
        pileup_extract_readseq  $bam_file ${pos_ori}-${pos_0} $input_basename ${tmp_path} >$output_fa 
        #echo "`date`_token_2"
        ${dir_tools}/gfClient 127.0.0.1 $Port ""  -nohead $output_fa $out_client_psl >/dev/null
        #echo "`date`_token_3"
        ${dir_scripts}/pslScore.pl $out_client_psl |awk 'BEGIN{OFS="\t"}{sub(/:[0-9]+-[0-9]+/,"",$4);print}'|sort  -k4,4 -k5nr,5 >$out_with_score_psl
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
distinguish_plus_minus_main_annotation_version(){
    #cd /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/test_ground/distinguish_strand_annotation
    ref_genome=${ref_genome}
    method=$1
    local input_file=$2
    local out_final_part_name=$3
    local in_bed=$4
    local tmp_work_path=`dirname $input_file`
    local input_file_name=`basename $input_file`
    local outname=`basename $tmp_work_path`
    test $method == "bmc" && {
        result_path=4-0-0Editing_sites 
        local inter_name=`echo $input_file_name|awk -F ".BQ20o6ES95v2" '{print $2}'`
        }
    test $method == "gatk" && { 
        result_path=5-0-0vcf_filter
        local inter_name=`echo $input_file_name|awk -F "_Haplo" '{print $2}'`
        }   
    #inter_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10.vcf";out_final_part_name="_HaplotypeCaller_Variants_2pass_SNP_de_commonSNP_hits10_annotationStrand.vcf"
    #/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/${ref_genome}/bmc/4-0-0Editing_sites/SRR111896_trimed/SRR111896_trimed.BQ20o6ES95v2.allvariants_de_commonSNP_HPB3_ER5_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_annotationStrand_recal

    local input_result=$input_file
    local out_aviput="${tmp_work_path}/${outname}${inter_name}_recal"
    #in_bed="${dep_path}/${ref_genome}/ref_all_6.bed"
    out_bed_minus=${in_bed}.minus.bed
	#"${dep_path}/${ref_genome}/ref_all_6_minus.bed"
    out_bed_plus=${in_bed}.plus.bed
	#"${dep_path}/${ref_genome}/ref_all_6_plus.bed"
    #out_1baseBed="${tmp_work_path}/${outname}${inter_name}_bed"
    #out_1baseBed_bothStrand="${tmp_work_path}/${outname}${inter_name}_bed_plus_minus"
    out_1baseBed_minus="${tmp_work_path}/${outname}${inter_name}_bed_minus"
    out_1baseBed_plus="${tmp_work_path}/${outname}${inter_name}_bed_plus"
    out_strand_result="${out_final_part_name}"
    out_strand_result_minus="${tmp_work_path}/${outname}${inter_name}_annotationStrandminus_recal"
    out_strand_result_plus="${tmp_work_path}/${outname}${inter_name}_annotationStrandplus_recal"
    both_strands_site="${tmp_work_path}/${outname}${inter_name}_both_strand_sites"

    #test -s $in_bed || ln -s ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all_6.bed $in_bed  ###### Fri Nov 22 10:59:48 CST 2019 fzc add 
    test -s $out_bed_minus||awk '$NF ~"-"'  $in_bed >$out_bed_minus
    test -s $out_bed_plus|| awk '$NF ~"+"'  $in_bed >$out_bed_plus
    test $method == "bmc" && ( test -e $out_aviput || awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),"het","50.0",$3}' $input_result >${out_aviput} )
    test $method == "gatk" && ( test -e $out_aviput ||${dir_scripts}/convert2annovar.pl -format vcf4 $input_result > ${out_aviput})
    #annotate_variation.pl ${out_aviput} `dirname $in_bed` --buildver ${ref_genome}  -bedfile `basename $in_bed` -dbtype bed -regionanno -out ${out_1baseBed_bothStrand}
    ${dir_scripts}/annotate_variation.pl ${out_aviput} `dirname $out_bed_minus` --buildver ${ref_genome}  -bedfile `basename $out_bed_minus` -dbtype bed -regionanno -out ${out_1baseBed_minus}
    ${dir_scripts}/annotate_variation.pl ${out_aviput} `dirname $out_bed_plus` --buildver ${ref_genome}  -bedfile `basename $out_bed_plus` -dbtype bed -regionanno -out ${out_1baseBed_plus}
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
    echo $method $outname
    #awk_base_change_percentage $method ${out_strand_result}
    rm  ${out_1baseBed_minus}.${ref_genome}_bed ${out_1baseBed_plus}.${ref_genome}_bed ${out_strand_result_plus} ${out_strand_result_minus} ${both_strands_site} $out_aviput
}

pileup_extract_readseq(){
    #sub function of prepare_for_blat_2013_Natmethods_version
    #ref_genome_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/dep_files/bwa_mem_index_${ref_genome}/${ref_genome}_all.fa
    pos=$2
    chr1=`echo $pos|cut -d":" -f1`
    pos0=`echo $pos|cut -d":" -f2|cut -d"-" -f1`
    bam_file=$1
    #tmp_path=~/tmp/blat_filter
    tmp_path=$4
    input=$3
    readid_file=${tmp_path}/${input}/tmp_mpileup_readid_`echo $pos|cut -d"-" -f1`_$(${dir_tools}/random)
    sample_sam_file=${tmp_path}/${input}/view_cluster_sample_`echo $pos|cut -d"-" -f1`_$(${dir_tools}/random).sam

    samtools view -F 1024 $bam_file $pos > ${sample_sam_file} &
    #samtools mpileup -Q 20 -O -r $pos -f $ref_genome_path --output-QNAME $bam_file |awk -F"\t" '{
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
psl_filter(){
    #sub function of prepare_for_blat_2013_Natmethods_version
    ###per site in bmc/gatk
    ###### Sun Oct 6 10:29:14 CST 2019
    ##update:
    #1. if ($4 in fa_title_array && count[$4]==2){ ->if ($4 in fa_title_array){\if (count[$4]==2){
    #2. count1,if (count[i]==1 && "'$chr1'"==$1 && $2 <= '$pos_0_0' &&'$pos_0_0' < $3){->count1_line_record[i]
    psl_with_score_file=$1
    pos_ori=$2 #chr1:367182
    chr1=`echo $pos_ori|cut -d":" -f1`
    pos_0_0=$(expr `echo $pos_ori|cut -d":" -f2` - 1)
    tmp_a=`awk -F "\t" '{
        count[$4]++
        #if ($4 in fa_title_array && count[$4]==2){###### Sun Oct 6 10:00:30 CST 2019 fzc
        if (count[$4]==1){
            count1_line_record[$4]=$0
            }
        if ($4 in fa_title_array){
            if (count[$4]==2){
            ratio=$5/fa_title_array[$4]
            if (ratio < 0.95){
                blat_survive[$4]++
            }
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
                if (count[i]==1 ){
                    split(count1_line_record[i],line_sep_part,"\t")
                    if ("'$chr1'"==line_sep_part[1] && line_sep_part[2] <= '$pos_0_0' &&'$pos_0_0' < line_sep_part[3]){
                        blat_survive[i]++
                    }
                }
            }
            #print length(blat_survive),length(count)
            if (2*length(blat_survive)>length(count)){
                print 1
            }else{print 0}
            }' $psl_with_score_file`
    echo $tmp_a
    
}

"$@"

