# -*- coding: UTF-8 -*-
'''
ClassName: file_name
@author: Yi-Chun Xiong
@version 1.0 
Create Time: 2019/09/30 10:41:18
Description: Summarize the number of all 12 types of RNA editing for each sample into a Tab-delimited text file.

'''

import sys
import re
import os


# name='custom_pipeline_our_previous'
# pwd='/Users/user/Documents/keti/teacher_chenjia/20190906/res-RNA/5_pipeline_compare/'+name+'/'
# fr_12_type_variants=open(pwd+'01-stats_result--12_kinds_variant.txt','r')
# fw_arranged_12_type_variants=open(pwd+'01-stats_result--12_kinds_variant--merged-new.txt','w')

fr_12_type_variants=open(sys.argv[1],'r')
fw_arranged_12_type_variants=open(sys.argv[2],'w')



dict_Alu={}
dict_Repetitive_non_Alu={}
dict_Nonrepetitive={}
dict_all={'dict_Alu':dict_Alu,'dict_Nonrepetitive':dict_Nonrepetitive,'dict_Repetitive_non_Alu':dict_Repetitive_non_Alu}

#variants_type=["A>G","T>C","C>T","G>A","C>A","G>T","C>G","G>C","A>T","T>A","A>C","T>G"]
variants_type=["A>C","A>G","A>U","C>A","C>G","C>U","G>A","G>C","G>U","U>A","U>C","U>G"]
fw_arranged_12_type_variants.write("sampleID"+"\t"+"total_variants"+"\tTotal_"+"\tTotal_".join(variants_type)+"\t"+"Alu_total"+"\t"+"Repetitive_non_Alu_total"+"\t"+"Nonrepetitive_total"+"\tAlu_"+"\tAlu_".join(variants_type)+"\tRepetitive_non_Alu_"+"\tRepetitive_non_Alu_".join(variants_type)+"\tNonrepetitive_"+"\tNonrepetitive_".join(variants_type)+'\n')

current_genomic_region=''
for oneline in fr_12_type_variants.readlines():
    if oneline.strip() != "":
        if re.match(r'\S+: Alu', oneline):
            
            if current_genomic_region!="":
                ## set 0 for the specified variant type in genomic region without variant
                for key in ['dict_Alu','dict_Repetitive_non_Alu','dict_Nonrepetitive']:
                    for one_type_variant in variants_type:
                        if not dict_all[key].has_key(one_type_variant):
                            dict_all[key][one_type_variant]=0
                            
                ## output total specified variant type 
                fw_arranged_12_type_variants.write(str(sum(dict_all['dict_Alu'].values())+sum(dict_all['dict_Repetitive_non_Alu'].values())+sum(dict_all['dict_Nonrepetitive'].values()))+"\t")
                for one_type_variant in variants_type:
                    each_variant_merged=dict_all['dict_Alu'][one_type_variant]+dict_all['dict_Nonrepetitive'][one_type_variant]+dict_all['dict_Repetitive_non_Alu'][one_type_variant]
                    fw_arranged_12_type_variants.write(str(each_variant_merged)+'\t')
                    
                ## output specified variant type in Alu, repetitive non-Alu, non-repetitive region
                fw_arranged_12_type_variants.write(str(sum(dict_all['dict_Alu'].values()))+"\t"+str(sum(dict_all['dict_Repetitive_non_Alu'].values()))+"\t"+str(sum(dict_all['dict_Nonrepetitive'].values()))+"\t")
                for key in ['dict_Alu','dict_Repetitive_non_Alu','dict_Nonrepetitive']:
                    for one_type_variant in variants_type:
                        if dict_all[key].has_key(one_type_variant):
                            fw_arranged_12_type_variants.write(str(dict_all[key][one_type_variant])+'\t')
                        else:
                            fw_arranged_12_type_variants.write('0\t')
                            dict_all[key][one_type_variant]=0
                fw_arranged_12_type_variants.write('\n')
                dict_all['dict_Alu'].clear()
                dict_all['dict_Nonrepetitive'].clear()
                dict_all['dict_Repetitive_non_Alu'].clear()
            
            fw_arranged_12_type_variants.write(oneline.split(":")[0]+'\t')
            current_genomic_region='dict_Alu'
            
        elif re.match(r'\S+: Nonrepetitive', oneline):
            current_genomic_region='dict_Nonrepetitive'
        
        elif re.match(r'\S+: Repetitive_non_Alu', oneline):
            current_genomic_region='dict_Repetitive_non_Alu'
        else:
            onelineArr=re.split(r'\s+',oneline)
            #print onelineArr[2]+"####"+onelineArr[1]+"****"
            if onelineArr[2] in variants_type:
                dict_all[current_genomic_region][onelineArr[2]]=int(onelineArr[1])
    
else:
    ## set 0 for the specified variant type in genomic region without variant
    for key in ['dict_Alu','dict_Repetitive_non_Alu','dict_Nonrepetitive']:
        for one_type_variant in variants_type:
            if not dict_all[key].has_key(one_type_variant):
                dict_all[key][one_type_variant]=0
                
    ## output total specified variant type 
    fw_arranged_12_type_variants.write(str(sum(dict_all['dict_Alu'].values())+sum(dict_all['dict_Repetitive_non_Alu'].values())+sum(dict_all['dict_Nonrepetitive'].values()))+"\t")
    for one_type_variant in variants_type:
        each_variant_merged=dict_all['dict_Alu'][one_type_variant]+dict_all['dict_Nonrepetitive'][one_type_variant]+dict_all['dict_Repetitive_non_Alu'][one_type_variant]
        fw_arranged_12_type_variants.write(str(each_variant_merged)+'\t')
    
    ## output specified variant type in Alu, repetitive non-Alu, non-repetitive region
    fw_arranged_12_type_variants.write(str(sum(dict_all['dict_Alu'].values()))+"\t"+str(sum(dict_all['dict_Repetitive_non_Alu'].values()))+"\t"+str(sum(dict_all['dict_Nonrepetitive'].values()))+"\t")
    for key in ['dict_Alu','dict_Repetitive_non_Alu','dict_Nonrepetitive']:
        for one_type_variant in variants_type:
            if dict_all[key].has_key(one_type_variant):
                fw_arranged_12_type_variants.write(str(dict_all[key][one_type_variant])+'\t')
            else:
                fw_arranged_12_type_variants.write('0\t')
                dict_all[key][one_type_variant]=0
    fw_arranged_12_type_variants.write('\n')
    dict_all['dict_Alu'].clear()
    dict_all['dict_Nonrepetitive'].clear()
    dict_all['dict_Repetitive_non_Alu'].clear()


fr_12_type_variants.close()
fw_arranged_12_type_variants.close()

