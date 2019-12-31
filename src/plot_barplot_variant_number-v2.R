library("ggplot2")

############################### edit start
args<-commandArgs(T)
RNA_off_target = read.table(paste(args[1],"/result_summarized_number_of_twelve_types_RNA_editing_events_of_all_samples.txt",sep=''),header = TRUE)
output_file=args[2]
outname_of_replicates=unlist(strsplit(args[3],','))
#RNA_off_target = read.table('/Users/user/Documents/keti/teacher_chenjia/20190906/res-RNA/5_pipeline_compare/newest_pipeline/01-stats_result--12_kinds_variant--merged-new.txt',header = TRUE)
############################### edit end

## extract Alu, repetitive non-Alu, non-repetitive data for plot
data_Alu=RNA_off_target[,18:29]
data_Alu=sapply(cbind(RNA_off_target[,2],data_Alu), as.numeric)
data_repetitiveNonAlu=RNA_off_target[,30:41]
data_repetitiveNonAlu=sapply(cbind(RNA_off_target[,2],data_repetitiveNonAlu), as.numeric)
data_nonrepetitive=RNA_off_target[,42:53]
data_nonrepetitive=sapply(cbind(RNA_off_target[,2],data_nonrepetitive), as.numeric)

variant_type=c("A>C","A>G","A>U","C>A","C>G","C>U","G>A","G>C","G>U","U>A","U>C","U>G")
dat=t(rbind(data_Alu,data_repetitiveNonAlu,data_nonrepetitive))
rownames(dat) <- c("total_variants",variant_type)

sample_ID_str=rep(RNA_off_target[,1],3)
matrix_for_barplot=rbind(c(rep(sapply(RNA_off_target[,1],as.numeric),3)),c(rep('Alu',nrow(RNA_off_target)),rep('repetitive non-Alu',nrow(RNA_off_target)),rep('non-repetitive',nrow(RNA_off_target))),dat)
rownames(matrix_for_barplot)[1:2]=c("sampleID",'genomic_region')
max_number=max(dat[-1,])

## data preparation for plot
########################### editing start: set sample ID for one treatment
#treatment_sample=c(1,8,9)
treatment_sample=outname_of_replicates
############################ edit end
variant_type_num=12
df_matrix_for_ggplot=data.frame(var_count=c(),var_type=c(),var_region=c())
for(sampleID in treatment_sample){
  if(sampleID!=0){
    #oneSample=matrix_for_barplot[,matrix_for_barplot[1,]==sampleID]
    oneSample=matrix_for_barplot[,sample_ID_str==sampleID]
    genomicRegion=oneSample[2,]
    oneSample=oneSample[-1:-3,]
    var_type=row.names(oneSample)

    oneSample=matrix(sapply(oneSample, as.numeric),variant_type_num,3)
    #print(oneSample)
    for(oneGenomicRegionType in 1:dim(oneSample)[2]){
      oneGenomicRegion=data.frame(var_count=oneSample[,oneGenomicRegionType],var_type=var_type,var_region=rep(genomicRegion[oneGenomicRegionType],variant_type_num) )
      df_matrix_for_ggplot=rbind(df_matrix_for_ggplot,oneGenomicRegion)
      
    }
  }
}

## ggplot: barplot, error bar, point plot
bar_width=0.75
p_variants=ggplot(df_matrix_for_ggplot,aes(x=var_type, y=var_count,fill=var_region))+
  stat_summary(fun.y=mean, geom="bar" ,width=bar_width,position=position_dodge())+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=.5,size=0.4,position=position_dodge(bar_width)) +
  geom_point(size=0.5,shape=21,position=position_jitterdodge(dodge.width = bar_width,jitter.width = 0.1, jitter.height = 0.1))+

  theme(axis.ticks.length=unit(.25, "cm"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ #,size = 0.75))+
  #theme(axis.text.y = element_blank(),axis.text.x = element_blank())+
  xlab("") + ylab("")+
  
  #scale_x_continuous(breaks=c(x_tick),limits = c(0.5,max(x_tick)),expand=c(0,0))+
  #scale_y_continuous(expand = c(0, 0),breaks = seq(0,max_number,Y_axis_breaks))+ 
  coord_cartesian(expand = c(0, 0),ylim = c(0,max_number))+
  scale_fill_manual(values=c("burlywood","burlywood4","#FF9797"))

pdf(output_file,width=8,height=2)
print(p_variants)
dev.off()

