library("ggplot2")
library("dplyr")
#source('/Users/user/Documents/workspace/R/RADAR_plot/Function/fun_Manhattan_plot.R')
source('./src/plot_fun_Manhattan_plot.R')

args<-commandArgs(T)
#args=c('/Users/user/Documents/workspace/R/RADAR_plot/tmp_merged_samples.vcf','/Users/user/Documents/workspace/R/RADAR_plot/tmp-mahattan.pdf',"20190906_1,20190906_2,20190906_4","#919191,#FF3F00,#FF3F00")
#args=c('/Users/user/Documents/workspace/R/RADAR_plot/tmp_merged_samples.vcf','/Users/user/Documents/workspace/R/RADAR_plot/tmp-mahattan.pdf',"20190906_2","#919191,#FF3F00")

#workspace='/Users/user/Documents/keti/teacher_chenjia/20191016/RNA/res/manhattan_plot/'
input_df_manhattan=args[1]
output_file=args[2]
outname_of_samples=args[3]
color_of_samples=args[4]
df_manhattan <- read.table(input_df_manhattan, header=FALSE)

###################### data format
#"chr","BP","SNV","ref_alt","P(editng ratio)","sampleID"
# chr1	6522535	C>T	23,4	0.148148	10
# chr1	52373059	C>T	38,4	0.0952381	10
# chr1	179842190	C>T	36,7	0.162791	10
# chr15	34980672	C>T	30,4	0.117647	10
# chr15	66336541	C>T	31,3	0.0882353	10

colnames(df_manhattan)=c("chr","BP","SNV","ref_alt","P","sample")
df_manhattan$CHR=substring(df_manhattan$chr,4)
df_manhattan$CHR[df_manhattan$CHR=="X"]=30
df_manhattan$CHR[df_manhattan$CHR=="Y"]=31
df_manhattan$CHR[df_manhattan$CHR=="M"]=32
df_manhattan$CHR=as.numeric(df_manhattan$CHR)

#sample=c(16,17,1,2,3,4,6,7,10,11,12,13,14,15)
#sampleColor=c("#919191","#919191","#FF3F00","#FF3F00","#FF3F00","#00BEF2","#00BEF2","#521B93","#E1CABA","#E1CABA","#E1CABA","#FFA9FF","#FFA9FF","#FFA9FF")
sample=unlist(strsplit(outname_of_samples,","))
sampleColor=unlist(strsplit(color_of_samples,","))

if(length(sample)>1){
  
  ############# all treatments
  plot_Manhattan=multiple_sample_Manhattan_plot(df_manhattan,sample,sampleColor)
  #pdf(paste(workspace,"pdf_manhattan_01_all_sample.pdf",sep=""),width=10,height=2)
  pdf(output_file,width=0.4*length(sample)+5,height=2)
  print(plot_Manhattan)
  dev.off()

}else{
  ########### one sample Manhattan plot
  #sampleID=c(6)
  #sampleColor=c("#FFA8FF","#FFD0FF")
  
  plot_Manhattan=one_sample_Manhattan_plot(df_manhattan,sample,sampleColor)
  #pdf(paste(workspace,"pdf_manhattan_04_one_sample_BEACON2_1st.pdf",sep=""),width=10,height=2)
  pdf(output_file,width=10,height=2)
  print(plot_Manhattan)
  dev.off()
}

########### one treatment
#sample=c(1,2,3)
#sampleColor=c("#FF3F00","#FF3F00","#FF3F00")
#plot_Manhattan=multiple_sample_Manhattan_plot(df_manhattan,sample,sampleColor)
#pdf(paste(workspace,"pdf_manhattan_03_",treat_name,".pdf",sep=""),width=6,height=2)
#print(plot_Manhattan)
#dev.off()
