
#chromosome_size=data.frame(CHR=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,30,31,32),BP=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569))
chromosome_size=data.frame(CHR=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,30),BP=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895))


## function for multiple sample Manhattan plot
## df_multiple_sample: data.frame including 
#### CHR:chromosome in number,
#### BP: chromosome coordinate
#### P: editing ratio or P-value in GWAS
#### sample: sample id
## sampleOrder: sample ids in df_multiple_sample$sample 
## sampleColor: the color for each sample
## one_sample_width: one sample width in Manhattan plot
## gap_between_sample: gap between samples in Manhattan plot
## dot_size: 
multiple_sample_Manhattan_plot <- function(df_multiple_sample,sampleOrder,sampleColor,sampleXaxisLabels=c(),one_sample_width=100,gap_between_sample=40,dot_size=1.25,output_jpeg=FALSE){
  
  
  #df_multiple_sample=df_multiple_sample[order(df_multiple_sample$sample,df_multiple_sample$CHR),]
  
  axis.set=data.frame()
  df=data.frame()
  x_axis=0
  for(sampleID in sampleOrder){
    oneSample=df_multiple_sample[df_multiple_sample$sample==sampleID,]
    oneSample=oneSample[order(oneSample$CHR),]
    ### calculate relative position of each dot in one sample according to the genomic coordinate.
    oneSample$BPcum <- NA
    s <- 0
    nbp <- c()
    for (i in unique(oneSample$CHR)){
      nbp[i] <- max(oneSample[oneSample$CHR == i,]$BP)
      oneSample[oneSample$CHR == i,"BPcum"] <- oneSample[oneSample$CHR == i,"BP"] + s
      s <- s + nbp[i]
    }
    
    ### compress one sample's dot into one_sample_width in the x axis
    oneSample$BPcum=one_sample_width*oneSample$BPcum/max(oneSample$BPcum)+x_axis
    ### calculate x axis coordinate of where x labs should be positioned
    axis.set=rbind(axis.set,data.frame(sample=c(sampleID),center=c((x_axis+x_axis+one_sample_width)/2)))
    ### add one sample to the data.frame for plot
    df=rbind(df,oneSample)
    x_axis=x_axis+one_sample_width+gap_between_sample
    
  }
  
  ### plot in the order of sampleOrder
  df=within(df,{
    sample=factor(sample,levels=sampleOrder)
  })
  if(length(sampleXaxisLabels)>0){
    axis.set$sample=sampleXaxisLabels
  }
  
  ### create dotplot of manhattan
  manhplot <- ggplot(df, aes(x=BPcum, y=P)) +
    geom_point(aes(color=as.factor(sample)), alpha = 0.75, size = dot_size) +
    scale_color_manual(values = sampleColor) +
    #scale_x_continuous(expand = c(0,gap_between_sample))+
    scale_x_continuous(expand = c(0,gap_between_sample),label = axis.set$sample, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2)) +
    coord_cartesian(ylim = c(0,1))+
    #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +p
    theme(axis.ticks.x=element_blank())+ #.length.x=unit(0, "cm"))+ #axis.ticks.length.y=unit(.25, "cm"),
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ #,size = 0.75))+
    #theme(axis.text.y = element_blank(),axis.text.x = element_blank())+
    xlab("") +ylab("")+ # ylab("edit ratio")  #+ ggtitle("this is title")
    theme(legend.position = 'none')
    #theme(axis.line= element_blank(),axis.ticks=element_blank())
    #scale_y_continuous(expand = c(0,0), breaks = seq(0,1,1))
    
    if(output_jpeg){
      manhplot=manhplot+theme(axis.line= element_blank(),axis.ticks=element_blank())
    }
  
    return(manhplot)
  
 }


## one sample manhattan plot by ggplot:  http://www.danielroelfs.com/coding/manhattan_plots/
one_sample_Manhattan_plot <- function(df_multiple_sample,sampleID,sampleColor,dot_size=1.25,output_jpeg=FALSE){
  ## extract editing ratio of required sample
  oneSample=df_multiple_sample[df_multiple_sample$sample==sampleID,]
  oneSample=oneSample[order(oneSample$CHR),]
  
  ## calculate the X coordinate for each variant
  nCHR <- length(unique(oneSample$CHR))
  oneSample$BPcum <- NA
  s <- 0
  nbp <- c()
  vline=c()
  
  coordinate_x_axis_chromosome=data.frame(CHR=c(),BPcum=c())
  coordinate_start=1
  for (i in unique(chromosome_size$CHR)){
    nbp[i] <- max(chromosome_size[chromosome_size$CHR == i,]$BP)
    oneSample[oneSample$CHR == i,"BPcum"] <- oneSample[oneSample$CHR == i,"BP"] + s
    s <- s + nbp[i]
    vline=c(vline,s)
    coordinate_x_axis_chromosome=rbind(coordinate_x_axis_chromosome,data.frame(CHR=c(i),BPcum=c(coordinate_start)))
    coordinate_x_axis_chromosome=rbind(coordinate_x_axis_chromosome,data.frame(CHR=c(i),BPcum=c(s)))
    coordinate_start=s
  }
  
  ## calculate coordinate for X axis lable
  # axis.set <- oneSample %>%
  #   group_by(CHR) %>%
  #   summarize(center = (max(BPcum) + min(BPcum)) / 2)
  x_axis_labels <- coordinate_x_axis_chromosome %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  
  x_axis_labels[x_axis_labels$CHR==30,1]='X'
  x_axis_labels[x_axis_labels$CHR==31,1]='Y'
  x_axis_labels[x_axis_labels$CHR==32,1]='M'

  ## create dotplot of manhattan
  manhplot <- ggplot(oneSample, aes(x=BPcum, y=P)) +
    #geom_vline(xintercept = vline)+
    geom_point(aes(color=as.factor(CHR)), alpha = 0.75, size = dot_size) +
    scale_color_manual(values = rep(sampleColor, nCHR)) +
    scale_x_continuous(label = x_axis_labels$CHR, breaks = x_axis_labels$center) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2)) +
    coord_cartesian(ylim = c(0,1))+
    #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
    
    theme(axis.ticks.x=element_blank())+ #length.x=unit(0, "cm"))+  #axis.ticks.length.y=unit(.25, "cm"),
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ #,size = 0.75))+
    #theme(axis.text.y = element_blank(),axis.text.x = element_blank())+
    xlab("") +ylab("") # ylab("edit ratio")#  #+ ggtitle("this is title")
    #theme(legend.position = 'none')
    #theme(axis.line= element_blank(),axis.ticks=element_blank())
    #scale_y_continuous(expand = c(0,0), breaks = seq(0,1,1))
  
  if(output_jpeg){
    manhplot=manhplot+theme(axis.line= element_blank(),axis.ticks=element_blank())
  }else{
    #manhplot=manhplot+geom_vline(xintercept = vline)
  }
    
  return(manhplot) 
  
  
}
