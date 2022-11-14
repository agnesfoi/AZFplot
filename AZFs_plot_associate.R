sample_read<-function(Input,type,filename){
  
  # read
  if (type=='web'){
    df2 <- fromJSON(Input)
    df2 <- as.data.frame(df2)}
  if (type=='raport'){
    df2 <- read.csv(filename)
    }
  
  # extract
  tem<-df2[c('mutation_type', 'chr_grch37_hg19','sample_id')]
  tem<-tem[which(tem$mutation_type=='ChrY'),]
  colnames(tem)[2]<-'infor'
  colnames(tem)[1]<-'chrom'
  
  #tem<-separate_rows(tem,infor,sep="\n")
  tem['chrom']<-gsub('C', 'c' ,tem$chrom)
  tem['infor']<-gsub(' ', '' ,tem$infor)
  tem['loc']<-gsub('.*:', '' ,tem$infor)
  tem['start']<-gsub('-.*', '' ,tem$loc)
  tem['end']<-gsub('.*-', '' ,tem$loc)
  tem['indel']<-gsub('[0-9]', '' ,tem$end)
  tem['start']<-gsub('[^0-9]', '' ,tem$start)
  tem['end']<-gsub('[^0-9]', '' ,tem$end)
  tem['indel']<-gsub('[^A-z]', '' ,tem$indel)
  tem<-tem[which(tem$indel=='Deletion' |tem$indel=='Duplication'),]
  
  sample<-tem[c('chrom', 'start', 'end', 'infor', 'sample_id', 'indel' )]
  
  # transform data
  
  if (type=='web'){
  sample<-sample[order(as.numeric(gsub('A','',sample$sample_id))),]
  sample<-sample[!duplicated(sample$sample_id),] } # 2nd version need to change
  sample<-sample[which(sample$indel=='Deletion'),]
  sample<-sample[is.na(sample$indel)==FALSE,]
  
  sample$start<-as.numeric(sample$start)
  sample$end<-as.numeric(sample$end) 
  sample_a<-sample[which(sample$end<=19567683),]
  sample<-sample[which(sample$end>=19567683),]
  
  # export
  return(list(sample, sample_a))
}







sampDel <- function(i,refer,sample){
  refer$col<-1:nrow(refer)
  refsamp<-refer[,c('col','Amplicon','start','end','str')]
  row.names(refsamp)<-1:nrow(refsamp)
  
  
  d1<- refer[which(sample$start[i]>=refer$start & sample$start[i]<=refer$end),]
  d2<- refer[which(sample$end[i]  >=refer$start & sample$end[i]  <=refer$end),]
  
  
  if (nrow(d1)==0 & nrow(d2)==0)
  { 
    print(paste('deletion start point of',sample$sample_id[i],'is located outside of reference region', sep=' '))
    return('Outlier')
  }
  
  
  else {
    

    if (nrow(d1)!=0 & nrow(d2)!=0)
    {
      
      del1<-factor(d1[,'Amplicon'])
      del2<-factor(d2[,'Amplicon'])
      
      if (length(del1)==1) { a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      if (length(del1)==2) { del1<-del1[1]; a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      if (length(del2)==1) { b<-refer[which(refer$Amplicon==toString(del2)),]$col }
      if (length(del2)==2) { del2<-del2[2]; b<-refer[which(refer$Amplicon==toString(del2)),]$col }
      
      
      
      if (del1!=del2){
        
        if (sample$start[i]==refer$start[a]){ Del1<-a 
        print(paste('deletion start point of',sample$sample_id[i],'is located at the beganning of amplicon',del1, sep=' '))}
        if (sample$start[i]==refer$end[a])  { Del1<-a+1
        print(paste('deletion start point of',sample$sample_id[i],'is located at the end of amplicon',del1, sep=' '))}
        if (sample$start[i] > refer$start[a] & sample$start[i] < refer$end[a]){ Del1<-a+1; refer$end[a]<-sample$start[i]
        print(paste('deletion start point of',sample$sample_id[i],'is located within amplicon',del1, sep=' '))}
        
        if (sample$end[i] > refer$start[b] & sample$start[i] < refer$end[b] & b-a>1){  
          
          if (sample$end[i] == refer$start[b]){ Del2<-b-1
          print(paste('deletion end point of',sample$sample_id[i],'is located at the beganning of amplicon',del2, sep=' '))}
          if (sample$end[i] == refer$end[b])  { Del2<-b
          print(paste('deletion end point of',sample$sample_id[i],'is located at the end of amplicon',del2, sep=' '))}
          if (sample$end[i] > refer$start[b] & sample$start[i] < refer$end[b] & b-a>1){ Del2<-b-1; refer$start[b]<-sample$end[i]
          print(paste('deletion end point of',sample$sample_id[i],'is located within amplicon',del2, sep=' '))}
          
          refer<-refer[-(Del1:Del2),]
          
        } else{
          
          refer$start[b]<-sample$end[i]
          refer<-refer
          print(paste('deletion end point of',sample$sample_id[i],'is located within amplicon',del2, sep=' '))}
        
      }
      
      
      if (del1==del2){
        
        refer[nrow(refer)+1,]<-refer[a,]
        refer$end[a]<-sample$end[i]
        refer$start[nrow(refer)]<-sample$start[i]
        refer$str[nrow(refer)]<-'*'
        
      }
      
      
    }
    
    
    
    
    
    if (nrow(d1)==0 & nrow(d2)!=0)
    { 
      if(sample$start[i]>=19567683)
      {
        
        d1<- refer[which(sample$start[i]<=refer$end),][1,]
        d1_1<- tail(refer[which(sample$start[i]>=refer$end),],1)
        
        del1<-factor(d1[,'Amplicon'])
        del1_1<-factor(d1_1[,'Amplicon'])
        del2<-factor(d2[,'Amplicon'])
        
        if (length(del1)==1) { a<-refer[which(refer$Amplicon==toString(del1)),]$col }  
        if (length(del2)==1) { b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        if (length(del2)==2) { del2<-del2[2]; b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        
        if (del1!=del2){
          
          if (sample$start[i]==refer$start[a]){ Del1<-a 
          print(paste('deletion start point of',sample$sample_id[i],'is located at the beganning of amplicon',del1, sep=' '))}
          if (sample$start[i]==refer$end[a])  { Del1<-a+1
          print(paste('deletion start point of',sample$sample_id[i],'is located at the end of amplicon',del1, sep=' '))}
          if (sample$start[i] > refer$start[a] & sample$start[i] < refer$end[a]){ Del1<-a+1; refer$end[a]<-sample$start[i]
          print(paste('deletion start point of',sample$sample_id[i],'is located within amplicon',del1, sep=' '))}
          if (sample$start[i] > refer$start[a-1] & sample$start[i] < refer$end[a]){ Del1<-a
          print(paste('deletion start point of',sample$sample_id[i],'is located within amplicons',del1,'and',del1_1, sep=' '))}
          
          if (sample$end[i] == refer$start[b]){ Del2<-b-1
          print(paste('deletion end point of',sample$sample_id[i],'is located at the beganning of amplicon',del2, sep=' '))}
          if (sample$end[i] == refer$end[b])  { Del2<-b
          print(paste('deletion end point of',sample$sample_id[i],'is located at the end of amplicon',del2, sep=' '))}
          if (sample$end[i] > refer$start[b] & sample$start[i] < refer$end[b]){ Del2<-b-1; refer$start[b]<-sample$end[i]
          print(paste('deletion end point of',sample$sample_id[i],'is located within amplicon',del2, sep=' '))}
          
          refer<-refer[-(Del1:Del2),]
          
        }
        
      }
      
      else{  
        
        del2<-factor(d2[,'Amplicon'])
        
        if (length(del2)==1) { b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        if (length(del2)==2) { del2<-del2[2]; b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        
        Del2<-b-1; refer$start[b]<-sample$end[i]
        refer<-refer[-(1:Del2),]
        
      }
      
    }
    
    
    
    
    
    if (nrow(d1)!=0 & nrow(d2)==0)
    {
      
      del1<-factor(d1[,'Amplicon'])
      
      if (length(del1)==1) { a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      if (length(del1)==2) { del1<-del1[1]; a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      
      Del1<-a+1; refer$end[a]<-sample$start[i]
      refer<-refer[-(Del1:nrow(refer)),]
    }
    
    
    
    
    
    levels(refer$Amplicon)<-c(levels(refer[nrow(refer)+1,]$Amplicon), 'del')
    refer[nrow(refer)+1,]$Amplicon<-'del'
    refer[nrow(refer),]$start<-sample[i,'start']
    refer[nrow(refer),]$end<-sample[i,'end']
    refer[nrow(refer),]$str<-'*'
    #levels(refer$coleur)<-c(levels(refer$coleur),paste(refer[nrow(refer),]$Amplicon,toString(sprintf("'%s'",'gray98')),sep = '='))
    #refer[nrow(refer),]$coleur<-paste(refer[nrow(refer),]$Amplicon,toString(sprintf("'%s'",'gray98')),sep ='=')
    refer$Amplicon<-factor(refer$Amplicon)
    
    samp<-refer
    return(samp)
    
  }
  
}









sampDel_a <- function (refer,sample){
  
  if (sample$start>=refer$start[1] & sample$end<=refer$end[1]) {
    
    refer$end[1]<-sample$start; refer$start[2]<-sample$end
    refer[3,]<-sample 
    
    }
  
  if (sample$start<=refer$start[1] & sample$end>=refer$start[1] & sample$end<=refer$end[1] ){
    refer<-refer[!duplicated(refer$start),]
    refer$start<-sample$end 
    refer[2,]<-sample
    
  }
  
  
  if (sample$end>=refer$end[1] & sample$start>=refer$start[1] & sample$start<=refer$end[1] ){
    refer<-refer[!duplicated(refer$start),]
    refer$end<-sample$start
    refer[2,]<-sample
    
  }
  
  if (sample$start<refer$start[1] & sample$end>=refer$end[1] ){
    refer<-sample
  }
  

  

  
  return(refer)
  }




  
  
  
  

sampDelweb <- function(i,refer,sample){
  refer$col<-1:nrow(refer)
  refsamp<-refer[,c('col','Amplicon','start','end','str')]
  row.names(refsamp)<-1:nrow(refsamp)
  
  
  d1<- refer[which(sample$start[i]>=refer$start & sample$start[i]<=refer$end),]
  d2<- refer[which(sample$end[i]  >=refer$start & sample$end[i]  <=refer$end),]
  
  
  if (nrow(d1)==0 & nrow(d2)==0)
  { 
   return('Outlier')
  }
  
  
  else {
    
    
    if (nrow(d1)!=0 & nrow(d2)!=0)
    {
      
      del1<-factor(d1[,'Amplicon'])
      del2<-factor(d2[,'Amplicon'])
      
      if (length(del1)==1) { a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      if (length(del1)==2) { del1<-del1[1]; a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      if (length(del2)==1) { b<-refer[which(refer$Amplicon==toString(del2)),]$col }
      if (length(del2)==2) { del2<-del2[2]; b<-refer[which(refer$Amplicon==toString(del2)),]$col }
      
      
      
      if (del1!=del2){
        
        if (sample$start[i]==refer$start[a]){ Del1<-a }
        if (sample$start[i]==refer$end[a])  { Del1<-a+1 }
        if (sample$start[i] > refer$start[a] & sample$start[i] < refer$end[a]){ Del1<-a+1; refer$end[a]<-sample$start[i] }
        
        if (sample$end[i] > refer$start[b] & sample$start[i] < refer$end[b] & b-a>1){  
          
          if (sample$end[i] == refer$start[b]){ Del2<-b-1 }
          if (sample$end[i] == refer$end[b])  { Del2<-b }
          if (sample$end[i] > refer$start[b] & sample$start[i] < refer$end[b] & b-a>1){ Del2<-b-1; refer$start[b]<-sample$end[i] }
          
          refer<-refer[-(Del1:Del2),]
          
        } else{
          
          refer$start[b]<-sample$end[i]
          refer<-refer}
        
      }
      
      
      if (del1==del2){
        
        refer[nrow(refer)+1,]<-refer[a,]
        refer$end[a]<-sample$end[i]
        refer$start[nrow(refer)]<-sample$start[i]
        refer$str[nrow(refer)]<-'*'
        
      }
      
      
    }
    
    
    
    
    
    if (nrow(d1)==0 & nrow(d2)!=0)
    { 
      if(sample$start[i]>=19567683)
      {
        
        d1<- refer[which(sample$start[i]<=refer$end),][1,]
        d1_1<- tail(refer[which(sample$start[i]>=refer$end),],1)
        
        del1<-factor(d1[,'Amplicon'])
        del1_1<-factor(d1_1[,'Amplicon'])
        del2<-factor(d2[,'Amplicon'])
        
        if (length(del1)==1) { a<-refer[which(refer$Amplicon==toString(del1)),]$col }  
        if (length(del2)==1) { b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        if (length(del2)==2) { del2<-del2[2]; b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        
        if (del1!=del2){
          
          if (sample$start[i]==refer$start[a]){ Del1<-a }
          if (sample$start[i]==refer$end[a])  { Del1<-a+1 }
          if (sample$start[i] > refer$start[a] & sample$start[i] < refer$end[a]){ Del1<-a+1; refer$end[a]<-sample$start[i] }
          if (sample$start[i] > refer$start[a-1] & sample$start[i] < refer$end[a]){ Del1<-a }
          
          if (sample$end[i] == refer$start[b]){ Del2<-b-1 }
          if (sample$end[i] == refer$end[b])  { Del2<-b }
          if (sample$end[i] > refer$start[b] & sample$start[i] < refer$end[b]){ Del2<-b-1; refer$start[b]<-sample$end[i] }
          
          refer<-refer[-(Del1:Del2),]
          
        }
        
      }
      
      else{  
        
        del2<-factor(d2[,'Amplicon'])
        
        if (length(del2)==1) { b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        if (length(del2)==2) { del2<-del2[2]; b<-refer[which(refer$Amplicon==toString(del2)),]$col }
        
        Del2<-b-1; refer$start[b]<-sample$end[i]
        refer<-refer[-(1:Del2),]
        
      }
      
    }
    
    
    
    
    
    if (nrow(d1)!=0 & nrow(d2)==0)
    {
      
      del1<-factor(d1[,'Amplicon'])
      
      if (length(del1)==1) { a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      if (length(del1)==2) { del1<-del1[1]; a<-refer[which(refer$Amplicon==toString(del1)),]$col }
      
      Del1<-a+1; refer$end[a]<-sample$start[i]
      refer<-refer[-(Del1:nrow(refer)),]
    }
    
    
    levels(refer$Amplicon)<-c(levels(refer[nrow(refer)+1,]$Amplicon), 'del')
    refer[nrow(refer)+1,]$Amplicon<-'del'
    refer[nrow(refer),]$start<-sample[i,'start']
    refer[nrow(refer),]$end<-sample[i,'end']
    refer[nrow(refer),]$str<-'*'
    refer$Amplicon<-factor(refer$Amplicon)
    samp<-refer
    return(samp)
    
  }
  
}










plot_AZFbc<-function(list,name,type) {
  
  #list<-unique(Filter(Negate(is.null),list))
  
  if (type=='pdf'){pdf(name, width = 40, height = length(list)*1.1, paper='special')}
  if (type=='png'){png(name, width = 2800, height = round(175*3*length(list)/7),units = "px")}
  
  
  plotTracks(list,extend.left=350000,size=0.5, extend.right =320000,
             # ideoTrack
             showBandId=TRUE,fontsize=20,fontcolor="gray40",cex.bands=1,showId=TRUE,
             # axisTrack
             labelPos="below",add53=TRUE,
             # others
             featureAnnotation="feature", groupAnnotation="gr", min.distance=1,
             cex.group=3, rotation.item = 0, cex.feature=1.2, rotation.item=0,
             showTitle=FALSE, title.width=0, collapse=FALSE, just.group="left", 
             fontsize.group=9,
             # markTrack
             'sY145/1192' = 'gray60', 'sY254/255'='gray60',
             # ampTrack
             yel3 ='lightgoldenrod2', yel4 = 'lightgoldenrod2', b5 = 'deepskyblue3', u0 = 'navajowhite2', b6 = 'deepskyblue3', u1 = 'navajowhite2', b1 = 'deepskyblue3', t1 = 'lightseagreen', u2 = 'navajowhite2', t2 = 'lightseagreen', b2 = 'deepskyblue3', u3 = 'navajowhite2', g1 = 'darkseagreen3', r1 = 'firebrick', r2 = 'firebrick', Gr1 = 'gray70', b3 = 'deepskyblue3', yel1 = 'lightgoldenrod2', g2 = 'darkseagreen3', r3 = 'firebrick', r4 = 'firebrick', g3 = 'darkseagreen3', yel2 = 'lightgoldenrod2', b4 = 'deepskyblue3', Gr2 = 'gray70',
             # geneTrack
             CDY2A='lightgoldenrodyellow', HSFY1='slategray2', BPY2='lightyellow', RBMY1E='slategray2', PRY2='lightcyan2', RBMY1J='lightyellow', RBMY1J='lightcyan2', BPY2='slategray2', BPY2='honeydew2', DAZ2='pink', DAZ2='gray99', DAZ2='lightgoldenrodyellow', BPY2B='honeydew2', BPY2B='pink', DAZ4='pink', DAZ4='honeydew2', BPY2B='lightgoldenrodyellow',
             # palTrack
             P1='firebrick2', P2='firebrick2', P3='firebrick2', P4='firebrick2',P5='firebrick2', 
             # regionTrack
             AZFa='firebrick2', AZFb='firebrick2', AZFbc='brown', AZFc='firebrick2',
             # sampTrack
             del='seashell3',
             AZFb_deletion='seashell3', AZFbc_deletion='seashell3', AZFc_deletion='seashell3', b1_b3='seashell3', b2_b3='seashell3', g1_g2='seashell3',r1_r3='seashell3',r2_r4='seashell3',
             #sY127='', sY134='', sY105='', sY121='', sY1224='', sY143='', sY1192='', sY153='', sY145='', 
             Deletion='gray60', del='gray60',缺失='gray60',
              )
  
  dev.off() 
  
}










plot_AZFa<-function(list,name,type) {
  
  list<-unique(Filter(Negate(is.null),list))
  
  if (type=='pdf'){pdf(name, width = 42, height = length(list)*1.2, paper='special')}
  if (type=='png'){png(name, width = 2750, height = 175*3,units = "px")}
  
  plotTracks(list,extend.left=0.05,extend.right=0.05,size=0.5,
             # ideoTrack
             showBandId=TRUE,fontsize=21,fontcolor="gray40",cex.bands=1,showId=TRUE,
             # axisTrack
             labelPos="below",add53=TRUE,
             # others
             featureAnnotation="feature",groupAnnotation="gr",min.distance=1, 
             cex.group=2.2, rotation.item = 0, cex.feature=1.2, rotation.item=0, 
             showTitle=FALSE, title.width=0, collapse=FALSE, just.group="left", 
             showFeatureId=TRUE,
             # regionTrack
             AZFa='firebrick2', AZFb='firebrick2', AZFbc='brown', AZFc='firebrick2',
             # geneTrack
             just.group="right",
             GYG2P1 = 'darkseagreen3',TTTY15 = 'navajowhite2',USP9Y = 'lightseagreen',DDX3Y = 'deepskyblue3',UTY = 'lightgoldenrod2',TMSB4Y = 'firebrick',VCY = 'lightblue3',VCY = 'lightblue3',NLGN4Y = 'paleturquoise4',
             # sampTrack
             AZFa='firebrick2', Deletion='gray60', del='gray60',缺失='gray60',
  )
  
  dev.off()
}







### for shiny



plot_AZFbc_shiny<-function(list) {
  
  #list<-unique(Filter(Negate(is.null),list))
  
  #png(name, width = round(700*4), height = round(175*3*length(list)/10),units = "px")
  
  plotTracks(list,extesamnd.right=0.1,size=0.5,
             # ideoTrack
             showBandId=TRUE,fontsize=20,fontcolor="gray40",cex.bands=1,
             # axisTrack
             labelPos="below",add53=TRUE,
             # others
             featureAnnotation="feature",groupAnnotation="gr",min.distance=1,
             cex.group=1.8, rotation.item = 0, cex.feature=0.9, rotation.item=0,
             showTitle=FALSE, collapse=FALSE,just.group="left",
             # ampTrack
             yel3 ='lightgoldenrod2',yel4 = 'lightgoldenrod2',b5 = 'deepskyblue3',u0 = 'navajowhite2',b6 = 'deepskyblue3',u1 = 'navajowhite2',b1 = 'deepskyblue3',t1 = 'lightseagreen',u2 = 'navajowhite2',t2 = 'lightseagreen',b2 = 'deepskyblue3',u3 = 'navajowhite2',g1 = 'darkseagreen3',r1 = 'firebrick',r2 = 'firebrick',Gr1 = 'gray70',b3 = 'deepskyblue3',yel1 = 'lightgoldenrod2',g2 = 'darkseagreen3',r3 = 'firebrick',r4 = 'firebrick',g3 = 'darkseagreen3',yel2 = 'lightgoldenrod2',b4 = 'deepskyblue3',Gr2 = 'gray70',
             # geneTrack
             CDY2A='lightgoldenrodyellow',HSFY1='slategray2',BPY2='lightyellow',RBMY1E='slategray2',PRY2='lightcyan2',RBMY1J='lightyellow',RBMY1J='lightcyan2',BPY2='slategray2',BPY2='honeydew2',DAZ2='pink',DAZ2='gray99',DAZ2='lightgoldenrodyellow',BPY2B='honeydew2',BPY2B='pink',DAZ4='pink',DAZ4='honeydew2',BPY2B='lightgoldenrodyellow',
             # palTrack
             P1='firebrick2',P2='firebrick2',P3='firebrick2',P4='firebrick2',P5='firebrick2', 
             # sampTrack
             del='gray98')
  
  #dev.off() 
  
}







plot_AZFa_shiny<-function(list) {
  
 # list<-unique(Filter(Negate(is.null),list))
  
  
 # png(name, width = 650*3, height = 175*3,units = "px")
  
  
  plotTracks(list,extend.right=0.1,size=0.5,
             # ideoTrack
             showBandId=TRUE,fontsize=20,fontcolor="gray40",cex.bands=1,
             # axisTrack
             labelPos="below",add53=TRUE,
             # others
             featureAnnotation="feature",groupAnnotation="gr",min.distance=1, 
             cex.group=1.5, rotation.item = 0, cex.feature=1, rotation.item=0,showFeatureId=FALSE,
             showTitle=FALSE, collapse=FALSE,
             # geneTrack
             just.group="right",
             GYG2P1 = 'darkseagreen3',TTTY15 = 'navajowhite2',USP9Y = 'lightseagreen',DDX3Y = 'deepskyblue3',UTY = 'lightgoldenrod2',TMSB4Y = 'firebrick',VCY = 'lightblue3',VCY = 'lightblue3',NLGN4Y = 'paleturquoise4',
             # sampTrack
             Deletion='gray95'
  )
  
  #dev.off()
  
}

