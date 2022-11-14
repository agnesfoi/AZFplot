###### Global variables

library(Gviz)
library(shiny)
library(tidyr)

load("AZFbc_plot_Reference.RData")
  
  
######  define the UI  ######

ui<-fluidPage(
  
  sidebarPanel(fileInput("file", "Choose CSV File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               width = 2),
  
    mainPanel(
    h1("AZFbc Deletion",align="left"),
    plotOutput('plotAZF_1'),
    plotOutput('plotAZF_2') 
  )
) 




####  define the Server code  ######

server <- function(input,output,session){
  
  output$plotAZF_1 <- renderImage({
    
    req(input$file)   

    tryCatch(
      {
        df2 <- read.csv(input$file$datapath)
        
        tem<-df2[c('mutation_type', 'chr_grch37_hg19','sample_id')]
        tem<-tem[which(tem$mutation_type=='ChrY'),]
        colnames(tem)[2]<-'infor'
        colnames(tem)[1]<-'chrom'
        
        tem<-separate_rows(tem,infor,sep="\n")
        tem['chrom']<-gsub('C', 'c' ,tem$chrom)
        tem['infor']<-gsub(' ', '' ,tem$infor)
        tem['loc']<-gsub('.*:', '' ,tem$infor)
        tem['start']<-gsub('-.*', '' ,tem$loc)
        tem['end']<-gsub('.*-', '' ,tem$loc)
        tem['indel']<-gsub('[0-9]', '' ,tem$end)
        tem['end']<-gsub('[^0-9]', '' ,tem$end)
        tem['indel']<-gsub('缺失', 'Deletion' ,tem$indel)
        tem['indel']<-gsub('重复', 'Duplication' ,tem$indel)
        tem['indel']<-gsub('[^A-z]', '' ,tem$indel)
        tem<-tem[which(tem$indel=='Deletion' |tem$indel=='Duplication'),]
        
        sample<-tem[c('chrom', 'start', 'end', 'infor', 'sample_id', 'indel' )]
        rm(df2,tem)
        sample<-sample[order(sample$sample_id),]
        sample<-sample[which(sample$indel=='Deletion'),]
        sample<-sample[is.na(sample$indel)==FALSE,]
        sample_a<-sample[which(sample$end<19567683),]
        sample<-sample[which(sample$end>=19567683),]
        
        
        if (nrow(sample)!=0){
          
          n= nrow(sample)
          
          for (i in 1:n){
            
            
            if (is.na(sample$indel[i])==FALSE && sample$indel[i]=="Deletion"){
              
              sampleI<-sample[i,]
              n<-nrow(sample[which(sample$sample_id==sampleI$sample_id),])
              
              
              if (n>1){
                
                sampleI<-sample[which(sample$sample_id==sampleI$sample_id),]
                sample <-sample[-((i+1):(n+1)),]
                
                for (a in 1:n){
                  if (a==1){samp<-sampDel(i=a, refer=ref_bc, sample=sampleI)
                  } else {
                    ref_mul<-samp ; samp<-sampDel(i=a, refer=ref_mul, sample=sampleI)}}
              } 
              
              if (n==1) { samp<-sampDel(i, refer=ref_bc, sample=sample) }
              
              
              
              if (samp!='Outlier'){
                
                if (samp[samp$Amplicon=='del',]$start >= 19567683)
                {
                  sampTrack_bc <-AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                                                 genome="hg19", col.title="black", name=paste('Sample',i),
                                                 fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                                 group=sample[i,]$sample_id)
                  feature(sampTrack_bc)<- as.character(samp$Amplicon) 
                  
                  #assign(paste0("sampTrack_bc"),temTrack_bc)
                  trackList_bc1[[4+i]]<-sampTrack_bc
                  rm(sampTrack_bc)
                  
                }
                
                
                if (samp[samp$Amplicon=='del',]$start < 19567683)
                {
                  sampTrack_bc <-AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                                                 genome="hg19", col.title="black", name=paste('Sample',i),
                                                 fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                                 group=sample[i,]$sample_id)
                  feature(sampTrack_bc)<- as.character(samp$Amplicon) 
                  
                  #assign(paste0("sampTrack_bc"),temTrack_bc)
                  trackList_bc2[[4+i]]<-sampTrack_bc
                  rm(sampTrack_bc)
                  
                }
                
                
              }
              
              if (samp=='Outlier') { 
                sample_a<-rbind(sample_a,sample[i,])
                print(paste0('sample ',sample[i,]$sample_id,' : please check for AZFa'))}
              
              
            }}}
        
        trackList_bc1<-unique(Filter(Negate(is.null),trackList_bc1))
        trackList_bc2<-unique(Filter(Negate(is.null),trackList_bc2))
        
        width<-1200
        height<-77*length(trackList_bc1)
        outfile<-tempfile(fileext = ".png")
        png(outfile, width=width, height=height)
        plot_AZFbc_shiny(list=trackList_bc1) 
        dev.off()
       
        list(
          src = outfile,
          contentType = "image/png",
          width = width,
          height = height,
          alt = NULL
        )
        
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
    
  
    
  } )

  output$plotAZF_2 <- renderImage({
    
    req(input$file)   
    
    tryCatch(
      {
        df2 <- read.csv(input$file$datapath)
        
        tem<-df2[c('mutation_type', 'chr_grch37_hg19','sample_id')]
        tem<-tem[which(tem$mutation_type=='ChrY'),]
        colnames(tem)[2]<-'infor'
        colnames(tem)[1]<-'chrom'
        
        tem<-separate_rows(tem,infor,sep="\n")
        tem['chrom']<-gsub('C', 'c' ,tem$chrom)
        tem['infor']<-gsub(' ', '' ,tem$infor)
        tem['loc']<-gsub('.*:', '' ,tem$infor)
        tem['start']<-gsub('-.*', '' ,tem$loc)
        tem['end']<-gsub('.*-', '' ,tem$loc)
        tem['indel']<-gsub('[0-9]', '' ,tem$end)
        tem['end']<-gsub('[^0-9]', '' ,tem$end)
        tem['indel']<-gsub('??????', 'Deletion' ,tem$indel)
        tem['indel']<-gsub('??????', 'Duplication' ,tem$indel)
        tem['indel']<-gsub('[^A-z]', '' ,tem$indel)
        tem<-tem[which(tem$indel=='Deletion' |tem$indel=='Duplication'),]
        
        sample<-tem[c('chrom', 'start', 'end', 'infor', 'sample_id', 'indel' )]
        rm(df2,tem)
        sample<-sample[order(sample$sample_id),]
        sample<-sample[which(sample$indel=='Deletion'),]
        sample<-sample[is.na(sample$indel)==FALSE,]
        sample_a<-sample[which(sample$end<19567683),]
        sample<-sample[which(sample$end>=19567683),]
        
        
        if (nrow(sample)!=0){
          
          n= nrow(sample)
          
          for (i in 1:n){
            
            
            if (is.na(sample$indel[i])==FALSE && sample$indel[i]=="Deletion"){
              
              sampleI<-sample[i,]
              n<-nrow(sample[which(sample$sample_id==sampleI$sample_id),])
              
              
              if (n>1){
                
                sampleI<-sample[which(sample$sample_id==sampleI$sample_id),]
                sample <-sample[-((i+1):(n+1)),]
                
                for (a in 1:n){
                  if (a==1){samp<-sampDel(i=a, refer=ref_bc, sample=sampleI)
                  } else {
                    ref_mul<-samp ; samp<-sampDel(i=a, refer=ref_mul, sample=sampleI)}}
              } 
              
              if (n==1) { samp<-sampDel(i, refer=ref_bc, sample=sample) }
              
              
              
              if (samp!='Outlier'){
                
                if (samp[samp$Amplicon=='del',]$start >= 19567683)
                {
                  sampTrack_bc <-AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                                                 genome="hg19", col.title="black", name=paste('Sample',i),
                                                 fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                                 group=sample[i,]$sample_id)
                  feature(sampTrack_bc)<- as.character(samp$Amplicon) 
                  
                  #assign(paste0("sampTrack_bc"),temTrack_bc)
                  trackList_bc1[[4+i]]<-sampTrack_bc
                  rm(sampTrack_bc)
                  
                }
                
                
                if (samp[samp$Amplicon=='del',]$start < 19567683)
                {
                  sampTrack_bc <-AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                                                 genome="hg19", col.title="black", name=paste('Sample',i),
                                                 fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                                 group=sample[i,]$sample_id)
                  feature(sampTrack_bc)<- as.character(samp$Amplicon) 
                  
                  #assign(paste0("sampTrack_bc"),temTrack_bc)
                  trackList_bc2[[4+i]]<-sampTrack_bc
                  rm(sampTrack_bc)
                  
                }
                
                
              }
              
              if (samp=='Outlier') { 
                sample_a<-rbind(sample_a,sample[i,])
                print(paste0('sample ',sample[i,]$sample_id,' : please check for AZFa'))}
              
              
            }}}
        
        trackList_bc1<-unique(Filter(Negate(is.null),trackList_bc1))
        trackList_bc2<-unique(Filter(Negate(is.null),trackList_bc2))
        
        width<-1200
        height<-77*length(trackList_bc2)
        outfile<-tempfile(fileext = ".png")
        png(outfile, width=width, height=height)
        plot_AZFbc_shiny(list=trackList_bc2) 
        dev.off()
        
        list(
          src = outfile,
          contentType = "image/png",
          width = width,
          height = height,
          alt = NULL
        )
        
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
    
    
    
  })

  
  session$onSessionEnded(stopApp)
  
} 



shinyApp(ui, server)


