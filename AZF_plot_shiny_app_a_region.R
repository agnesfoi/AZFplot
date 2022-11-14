###### Global variables

library(Gviz)
library(shiny)
library(tidyr)

load("AZFa_plot_Reference.RData")


###### define the UI


ui<-fluidPage(
  
  sidebarPanel(fileInput("file", "Choose CSV File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               width = 2),
  
  
  mainPanel(
    h1("AZFa Deletion",align="left"),
    plotOutput('plotAZF')
    
  )
  
  
) 




#### define the Server code

server <- function(input, output, session){

   
  output$plotAZF <- renderImage({
    
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
        
        
        
        trackList_a<-list(ideoTrack_a, axisTrack_a, geneTrack_a)     
        
        if (nrow(sample_a)!=0) {
          
          for (i in 1:nrow(sample_a)){
            
            samp<-sample_a[i,]
            
            sampTrack_a <- AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                                           genome="hg19", col.title="black", name='indel',
                                           fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                           group=sample_a[i,]$sample_id)
            feature(sampTrack_a)<- as.character(samp$indel)
            
            trackList_a[[3+i]]<-sampTrack_a
            rm(sampTrack_a)
          }
          
          
        }
        
        trackList_a<-unique(Filter(Negate(is.null),trackList_a))
        
        width<-1200
        height<-77*length(trackList_a)
        outfile<-tempfile(fileext = ".png")
        png(outfile, width=width, height=height)
        plot_AZFa_shiny(list=trackList_a) 
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
  
} 





shinyApp(ui, server)


