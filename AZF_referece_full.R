########### 
########### 
########### This document is for creating reference gviz Track objects, for AZF ploting
########### 
########### ### reference includes AZFa and AZFbc region ###
########### 
########### AZFa_ref_del.png : classic deletion in AZFa region
########### AZFa_ref_no_deletion.png : example of sample without deletion in AZFa region
########### AZFa_sample.pdf : example of sample with deletion in AZFa region
########### 
########### AZFbc_ref_del.png : classic deletion in AZFbc region
########### AZFbc_ref_no_deletion.png : example of sample without deletion in AZFbc region
########### AZFbc_sample_1/2.pdf : example of sample with deletion in AZFbc region
########### 
########### ### each figure is composed by several gviz Track objects ###
########### IdeogramTrack : chromosome Ideogram plot [ideoTrack]
########### dataTrack : region zoom plot (blue trangle line [dataTrack]) + separate line between reference and sample [sepTrack]
########### AxisTrack : genome location plot ([axisTrack] gray line with location in mb)
########### AnnotationTrack (annotation objects plot) : region ([regTrack] AZFa or bc); makers([markTrack] AZFa or bc); amplicons([ampTrack] AZFbc); genes([geneTrack] AZFa)
########### 
########### ### each xxTrack object is a list, containning several xxTrack objects ###
########### eg. markTrack : list of 3 annotationTracks markTrack[[1]], markTrack[[2]], markTrack[[3]], named markTrack_a (markers in AZFa) + markTrack_bc (markers in AZFbc) + markTrack_abc (markers in AZFabc) respectively.
########### special eg. trackList_a contains all xxTrack objects for plotting sample deletion in AZFa region
########### special eg. trackList_bc1/2 contains all xxTrack objects ([IdeogramTrack, dataTrack, AxisTrac, regTrack, markTrack, nosampTrack, sepTrack]) for plotting sample deletion in AZFbc region
########### 
########### 



### load library

# http://breakpointmapper.wi.mit.edu/mapper.html
# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# https://www.rdocumentation.org/packages/Gviz/versions/1.16.3/topics/settings
# https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html

Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin"

library(Gviz)
library(GenomicRanges)
#library(showtext)
#library(GenomicFeatures)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg19)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db)
#library(CNVRanger)
#library(dplyr)
#library(ucscChromosomeNames)



#########################################################################
#####################      AZFa & AZFbc plot        #####################
#########################################################################

#source('Sample_detection.R')
load("Reference_AZFbc/AZFbc_plot_Reference.RData") # load Rdata reference





######################## create reference plot ########################

### plot abc
plot_AZFbc(list=list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], dataTrack[[1]], geneTrack[[2]], dataTrack[[2]], 
                     markTrack[[1]], ampTrack[[1]], ampTrack_del[[1]],ampTrack_del[[2]],ampTrack_del[[3]],ampTrack_del[[4]],
                     ampTrack_del[[5]],ampTrack_del[[6]],ampTrack_del[[7]],ampTrack_del[[8]]), name='test_abc.png',type='png')


### plot only bc
#trackList_bc1<-list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], ampTrack_del[[1]],ampTrack_del[[2]],ampTrack_del[[3]],ampTrack_del[[4]],
#                   ampTrack_del[[5]],ampTrack_del[[6]],ampTrack_del[[7]],ampTrack_del[[8]])
trackList_bc1<-list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], sepTrack[[2]])
plot_AZFbc(list=trackList_bc1, name='AZFbc_sample_1.pdf',type='pdf')


#trackList_bc2<-list(ideoTrack[[2]], dataTrack[[4]], axisTrack[[2]], regTrack[[2]], markTrack[[2]],  ampTrack[[1]], ampTrack_del[[1]],ampTrack_del[[2]],ampTrack_del[[3]],ampTrack_del[[4]],
#                    ampTrack_del[[5]],ampTrack_del[[6]],ampTrack_del[[7]],ampTrack_del[[8]])
trackList_bc2<-list(ideoTrack[[2]], dataTrack[[4]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], sepTrack[[2]])

plot_AZFbc(list=trackList_bc2, name='AZFbc_sample_2.pdf',type='pdf')


# neg ref bc
trackList_bc<-list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], sepTrack[[2]], ampTrack[[2]])
plot_AZFbc(list=trackList_bc, name='AZFbc_ref_no_deletion.png', type='png')

# reference deletion bc
trackList_bc<-list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], sepTrack[[2]], ampTrack[[1]], ampTrack_del[[1]],
                   ampTrack_del[[2]],ampTrack_del[[3]],ampTrack_del[[4]],ampTrack_del[[5]],ampTrack_del[[6]],ampTrack_del[[7]],
                   ampTrack_del[[8]],ampTrack_del[[9]],ampTrack_del[[10]],ampTrack_del[[11]],ampTrack_del[[12]],ampTrack_del[[13]],
                   ampTrack_del[[14]],ampTrack_del[[16]],ampTrack_del[[17]],ampTrack_del[[18]])

plot_AZFbc(list=trackList_bc,name='AZFbc_ref_del.png', type='png')#,extend.left=500000,cex.group=3,cex.feature=1.2)

.



### plot only a 
trackList_a<-list(ideoTrack[[1]], dataTrack[[1]], axisTrack[[1]], regTrack[[1]], markTrack[[1]], geneTrack[[1]], sepTrack[[1]])
plot_AZFa(list=trackList_a, name='AZFa_sample.pdf',type='pdf')

trackList_a<-list(ideoTrack[[1]], dataTrack[[1]], axisTrack[[1]], regTrack[[1]], markTrack[[1]], geneTrack[[1]], sepTrack[[1]] , regTrack[[4]])
plot_AZFa(list=trackList_a, name='AZFa_ref_no_deletion.png',type='png')


# reference deletion a
trackList_a<-list(ideoTrack[[1]], dataTrack[[1]], axisTrack[[1]], regTrack[[1]], markTrack[[1]], geneTrack[[1]], sepTrack[[1]] , ampTrack_del[[19]], ampTrack_del[[20]])
plot_AZFa(list=trackList_a, name='AZFa_ref_del.png',type='png')












######################## Import data for Track input ########################


# import

sample <- read.delim("AZF20190620.txt")
sample$'sample_ID'<-paste('sample', 1:nrow(sample))
sample$Indel<-gsub(".* ","",sample$infor)

# read region
st<-read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_region_abc1.txt")

# read deletion refer
#ref_bc <- read.delim("reference/AZF_reference_20190620.txt")
#ref_bc<-ref_bc[order(ref_bc$start),]
#ref_a <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFa/hg19_refGene_AZFa.txt")

# read genes
ref_a_gene <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFa/hg19_refGene_AZFa_gene.txt")
ref_a_gene<-ref_a_gene[which(ref_a_gene$end<=15592550),]
ref_a_sim<-ref_a_gene[!duplicated(ref_a_gene$Amplicon),]  # 2nd version
ref_a_sim<-ref_a_sim[which(ref_a_sim$Amplicon=='USP9Y' | ref_a_sim$Amplicon=='DDX3Y'),]

ref_gen_bc <- read.delim("chrY_reference_gene.txt")
ref_gen_a <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFa/hg19_refGene_AZFa_gene.txt")

# read Palindrome file
refpal_bc <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_pal_bc.txt")

# read marker location
Ma <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_markers_20201020.txt")




# creat reference palidrome

#tem<-ref[,c('Palindrome','start','end')]
refpal<-ref_bc[1,c('Palindrome','start','end')]
for (i in levels(ref_bc$Palindrome)){
  refpal[i,1]<-i
  refpal[i,2]<-min(ref_bc[ref_bc$Palindrome==i,]$start)
  refpal[i,3]<-max(ref_bc[ref_bc$Palindrome==i,]$end)}
refpal<-refpal[-(1:2),]
refpal$Palindrome<-factor(refpal$Palindrome)
refpal$color<-rep('firebrick2',nrow(refpal))




# creat reference amplicon with standard deletion type

#creat deletion location (nearby) for each type 
name<-'b2'
location<-refer[which(refer$Amplicon==name),]$start+((refer[which(refer$Amplicon==name),]$end-refer[which(refer$Amplicon==name),]$start)*0.2)
print(location)
(location-refer[which(refer$Amplicon==name),]$start)/(refer[which(refer$Amplicon==name),]$end-location)

# transform into amplicon location
for (i in 1:nrow(sample)){
  print(i)
  sampleI<-sample[i,]
  samp<-sampDelweb(i=1, refer=ref_bc, sample=sampleI)
  write.table(samp,paste0('reference_deletion_',unique(sampleI$sample_id),'.txt'),sep='\t',quote=FALSE,row.names = FALSE)}

i=5
sampleI<-sample[i,]
sampleI<-sample[which(sample$sample_id==sampleI$sample_id),]
repn<-nrow(sample[which(sample$sample_id==sampleI$sample_id),])
for (a in 1:repn){
  if (a==1){samp<-sampDelweb(i=a, refer=ref_bc, sample=sampleI)
  } else {
    ref_mul<-samp ; samp<-sampDelweb(i=a, refer=ref_mul, sample=sampleI)}}




# creat reference markers

ma <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_markers_20201020.txt")

library(tidyr)
#ma<-separate_rows(ma,location,sep=',') 
#ma$start<-gsub('.*:','',ma$location)
#ma$end<-gsub('.*-','',ma$start)
#ma$start<-gsub('-.*$','',ma$start
ma <- ma %>% separate(location, c('start','end'),'-')
write.table(ma,'C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_markers_20201020.txt',sep='\t',row.names = FALSE)




# creat AZFa deletion refer

sampleI <- sample_a[2,]
ref_a<-sampleI
ref_a$start=14517914; ref_a$end=16955848
ref_a$infor='region AZFa deletion reference'
ref_a
ref_a[2,]<-ref_a[1,]
ref_a$sample_id<-'AZFa'
ref_a$indel<-'NA'
rownames(ref_a)<-1:nrow(ref_a)



# creat reference gene (via annotation, not used here)

# a
ref_a_gene <-ref_a_gene[!duplicated(ref_a_gene$strans_start),]
ref_a_gene[7,5]<-16168838
ref_a_gene<-ref_a_gene[-8,]

ref_a_gene$coleur<-rep(as.factor('color'),nrow(ref_a_gene))
for (i in 1:nrow(ref_a_gene)){
  levels(ref_a_gene$coleur)<-c(levels(ref_a_gene$coleur),paste(ref_a_gene$Amplicon[i],toString(sprintf("'%s'", ref_a_gene$color[i])),sep = ' = '))
  ref_a_gene$coleur[i]<-paste(ref_a_gene$Amplicon[i],toString(sprintf("'%s'", ref_a_gene$color[i])),sep = ' = ')}

coleur<-paste(unlist(unique(ref_a_gene$coleur)), collapse=',') # add color to each feature
coleur 

ref_a_gene<-ref_a_gene[,c(2:5,11,13)]
colnames(ref_a_gene)[3:4]<-c('start','end')
ref_a_gene$col<-1:nrow(ref_a_gene)

write.table(ref_a_gene,'hg19_refGene_AZFa_gene.txt',sep='\t',row.names = FALSE)

# bc
tem<-ref_bc[,4:7]
tem$chr<-rep('chrY',nrow(tem))

nrow(ref_bc) #variants in region
gr <- GRanges(seqnames=Rle(tem$chr),
              ranges = IRanges(min(tem$start), end=max(tem$end)),
              strand = Rle(strand(c(rep("*", length(tem$chr)))))) 

genes(TxDb.Hsapiens.UCSC.hg19.knownGene) # all database build 37 / mm19
Genes<-subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr) # match

columns(org.Hs.eg.db)
symbol<-AnnotationDbi::select(org.Hs.eg.db, keys=Genes$gene_id, columns = c("SYMBOL","ENSEMBL","ACCNUM"), keytype="ENTREZID")
symbol<-a.symbol[!duplicated(a.symbol$SYMBOL), ]

refgen<-merge(Genes,symbol,by.x='gene_id',by.y='ENTREZID',all=TRUE)
write.table(refgen,'chrY_reference_gene.txt',sep="\t")













######################## create Track object ########################



###### IdeogramTrack (chromosome Ideogram plot)  * Reference * ######

ideoTrack_bc<- IdeogramTrack(genome='hg19', chromosome='chrY',fontcolor="black",
                             showBandId=TRUE, cex.bands=0.5) 
plotTracks(ideoTrack_bc,
           from=st[2,2], to=st[4,3], 
           showBandId=TRUE, fontsize=10,fontcolor="gray43",cex.bands=0.8)    

ideoTrack<-list()
ideoTrack[[1]]<-ideoTrack_a; ideoTrack[[2]]<-ideoTrack_bc
names(ideoTrack[[1]])<-'ideoTrack_a'; names(ideoTrack[[2]])<-'ideoTrack_bc'; 
rm(ideoTrack_a, ideoTrack_bc)





###### dataTrack (region zoom plot)  * Reference * ######


medposichrY_a<-(st[1,2]+st[1,3])/2-225000 #150000   # 250000 for negative reference AZFa
probes_a <- GRanges(seqnames="chrY", ranges=IRanges(start=c(st[1,2], medposichrY_a,  st[1,3]), end=c(st[1,2], medposichrY_a, st[1,3]+10)),hist=c(0,1,0))
dataTrack_a1 <- DataTrack(range=probes_a, genome="hg19", type="l", chromosome="chrY", size=0.2, col='lightskyblue')


medposichrY_a<-(st[1,2]+st[1,3])/2-100000
probes_a <- GRanges(seqnames="chrY", ranges=IRanges(start=c(st[1,2], medposichrY_a1,  st[1,3]), end=c(st[1,2], medposichrY_a1, st[1,3]+10)),hist=c(0,1,0))
dataTrack_a2 <- DataTrack(range=probes_a, genome="hg19", type="l", chromosome="chrY", size=0.2, col='lightskyblue')


medposichrY_bc2<-(st[2,2]+st[4,3])/2-3500000
probes_bc2 <- GRanges(seqnames="chrY", ranges=IRanges(start=c(st[2,2], medposichrY_bc2,  st[4,3]), end=c(st[2,2], medposichrY_bc2, st[4,3]+10)),hist=c(0,1,0))
dataTrack[[4]] <- DataTrack(range=probes_bc2, genome="hg19", type="l", chromosome="chrY", size=0.2, col='lightskyblue')
names(dataTrack[[4]])<-'dataTrack_bc2'

medposichrY_abc<-(st[1,2]+st[4,3])/2-1950000
probes_abc <- GRanges(seqnames="chrY", ranges=IRanges(start=c(st[1,2], medposichrY_abc,  st[4,3]), end=c(st[1,2], medposichrY_abc, st[4,3]+10)),hist=c(0,1,0))
dataTrack_abc <- DataTrack(range=probes_abc, genome="hg19", type="l", chromosome="chrY", size=0.2, col='lightskyblue')


dataTrack <- list(dataTrack_a1, dataTrack_a2, dataTrack_bc, dataTrack_abc)
dataTrack[[1]]<-dataTrack_a1
names(dataTrack[[1]]) <- 'dataTrack_a1'; names(dataTrack[[2]])<-'dataTrack_a2'; names(dataTrack[[3]])<-'dataTrack_bc';names(dataTrack[[4]])<-'dataTrack_abc'






###### (Genome) AxisTrack (genome location plot)  * Reference * ######

axisTrack_a <- GenomeAxisTrack(range=IRanges(start=st[1,2], end=st[1,3]), names='AZFa', cex=1,lwd=2) 

axisTrack_bc <- GenomeAxisTrack(range=IRanges(start=st[2,2], end=st[4,3]), names='AZFbc', cex=1,lwd=2)

axisTrack_abc <- GenomeAxisTrack(range=IRanges(start=st[1,2], end=st[4,3]), names='AZFabc', cex=1,lwd=2) 
plotTracks (axisTrack_bc, from=st[2,2], to=st[4,3],labelPos="below",add53=TRUE) 



# plot ideo + axis (checking step)
trackList_bc<-list(ideoTrack_bc,dataTrack_bc, axisTrack_bc)
plotTracks(trackList_bc,extend.right=0.1,size=0.2,
               showBandId=TRUE,fontsize=10,fontcolor="gray43",cex.bands=0.8,
               labelPos="below",add53=TRUE)



axisTrack<-list()
axisTrack[[1]]<-axisTrack_a; axisTrack[[2]]<-axisTrack_bc; axisTrack[[3]]<-axisTrack_abc
names(axisTrack[[1]])<-'axisTrack_a' ; names(axisTrack[[2]])<-'axisTrack_bc' ; names(axisTrack[[3]])<-'axisTrack_abc'
rm(axisTrack_a, axisTrack_bc, axisTrack_abc)






###### AnnotationTrack (annotation objects plot)  * Reference * ######
library(stringi)
neme<-'区域'
stri_enc_mark(name)
name <- stri_encode(name, "", "UTF-8")
all(stri_enc_isutf8(name))


### 1. AZFabc region  
regTrack_a<-AnnotationTrack(start=st[1,]$start, end=st[1,]$end, chromosome='chrY', strand="*",
                             genome="hg19", col.title="black",name='AZFa Region', shape='smallArrow',
                             fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                             group=rep('AZFa Region',nrow(st[1,])))
feature(regTrack_a)<- as.character(st[1,]$region)



regTrack_bc<-AnnotationTrack(start=st[2:4,]$start, end=st[2:4,]$end, chromosome='chrY', strand="*",
                               genome="hg19", col.title="black",name=name, shape='smallArrow',
                               fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                               group=rep(name,nrow(st[2:4,])))
feature(regTrack_bc)<- as.character(st[2:4,]$region)



regTrack_abc <- AnnotationTrack(start=st$start, end=st$end, chromosome='chrY', strand="*",
                                genome="hg19", col.title="black",name='Region', shape='smallArrow',
                                fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                group=rep('Region',nrow(st)))
feature(regTrack_abc)<- as.character(st$region)


regTrack_anosamp<-AnnotationTrack(start=st[1,]$start, end=st[1,]$end, chromosome='chrY', strand="*",
                            genome="hg19", col.title="black",name='Sample', shape='smallArrow',
                            fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                            group=rep('Sample',nrow(st[1,])))
feature(regTrack_anosamp)<- as.character(st[1,]$region)

regTrack<-list()
regTrack[[1]]<-regTrack_a; names(regTrack[[1]])<-'regTrack_a'
regTrack[[2]]<-regTrack_bc; names(regTrack[[2]])<-'regTrack_bc'
regTrack[[3]]<-regTrack_abc; names(regTrack[[3]])<-'regTrack_abc'
regTrack[[4]]<-regTrack_anosamp; names(regTrack[[4]])<-'regTrack_anosamp'
rm('regTrack_a','regTrack_bc','regTrack_abc', 'regTrack_anosamp')


# plot (checking step)
plotTracks(regTrack_bc, 
           featureAnnotation="feature",groupAnnotation="gr",
           cex.feature=0.6,collapse=FALSE,cex.group=0.8,
           AZFa='firebrick2',AZFb='firebrick2',AZFbc='firebrick2',AZFc='firebrick2',
)


# for color
for (i in 1:nrow(st)){
  st$coleur[i]<-paste(st$region[i],toString(sprintf("'%s'", st$color[i])),sep = '=')
}
coleur<-paste(unlist(st$coleur), collapse=',') # add color to each feature
coleur







### 2. markers 

ma<-Ma[which(Ma$region=='AZFb' | Ma$region=='AZFb+c'| Ma$region=='AZFc'),]
#ma<-ma[!duplicated(ma$marker),
markTrack_bc<- AnnotationTrack(start=ma$start, end=ma$end, chromosome='chrY', strand="*",
                            genome="hg19", col.title="black",name='Marker', shape='smallArrow',
                            fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                            group=rep('Marker',nrow(ma)))
feature(markTrack_bc)<- as.character(ma$marker)



ma<-Ma[which(Ma$region=='AZFa'),]
markTrack_a<- AnnotationTrack(start=ma$start, end=ma$end, chromosome='chrY', strand="*",
                                genome="hg19", col.title="black",name='Marker', shape='smallArrow',  
                                fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                group=rep('Marker',nrow(ma)))
feature(markTrack_a)<- as.character(ma$marker)



ma<-Ma
markTrack_abc<- AnnotationTrack(start=ma$start, end=ma$end, chromosome='chrY', strand="*",
                               genome="hg19", col.title="black",name='Marker', shape='smallArrow',
                               fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                               group=rep('Marker',nrow(ma)))
feature(markTrack_abc)<- as.character(ma$marker)



markTrack<-list()
markTrack[[1]]<-markTrack_a; markTrack[[2]]<-markTrack_bc; markTrack[[3]]<-markTrack_abc
names(markTrack[[1]])<-'markTrack_a'; names(markTrack[[2]])<-'markTrack_bc';names(markTrack[[3]])<-'markTrack_abc';



# plot (checking step)
plotTracks(markTrack_bc, 
           featureAnnotation="feature",groupAnnotation="gr",
           cex.feature=0.6,collapse=FALSE,cex.group=0.8,rotation.item=90,
           AZFa='firebrick2',AZFb='firebrick2',AZFbc='firebrick2',AZFc='firebrick2',
)










### 3. palindromes  

palTrack_bc <- AnnotationTrack(start=refpal$start, end=refpal$end, chromosome='chrY', strand="*",
                             genome="hg19", col.title="black",name='Palindromes', shape='smallArrow',
                             fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                             group=rep('Palindromes',nrow(refpal)))
feature(palTrack_bc)<- as.character(refpal$Palindrome)


# plot (checking step)
plotTracks(palTrack_bc, 
           featureAnnotation="feature",groupAnnotation="gr",
           cex.feature=0.6,collapse=FALSE,cex.group=0.8,
           P1='firebrick2',P2='firebrick2',P3='firebrick2',P4='firebrick2',P5='firebrick2',
)


# for color
for (i in 1:nrow(refpal)){
  refpal$coleur[i]<-paste(refpal$Palindrome[i],toString(sprintf("'%s'", refpal$color[i])),sep = '=')
}
coleur<-paste(unlist(refpal$coleur), collapse=',') # add color to each feature
coleur







### 4. amplicons  

## creat reference amplicon

ampTrack_bc <- AnnotationTrack(start=ref_amp_bc$start, end=ref_amp_bc$end, chromosome='chrY', strand=ref_amp_bc$str,
                          genome="hg19", col.title="black",name='Reference',
                          fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                          group=rep('Reference',nrow(ref_amp_bc)))

feature(ampTrack_bc)<- as.character(ref_amp_bc$Amplicon)

# plot (checking step)
plotTracks(ampTrack_bc, groupAnnotation="gr",
           featureAnnotation="feature",
           cex.feature=0.6,collapse=FALSE,
           yel3 = 'lightgoldenrod2',yel4 = 'lightgoldenrod2',b5 = 'deepskyblue3',u0 = 'navajowhite2',b6 = 'deepskyblue3',u1 = 'navajowhite2',b1 = 'deepskyblue3',t1 = 'lightseagreen',u2 = 'navajowhite2',t2 = 'lightseagreen',b2 = 'deepskyblue3',u3 = 'navajowhite2',g1 = 'darkseagreen3',r1 = 'firebrick',r2 = 'firebrick',Gr1 = 'gray99',b3 = 'deepskyblue3',yel1 = 'lightgoldenrod2',g2 = 'darkseagreen3',r3 = 'firebrick',r4 = 'firebrick',g3 = 'darkseagreen3',yel2 = 'lightgoldenrod2',b4 = 'deepskyblue3',Gr2 = 'gray99'
)


# for color
for (i in 1:nrow(ref)){
  levels(ref$coleur)<-c(levels(ref$coleur),paste(ref$Amplicon[i],toString(sprintf("'%s'", ref$color[i])),sep = ' = '))
  ref$coleur[i]<-paste(ref$Amplicon[i],toString(sprintf("'%s'", ref$color[i])),sep = ' = ')}
coleur<-paste(unlist(ref$coleur), collapse=',') # add color to each feature
coleur 


## creat reference + no sample

nosampTrack_bc <- AnnotationTrack(start=ref_amp_bc$start, end=ref_amp_bc$end, chromosome='chrY', strand=ref_amp_bc$str,
                               genome="hg19", col.title="black",name='Sample',
                               fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                               group=rep('Sample',nrow(ref_amp_bc)))

feature(nosampTrack_bc)<- as.character(ref_amp_bc$Amplicon)
names(nosampTrack_bc)<-'nosample'


ampTrack<-list()
ampTrack[[1]]<-ampTrack_bc; names(ampTrack[[1]])<-'ampTrack_bc';
ampTrack[[2]]<-nosampTrack_bc; names(ampTrack[[1]])<-'nosampTrack_bc'



## creat reference amplicon + reference deletion type
# method 1
setwd("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_deletion")

ampTrack_del<-list()

for (i in 1:length(list.files())){
  
file<-list.files()[i]
del <- read.delim(file)
name<- gsub('reference_deletion_','',file)
name<- gsub('.txt','',name)

tem <- AnnotationTrack(start=del$start, end=del$end, chromosome='chrY', strand=del$str,
                               genome="hg19", col.title="black",name=name,
                               fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                               group=rep(name,nrow(del)))
feature(tem)<- as.character(del$Amplicon)
ampTrack_del[[i]]<-tem
names(ampTrack_del[[i]])<-paste0('ampTrack_del_',name)
}


# method 2
Del <- read.csv("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_deletion_all_gcq_2020107.csv")
Del<-Del[which(Del$sample_id=='AZFa'|Del$sample_id=='AZFa_partial'),]
Del<-Del[order(Del$start),]
rownames(Del)<-1:nrow(Del)

rm(ampTrack_del)
ampTrack_del<-list()


for (i in 1:nrow(Del)){

  del <- Del[i,]
  #del <- Del[which(Del$sample_id=='b2/b3'),]
  name<-toString(unique(del$sample_id))
  len<-del$end-del$start
  name<-paste0(name,'(',round(len/1000000,digits=2),'Mb)')
  #name<-'b2/b3(2.23Mb)'
  
  tem <- AnnotationTrack(start=del$start, end=del$end, chromosome='chrY', strand='*',
                         genome="hg19", col.title="black",name=name,
                         fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                         group=rep(name,nrow(del)))
  
  #feature(tem)<- as.character(del$Amplicon)
  feature(tem)<-'缺失'
  #feature(tem)<-'Deletion'
  ampTrack_del[[i]]<-tem
  names(ampTrack_del[[i]])<-paste0('ampTrack_del_',name)
  
}

rm(Del,del)
#ampTrack_del[[19]]<-NULL












### 5. gene Tracks 

#ref$gene_res<-gsub(",.*","",ref$gene)
write.table(ref,'AZF_reference_20190620.txt',sep="\t",row.names = FALSE)

geneTrack_a <- AnnotationTrack (start=ref_a_gene$start, end=ref_a_gene$end, chromosome='chrY', strand=ref_a_gene$str,
                                genome="hg19", col.title="black", name='Gene',
                                # shape='box',
                                fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                group=rep('Gene',nrow(ref_a_gene)))
feature(geneTrack_a)<- as.character(ref_a_gene$Amplicon)




geneTrack_a_sim <- AnnotationTrack (start=ref_a_sim$start, end=ref_a_sim$end, chromosome='chrY', strand=ref_a_sim$str,
                                    genome="hg19", col.title="black", name='AZFa',
                                    fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                    group=rep('Gene',nrow(ref_a_sim)))
feature(geneTrack_a_sim)<- as.character(ref_a_sim$Amplicon)


geneTrack_bc <- AnnotationTrack(start=gene$start, end=gene$end, chromosome='chrY', strand=gene$strand,
                             genome="hg19", col.title="black",name='Genes', shape='box',
                             fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                             #group=rep('Genes',nrow(ref)),
                             group=gene$SYMBOL)
feature(geneTrack_bc)<- as.character(gene$SYMBOL)



geneTrack<-list()
geneTrack[[1]]<-geneTrack_a; geneTrack[[2]]<-geneTrack_a_sim
names(geneTrack[[1]])<-'geneTrack_a'; names(geneTrack[[2]])<-'geneTrack_a_sim'
rm(geneTrack_a, geneTrack_a_sim)






# plot (checking step)
plotTracks(geneTrack_a, featureAnnotation="feature",
           groupAnnotation="gr",
           cex.feature=0.6, collapse=FALSE,
           just.group="above",showFeatureId=FALSE,
           GYG2P1 = 'darkseagreen3',TTTY15 = 'navajowhite2',USP9Y = 'lightseagreen',DDX3Y = 'deepskyblue3',UTY = 'lightgoldenrod2',TMSB4Y = 'firebrick',VCY = 'lightblue3',VCY = 'lightblue3',NLGN4Y = 'paleturquoise4',
)



plotTracks(geneTrack_bc, 
           #featureAnnotation="feature",
           groupAnnotation="gr",rotation.item=90,cex.group=1,
           cex.feature=0.6,collapse=FALSE,
           #CDY2A='lightgoldenrodyellow',HSFY1='slategray2',BPY2='lightyellow',RBMY1E='slategray2',PRY2='lightcyan2',RBMY1J='lightyellow',RBMY1J='lightcyan2',BPY2='slategray2',BPY2='honeydew2',DAZ2='pink',DAZ2='gray99',DAZ2='lightgoldenrodyellow',BPY2B='honeydew2',BPY2B='pink',DAZ4='pink',DAZ4='honeydew2',BPY2B='lightgoldenrodyellow'
)


# for color
#ref$gene_res<-gsub(",.*","",ref$gene)

for (i in 1:nrow(ref)){
  if (ref$gene_res[i]!=""){
    levels(ref$coleur)<-c(levels(ref$coleur),paste(ref$gene_res[i],toString(sprintf("'%s'", ref$color_gene[i])),sep ='='))
    ref$coleur_gene[i]<-paste(ref$gene_res[i],toString(sprintf("'%s'", ref$color_gene[i])),sep ='=')
  }}
coleur_gene<-paste(unlist(unique(na.omit(ref$coleur_gene))), collapse=',') # add color to each feature
coleur_gene












### 6. Palindrome + Amplicon + Gene

trackList_bc<-list(palTrack_bc, ampTrack[[1]], geneTrack_bc)

# plot (checking step)
plotTracks(trackList_bc,extend.right=0.1,size=0.2,
           featureAnnotation="feature",groupAnnotation="gr",
           cex.feature=0.6,collapse=FALSE,showTitle=FALSE,
           transcriptAnnotation="symbol",
           P1='firebrick2',P2='firebrick2',P3='firebrick2',P4='firebrick2',P5='firebrick2',
           yel3 = 'lightgoldenrod2',yel4 = 'lightgoldenrod2',b5 = 'deepskyblue3',u0 = 'navajowhite2',b6 = 'deepskyblue3',u1 = 'navajowhite2',b1 = 'deepskyblue3',t1 = 'lightseagreen',u2 = 'navajowhite2',t2 = 'lightseagreen',b2 = 'deepskyblue3',u3 = 'navajowhite2',g1 = 'darkseagreen3',r1 = 'firebrick',r2 = 'firebrick',Gr1 = 'gray99',b3 = 'deepskyblue3',yel1 = 'lightgoldenrod2',g2 = 'darkseagreen3',r3 = 'firebrick',r4 = 'firebrick',g3 = 'darkseagreen3',yel2 = 'lightgoldenrod2',b4 = 'deepskyblue3',Gr2 = 'gray99'
)










### 7. seperation line 

ma<-Ma[which(Ma$region=='AZFb' | Ma$region=='AZFb+c'| Ma$region=='AZFc'),]
ma<-Ma[which(Ma$region=='AZFa'),]


probes_bc <- GRanges(seqnames="chrY", ranges=IRanges(start=c(min(ma$start),max(ma$start)), end=c(min(ma$end),max(ma$end))),hist=c(0,0))
sepTrack_bc <- DataTrack(range=probes_bc, genome="hg19", type="l", chromosome="chrY", size=0.5, col='lightskyblue')

probes_a <- GRanges(seqnames="chrY", ranges=IRanges(start=c(min(ma$start),max(ma$start)), end=c(min(ma$end),max(ma$end))),hist=c(0,0))
sepTrack_a <- DataTrack(range=probes_a, genome="hg19", type="l", chromosome="chrY", size=0.5, col='lightskyblue')


sepTrack<-list()
sepTrack[[1]]<-sepTrack_a; sepTrack[[2]]<-sepTrack_bc
names(sepTrack[[1]])<-'sepTrack_a'; names(sepTrack[[2]])<-'sepTrack_bc'








########################   plot samples   ######################## 

###### AnnotationTrack plot * Sample bc * ######



for (i in 40:nrow(sample)){

  if (sample$Indel[i]=="Deletion"){
    
  samp<-sampDel(i, refer=ref, sample=sample)

  sampTrack_bc <- AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                            genome="hg19", col.title="black", name=paste('Sample',i),
                            fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                            group=rep(paste('Sample',i),nrow(samp)))

   feature(sampTrack_bc)<- as.character(samp$Amplicon)

   plotTracks(sampTrack_bc, groupAnnotation="gr",
           featureAnnotation="feature",
           cex.feature=0.6,collapse=FALSE,
           yel3 = 'lightgoldenrod2',yel4 = 'lightgoldenrod2',b5 = 'deepskyblue3',u0 = 'navajowhite2',b6 = 'deepskyblue3',u1 = 'navajowhite2',b1 = 'deepskyblue3',t1 = 'lightseagreen',u2 = 'navajowhite2',t2 = 'lightseagreen',b2 = 'deepskyblue3',u3 = 'navajowhite2',g1 = 'darkseagreen3',r1 = 'firebrick',r2 = 'firebrick',Gr1 = 'gray99',b3 = 'deepskyblue3',yel1 = 'lightgoldenrod2',g2 = 'darkseagreen3',r3 = 'firebrick',r4 = 'firebrick',g3 = 'darkseagreen3',yel2 = 'lightgoldenrod2',b4 = 'deepskyblue3',Gr2 = 'gray99',
           deletion='gray98'
           )
  }}







### AnnotationTrack plot * Sample a * ######


samp<-sample[i,]


sampTrack_a <- AnnotationTrack(start=samp$start, end=samp$end, chromosome='chrY', strand=samp$str,
                               genome="hg19", col.title="black", name=paste('Sample',i),
                               fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                               group=rep(paste('Sample',i),nrow(samp)))

feature(sampTrack_a)<- as.character(samp$Indel)

plotTracks(sampTrack_a, groupAnnotation="gr",
           featureAnnotation="feature",
           cex.feature=0.6,collapse=FALSE,
           Deletion='gray98'
)


















########################    other informations   ######################## 



###### GeneRegionTrack ######


# (not use here) build GeneRegionTrack objects from TxDbs (gene model information from online sources) 

samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite", package="GenomicFeatures")
txdb <- loadDb(samplefile)

GeneRegionTrack(txdb)
txTr <- GeneRegionTrack(txdb, chromosome="chrY", start=19567683, end=20063140)
feature(txTr)
plotTracks(txTr)


# build GeneRegionTrack objects from ensembl

bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", name="ENSEMBL", symbol="TTTY15", biomart=bm)


for (i in 1:nrow(ref_a)){
  geneTrack_a <- BiomartGeneRegionTrack(genome="hg19", name="ENSEMBL", filter=list(with_refseq_mrna=TRUE),
                                        chromosome='chrY', start=min(ref_a$strans_start), end=max(ref_a$strans_end), biomart=bm) # name="AZFa"
  
  feature(geneTrack_a)<- as.character(ref_a$gene_ID)
  
  plotTracks(list(ideoTrack,axisTrack_a,geneTrack_a), 
             transcriptAnnotation="symbol", collapseTranscripts="longest")
  
  levels(ref_a$gene)<-c(levels(ref_a$gene), toString(unique(symbol(geneTrack_a))))
  ref_a$gene[i]=toString(unique(symbol(geneTrack_a)))
}








###### SequenceTrack ######

sTrack <- SequenceTrack(Hsapiens)
sTrack

plotTracks(sTrack, chromosome=1, from=20000, to=20050, add53=TRUE, complement=TRUE)
plotTracks(sTrack, chromosome=1, from=20000, to=20100)





### AlignmentsTrack

afrom <- 44945200
ato <- 44947200
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata", "snps.bam"), isPaired=TRUE)
plotTracks(c(alTrack, sTrack), chromosome="chr21", from=44946590, to=44946660)




# DataTrack ?

data(sample)
dTrack <- DataTrack(start=seq(min(st), max(ed), len=10), width=18000,
                    data=matrix(runif(40), nrow=4), genome="hg19", chromosome='chrY',
                    type="histogram", name="data track bar")
plotTracks(dTrack)






### Highlighting

ht <- HighlightTrack(trackList=list(atrack, grtrack, dtrack), start=c(26705000, 26720000), width=7000, chromosome=7)
plotTracks(list(itrack, gtrack, ht), from = lim[1], to = lim[2])









### UCSC data

from <- 25828122
to <- 25840710

knownGenes <- UcscTrack(genome="mm9", chromosome=chr, track="knownGene", from=from, to=to,
                        trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                        symbol="name", transcript="name", strand="strand", fill="#8282d2", name="UCSC Genes")

refGenes <- UcscTrack(genome="mm9", chromosome=chr, track="xenoRefGene", from=from, to=to,
                      trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                      symbol="name2", transcript="name", strand="strand", fill="#8282d2",
                      stacking="dense", name="Other RefSeq")

ensGenes <- UcscTrack(genome="mm9", chromosome=chr, track="ensGene", from=from, to=to,
                      trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                      symbol="name2", transcript="name", strand="strand", fill="#960000",
                      name="Ensembl Genes")

snpLocations <-  UcscTrack(genome="mm9", chromosome=chr, track="snp128", from=from, to=to,
                           trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",
                           feature="func", strand="strand", shape="box", stacking="dense", fill="black",
                           name="SNPs")

conservation <- UcscTrack(genome="mm9", chromosome=chr, track="Conservation", table="phyloP30wayPlacental",
                          from=from, to=to, trackType="DataTrack", start="start", end="end", data="score",
                          type="hist", window="auto", col.histogram="darkblue", fill.histogram="darkblue",
                          ylim=c(-3.7, 4), name="Conservation")

gcContent <- UcscTrack(genome="mm9", chromosome=chr, track="GC Percent", table="gc5Base",
                       from=from, to=to, trackType="DataTrack", start="start", end="end", data="score",
                       type="hist", window=-1, windowSize=1500, fill.histogram="black", col.histogram="black",
                       
                       ylim=c(30, 70), name="GC Percent")

axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome="mm9", chromosome=chr)

plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, ensGenes, cpgIslands,
                gcContent, conservation, snpLocations), from=from, to=to, showTitle=FALSE)