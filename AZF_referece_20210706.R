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

library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg19)



#########################################################################
#####################      AZFa & AZFbc plot        #####################
#########################################################################


load("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/AZF_plot_raport/AZF_plot_20201019/AZF_plot_Reference_raport_20201021_gcq.RData") # load Rdata reference





######################## Import data for Track input ########################

# read region
st<-read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_region_abc.txt")

# read genes
ref_a_gene <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFa/hg19_refGene_AZFa_gene.txt")
ref_a_gene<-ref_a_gene[which(ref_a_gene$end<=15592550),]
ref_a_sim<-ref_a_gene[!duplicated(ref_a_gene$Amplicon),]  # 2nd version
ref_a_sim<-ref_a_sim[which(ref_a_sim$Amplicon=='USP9Y' | ref_a_sim$Amplicon=='DDX3Y'),]

ref_gen_bc <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/chrY_reference_gene.txt")
ref_gen_a <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFa/hg19_refGene_AZFa_gene.txt")

# read marker location
Ma <- read.delim("C:/Users/Jingz/OneDrive/CarrierGene/2019.09.02-AZF.Plot (done)/Reference_AZFbc/reference_markers_20201020.txt")









######################## create Track object ########################



###### (Genome) AxisTrack (genome location plot)  * Reference * ######

axisTrack_a <- GenomeAxisTrack(range=IRanges(start=st[1,2], end=st[1,3]), names='AZFa', cex=1,lwd=2) 
axisTrack_bc <- GenomeAxisTrack(range=IRanges(start=st[2,2], end=st[4,3]), names='AZFbc', cex=1,lwd=2)
axisTrack_abc <- GenomeAxisTrack(range=IRanges(start=st[1,2], end=st[4,3]), names='AZFabc', cex=1,lwd=2) 



axisTrack<-list()
axisTrack[[1]]<-axisTrack_a; axisTrack[[2]]<-axisTrack_bc; axisTrack[[3]]<-axisTrack_abc
names(axisTrack[[1]])<-'axisTrack_a' ; names(axisTrack[[2]])<-'axisTrack_bc' ; names(axisTrack[[3]])<-'axisTrack_abc'
rm(axisTrack_a, axisTrack_bc, axisTrack_abc)






###### AnnotationTrack (annotation objects plot)  * Reference * ######
library(stringi)
name<-'区域'
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
rm(regTrack_a,regTrack_bc,regTrack_abc, regTrack_anosamp)











### 2. markers 

ma<-Ma[which(Ma$region=='AZFb' | Ma$region=='AZFb+c'| Ma$region=='AZFc'),]
markTrack_bc<- AnnotationTrack(start=ma$start_m, end=ma$end_m, chromosome='chrY', strand="*",
                            genome="hg19", col.title="black",name='Marker', shape='smallArrow',
                            fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                            group=rep('Marker',nrow(ma)))
feature(markTrack_bc)<- as.character(ma$marker)



ma<-Ma[which(Ma$region=='AZFa'),]
markTrack_a<- AnnotationTrack(start=ma$start_m, end=ma$end_m, chromosome='chrY', strand="*",
                                genome="hg19", col.title="black",name='Marker', shape='smallArrow',  
                                fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                                group=rep('Marker',nrow(ma)))
feature(markTrack_a)<- as.character(ma$marker)


markTrack<-list()
markTrack[[1]]<-markTrack_a; markTrack[[2]]<-markTrack_bc
names(markTrack[[1]])<-'markTrack_a'; names(markTrack[[2]])<-'markTrack_bc'

rm(markTrack_bc,markTrack_a)









### 3. amplicons  

## creat reference amplicon

ampTrack_bc <- AnnotationTrack(start=ref_amp_bc$start, end=ref_amp_bc$end, chromosome='chrY', strand=ref_amp_bc$str,
                          genome="hg19", col.title="black",name='Reference',
                          fontcolor.feature=1,cex.feature=1,region="grey",lty.grid="times",fontsize=15,
                          group=rep('Reference',nrow(ref_amp_bc)))

feature(ampTrack_bc)<- as.character(ref_amp_bc$Amplicon)






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
rm(ampTrack_bc,nosampTrack_bc)








### 4. gene Tracks 


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



geneTrack<-list()
geneTrack[[1]]<-geneTrack_a; geneTrack[[2]]<-geneTrack_a_sim
names(geneTrack[[1]])<-'geneTrack_a'; names(geneTrack[[2]])<-'geneTrack_a_sim'
rm(geneTrack_a, geneTrack_a_sim)








### creat trackList for plot
trackList_a<-list(ideoTrack[[1]], dataTrack[[1]], axisTrack[[1]], regTrack[[1]], markTrack[[1]], geneTrack[[1]], sepTrack[[1]])
trackList_bc1<-list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], sepTrack[[2]])
trackList_bc2<-list(ideoTrack[[2]], dataTrack[[4]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], sepTrack[[2]])

### creat neg ref bc
trackList_bc<-list(ideoTrack[[2]], dataTrack[[3]], axisTrack[[2]], regTrack[[2]], markTrack[[2]], ampTrack[[1]], sepTrack[[2]], ampTrack[[2]])






######################## save  ########################

rm(st, ref_a_gene, ref_a_sim, ref_gen_bc, ref_gen_a, Ma,ma, name)
save.image("AZF_plot_Reference_raport_202110.Rdata") # save new data into .Rdata reference




