library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")


library("DESeq2")
library("Matrix")
library("caret")


args = commandArgs(trailingOnly=TRUE)


convert_0_1vector<-function(X){
  if(all(X==0)){
    Y = X
  } else {
    Y<-(X-min(X))/(max(X)-min(X))
  }
  
  Y
}


##if the features are along columns, then you have to use mode = 'col'
convert_0_1matrix<-function(X,featureLoc){
  if(featureLoc=='col') { ## features are along column, we need normalize along column
    x = 2
    temp = apply(X,x,convert_0_1vector)
    ans = temp
  }else { ## features are along row, we need normalize along row
    x = 1
    temp = apply(X,x,convert_0_1vector)
    ans = t(temp)
  }
  ans
}



calculateFreq<-function(X){
  x<-unique(X)
  freq<-rep(0,length(x))
  for(i in 1:length(x)){
    freq[i]<-length(which(X==x[i]))
  }
  ans<-list(x,freq)
  ans
}




checkForNAsInMatrices<-function(A){
  ##check for NAs
  na_list<-c()
  for(i in 1:nrow(A)){
    ind1<-which(is.na(A[i,])==T)
    if(length(ind1)>0)  na_list<-c(na_list,paste0(i,",",ind1))
  }
  na_list
}




##

brafDEG<-read.csv(file = paste0("parallel_8_NvsBraf.like_GA_deGenes_SolnSize=273_2026-05-02.csv"), header = T)
braf_classifier<-readRDS(file = paste0("parallel_8_NvsBraf.like_GA_Classifier_SolnSize=273_2026-05-02.rds"))

rasDEG<-read.csv(file = paste0("parallel_8_NvsRas.like_GA_deGenes_SolnSize=986_2026-05-02.csv"), header = T)
ras_classifier<-readRDS(file = paste0("parallel_8_NvsRas.like_GA_Classifier_SolnSize=986_2026-05-02.rds"))

retDEG<-read.csv(file = paste0("parallel_8_Nvs8RET.fusion_GA_deGenes_SolnSize=58_2026-05-02.csv"), header = T)
ret_classifier<-readRDS(file = paste0("parallel_8_Nvs8RET.fusion_GA_Classifier_SolnSize=58_2026-05-02.rds"))

nvtDEG<-read.csv(file = paste0("parallel_9_NvsT_GA_deGenes_SolnSize=4382_2026-05-06.csv"), header = T)
nvt_classifier<-readRDS(file = paste0("parallel_9_NvsT_GA_Classifier_SolnSize=4382_2026-05-06.rds"))


tmp = union(brafDEG$deGenes,rasDEG$deGenes)
tmp2 = union(tmp,retDEG$deGenes)
all4DEG = union(tmp2,nvtDEG$deGenes)

sampleName = args[1]
sample<-Load10X_Spatial(st_path,filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")

sample <- SCTransform(sample, assay = "Spatial", verbose = FALSE)
sample <- RunPCA(sample, assay = "SCT", verbose = FALSE)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
sample <- FindClusters(sample, verbose = FALSE)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:30)
Idents(sample)<-sample@meta.data$seurat_clusters

process2<-as.matrix(GetAssayData(sample@assays$Spatial))
rSums<-rowSums(process2)
tbr<-which(rSums==0)
process3<-process2[-tbr,]

count_stand0 = process3

##first level

deg2<-intersect(all4DEG,rownames(count_stand0))

if(length(all4DEG)==length(deg2)){
  ##
  print("All Same. Nothing to do")
  count_stand1 = count_stand0
  ##
}else {
  #
  foundIndx<-which(all4DEG %in% rownames(count_stand0))
  notFound<-all4DEG[-foundIndx]
  
  toBeAdded<-matrix(0, nrow = length(notFound), ncol = ncol(count_stand0))
  rownames(toBeAdded)<-notFound
  colnames(toBeAdded)<-colnames(count_stand0)
  count_stand1 = rbind(count_stand0,toBeAdded)
  #
}


#sanity check
deg2<-intersect(all4DEG,rownames(count_stand1))
length(all4DEG)==length(deg2)
deg2Indx<-which(rownames(count_stand1) %in% deg2)


#GE2<-as.matrix(count_stand1[deg2Indx,])

GE2<-as.matrix(count_stand1)
dim(GE2)
GE3<-convert_0_1matrix(GE2,'row')

pcaInput<-t(GE3)

dim(pcaInput)

firstLevelInput<-as.data.frame(pcaInput)

rn<-rownames(firstLevelInput)
rn2=gsub(".1","-1",rn)

#set.seed(56)

braf_pred<-predict(braf_classifier,firstLevelInput)
braf_prob<-predict(braf_classifier,firstLevelInput,type="prob")
ras_pred<-predict(ras_classifier,firstLevelInput)
ras_prob<-predict(ras_classifier,firstLevelInput,type="prob")
ret_pred<-predict(ret_classifier,firstLevelInput)
ret_prob<-predict(ret_classifier,firstLevelInput,type="prob")
nvt_pred<-predict(nvt_classifier,firstLevelInput)
nvt_prob<-predict(nvt_classifier,firstLevelInput,type="prob")

Pred_firstLevelOutputProb=data.frame(Barcode=rn2)
Pred_firstLevelOutputProb$BRAF_braf = as.vector(braf_prob$Braf.like)
Pred_firstLevelOutputProb$BRAF_normal = as.vector(braf_prob$Normal)
Pred_firstLevelOutputProb$RAS_ras = as.vector(ras_prob$Ras.like)
Pred_firstLevelOutputProb$RAS_normal = as.vector(ras_prob$Normal)
Pred_firstLevelOutputProb$RET_ret = as.vector(ret_prob$RET.fusion)
Pred_firstLevelOutputProb$RET_normal = as.vector(ret_prob$Normal)
Pred_firstLevelOutputProb$Tumor = as.vector(nvt_prob$Tumor)
Pred_firstLevelOutputProb$Normal = as.vector(nvt_prob$Normal)

Pred_firstLevelOutputProb$Decision=""
#Pred_firstLevelOutputProb$Details=""


getLabel<-function(a){
  ans=NULL
  m = max(a)
  if(a[4]==m) {
    ans = "Other Tumor Type"
  }else if(a[1]==m) {
    ans = "Braf.like"
  }else if(a[2]==m) {
    ans = "Ras.like"
  }else if(a[3]==m) {
    ans = "RET.fusion"
  }else if(a[1]==m & a[2]==m){
    ans = "Braf.like-Ras.like"
  }else if(a[1]==m & a[3]==m){
    ans = "Braf.like-RET.fusion"
  }else if(a[2]==m & a[3]==m){
    ans = "Ras.like-RET.fusion"
  }else if(a[1]==m & a[2]==m & a[3]==m){ ## just for sake of completion. We do not hope this case will arise
    ans = "Braf.like-Ras.like-RET.fusion"
  }
  ans
}

## if Pred_firstLevelOutputProb$Normal[i]>0.5 decision = Normal
## if Pred_firstLevelOutputProb$Tumor[i]>0.5 decision = get label of majority from amongst the subtypes
## if Pred_firstLevelOutputProb$Normal[i]==0.5
##    if all 3 subtype normals > 0.5, then decision Normal
##    else lean on majority subtype tumor

for(i in 1:nrow(Pred_firstLevelOutputProb)){
 
  A = unlist(Pred_firstLevelOutputProb[i,])
  N = A[c(3,5,7,9)]
  T = A[c(2,4,6,8)]
  g=length(which(N>0.5)) ##how many greater than 0.5
  g1=length(which(T>0.5)) ##how many greater than 0.5
  
  ##annotating normals
  if(Pred_firstLevelOutputProb$Normal[i]>0.5){
   Pred_firstLevelOutputProb$Decision[i] = "Normal"
   
  }else if(Pred_firstLevelOutputProb$Tumor[i]>0.5){
   ## now decide the subtype
   lbl = getLabel(T)
   Pred_firstLevelOutputProb$Decision[i] = lbl
   
  }else if(Pred_firstLevelOutputProb$Normal[i]==0.5){

    if(Pred_firstLevelOutputProb$BRAF_normal[i]>0.5 & Pred_firstLevelOutputProb$RAS_normal[i]>0.5 & Pred_firstLevelOutputProb$RET_normal[i]>0.5){
      Pred_firstLevelOutputProb$Decision[i] = "Normal"
    
    }else {
      
      lbl = getLabel(T)
      Pred_firstLevelOutputProb$Decision[i] = lbl
      
    }
  }
  
}



#table(Pred_firstLevelOutputProb$Decision)


write.csv(Pred_firstLevelOutputProb, file = paste0(sampleName,"_RFGA_FinalAnnotation.csv"), row.names = F)

#### drawing the figure #####

sanity = all(Pred_firstLevelOutputProb$Barcode==rownames(sample@meta.data))
sanity

if(sanity){
  ##do nothing
  data = Pred_firstLevelOutputProb
}else{
  indx1 = c()
  for(i in 1:length(rownames(sample@meta.data))){
    tmp = which(Pred_firstLevelOutputProb$Barcode==rownames(sample@meta.data)[i])
    indx1 = c(indx1,tmp)
  }
  data = Pred_firstLevelOutputProb[indx1,]
  all(data$Barcode==rownames(sample@meta.data))
  
}



# X = table(data$Decision)
# 
# names(X)
# 
# X/nrow(data)
# 
# sum(X/nrow(data))
# 
# 100*X/nrow(data)
# 
# round(100*X/nrow(data), digits = 3)
# 
# Y = as.data.frame(round(100*X/nrow(data), digits = 3))
# colnames(Y)=c("Type","Perc")
# 
# write.csv(Y, file = paste0("spatialData_firstLevel/",sampleName,"_RFGA_usingSubTypeClassifiers_Freq.csv"), row.names = F)

sample@meta.data$RFGA_prediction = data$Decision

sample@meta.data$RFGA_prediction2 = data$Decision
indx2 = which(sample@meta.data$RFGA_prediction2=="Normal")
sample@meta.data$RFGA_prediction2[-indx2] = "Tumor"


### finding pt.size.factor value by trial. Turns out 3000, is the correct value.
SpatialDimPlot(sample, label = FALSE, pt.size.factor = 3000, label.size = 3)

mycol2 = c("Tumor"="red", "Normal"="darkgreen")
mycol = c("Other Tumor Type"="red", "Normal"="darkgreen", "Braf.like"="yellow", "Ras.like"="cyan","Braf.like-Ras.like"="gold", "Braf.like-RET.fusion"="orange", "Ras.like-RET.fusion"="magenta","Braf.like-Ras.like-RET.fusion"="pink", "RET.fusion"="navyblue")

tiff(filename = paste0(sampleName,"_RFGA_Normal_Tumor_Annotation.tiff"), width = 4, height = 4, units = 'in', res = 200)
SpatialDimPlot(sample, label = FALSE, label.size = 2, pt.size.factor = 3, cols = mycol2, group.by = 'RFGA_prediction2')
dev.off()

tiff(filename = paste0(sampleName,"_RFGA_Tumor_SubType_Annotation.tiff"), width = 4, height = 4, units = 'in', res = 200)
SpatialDimPlot(sample, label = FALSE, label.size = 2, pt.size.factor = 3, cols = mycol, group.by = 'RFGA_prediction')
dev.off()


