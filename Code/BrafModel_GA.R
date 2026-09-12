### 09/04/2026
### implementing Genetic Algorithm: between Tumor SubType and normal #### 


# this time we implement parallel processing
# also weightage calculate is 1/3 for each score

### Same as 7_runGA except here we add more components to the objective function


library(dplyr)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")


#library("DESeq2")
library("Matrix")
library("caret")
library("ranger")
library(parallel)
library(doParallel)


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

### given a Truth and Prediction, compare each element by id(barcode). If match reward +1, if error reward -1. return fraction as score 
spatialValidation<-function(Truth,Prediction){ ##size of Prediction is more than that of Truth. First delete the extra spots then reorder and compare
  values = list()
  
  ## remove the extra spots from Prediction
  index2=which(Prediction$Barcode %in% Truth$X)
  pred = Prediction[index2,]
  
  
  ## separate out the tumors and normals
  tIndx1 = which(Truth$value=="Tumor")
  Truth_tumor = Truth[tIndx1,]
  Truth_normal = Truth[-tIndx1,]
  
  tIndx2 = which(pred$Predicted=="Normal")
  pred_tumor = pred[-tIndx2,]
  pred_normal = pred[tIndx2,]
  
  
  correctly_classified_truthNormal = length(intersect(Truth_normal$X,pred_normal$Barcode))
  incorrectly_classified_truthNormal = length(intersect(Truth_normal$X,pred_tumor$Barcode))
  correctly_classified_truthTumor = length(intersect(Truth_tumor$X,pred_tumor$Barcode))
  incorrectly_classified_truthTumor = length(intersect(Truth_tumor$X,pred_normal$Barcode))
  
  ## if something is predicted to be braf/ras/ret, it should not fall on a normal spot
  
  print(paste0("Sanity: ",(correctly_classified_truthTumor+correctly_classified_truthNormal+incorrectly_classified_truthNormal+incorrectly_classified_truthTumor)==nrow(Truth)))
  
  epsilon=0.000000
  
  print(paste0("#correctly_classified_truthNormal: ",correctly_classified_truthNormal," #incorrectly_classified_truthNormal: ",incorrectly_classified_truthNormal," #correctly_classified_truthTumor: ",correctly_classified_truthTumor," #incorrectly_classified_truthTumor: ",incorrectly_classified_truthTumor))
  
  values[[1]] = correctly_classified_truthNormal/(nrow(Truth_normal)+epsilon)
  values[[2]] = incorrectly_classified_truthNormal/(nrow(Truth_normal)+epsilon)
  values[[3]] = correctly_classified_truthTumor/(nrow(Truth_tumor)+epsilon)
  values[[4]] = incorrectly_classified_truthTumor/(nrow(Truth_tumor)+epsilon)
  
  values
}


createChromosome<-function(UnusedGenes,size){
  ans = list()
  indx<-sample(c(1:length(UnusedGenes)),size)
  chrom<-UnusedGenes[indx]
  ans[[1]] = chrom
  ans[[2]] = UnusedGenes[-indx]
  ans
}


##
checkChromosomeSimilarity<-function(chrom1,chrom2){
  ans<-F
  if(length(chrom1)!=length(chrom2)){
    ans=F
  } else {
    common<-intersect(chrom1,chrom2)
    if(length(common)==length(chrom1) & length(chrom2)==length(chrom1))  ans<-T
  }
  ans
}



areThereDuplicates<-function(chrom){
  # indx<-which(duplicated(chrom))
  # ans = F
  # if(length(indx)) ans = T
  # ans
  any(duplicated(chrom))
}


## return -1 if not present, else the index
isAlreadyPresent<-function(listOfChromosomes,chrom){
  ans=-1
  if(length(listOfChromosomes)==0){
    ans=-1
  } else {
    for(i in 1:length(listOfChromosomes)){
      if(checkChromosomeSimilarity(chrom,listOfChromosomes[[i]])){
        ans=i
        i<-lengths(listOfChromosomes)+1 ##break out of the loop
      }
    }
  }
  ans
}

## return -999999 if not present, else the stored score
getChromosomeFitnessScoreFromList<-function(listOfChromosomes,listOfChromosomeScores,chrom){
  indx = isAlreadyPresent(listOfChromosomes,chrom)
  ans = -999999
  if(indx==-1){
    ##
  }else{
    ans = listOfChromosomeScores[[indx]]
  }
  ans
}


## different length chromosomes. random separate crossover point on each chrom. Swap the tails
singlePntcrossOver<-function(chrom1, chrom2){
  
  cpoint1<-sample(c(2:length(chrom1)),1) ##cannot start from 1, because the break will happen before cpoint
  temp11<-chrom1[c(1:(cpoint1-1))]
  temp12<-chrom1[c(cpoint1:length(chrom1))]
  
  cpoint2<-sample(c(2:length(chrom2)),1) ##cannot start from 1, because the break will happen before cpoint
  temp21<-chrom2[c(1:(cpoint2-1))]
  temp22<-chrom2[c(cpoint2:length(chrom2))]
  
  offspring<-list()
  offspring[[1]]<-c(temp11,temp22)
  offspring[[2]]<-c(temp21,temp12)
  offspring
}

#chrom1<-c(1,2,3,4,5,6,7,8,9,10)
#chrom2<-c(11,22,33,44,55,66,77,88,99,100)

## different length chromosomes. random separate crossover points on each chrom. Swap the intermediate regions
doublePntcrossOver<-function(chrom1, chrom2){
  
  repeat {
    cpoint11<-sample(c(2:(length(chrom1)-0)),1)
    cpoint12<-sample(c(2:(length(chrom1)-0)),1)
    ##region after smaller but before larger
    if(cpoint11<cpoint12) {
      smaller=cpoint11
      larger=cpoint12
    } else if(cpoint11>cpoint12) {
      smaller=cpoint12
      larger=cpoint11
    }
    
    if(cpoint11!=cpoint12)  break
  }
  
  temp11<-chrom1[c(1:(smaller-1))]
  temp12<-chrom1[c(smaller:(larger-1))]
  temp13<-chrom1[c(larger:length(chrom1))]
  
  repeat {
    cpoint21<-sample(c(2:(length(chrom2)-0)),1)
    cpoint22<-sample(c(2:(length(chrom2)-0)),1)
    ##region after smaller but before larger
    if(cpoint21<cpoint22) {
      smaller=cpoint21
      larger=cpoint22
    } else if(cpoint21>cpoint22) {
      smaller=cpoint22
      larger=cpoint21
    }
    
    if(cpoint21!=cpoint22) break
  }
  
  temp21<-chrom2[c(1:(smaller-1))]
  temp22<-chrom2[c(smaller:(larger-1))]
  temp23<-chrom2[c(larger:length(chrom2))]
  
  offspring<-list()
  offspring[[1]]<-c(temp11,temp22,temp13)
  offspring[[2]]<-c(temp21,temp12,temp23)
  offspring
  
}



##
mutate<-function(yetToBeUsedGenes,chrom){
  
  mutationPoint<-sample(c(1:length(chrom)),1)
  e<-sample(c(1:length(yetToBeUsedGenes)),1)
  chrom[mutationPoint]<-yetToBeUsedGenes[e]
  
  answers<-list()
  answers[[1]]<-chrom
  answers[[2]]<-yetToBeUsedGenes[-e]
  answers
}





rouletteSelection<-function(partialSums,S){
  ans=-1
  
  print(paste0("S: ",S))
  
  if(S<=partialSums[1]) {
    ans = 1
  }else {
    for(i in 2:length(partialSums)){
      if(S>partialSums[i-1] && S<=partialSums[i]){
        ans = i
        i = length(partialSums)+999 ##break
      }
    }
  }
  ans
}



##
getChromosomeFitnessScore<-function(chrom){
  
  deGenes<-chrom
  ###deGenes = c("CDH13","TLE1","NDUFA13","ZNF471")
  degenesIndx<-which(rownames(count_stand2) %in% deGenes)
  GE = count_stand2[degenesIndx,]
  pcaInput<-t(GE)
  rownames(pcaInput)<-colnames(GE)
  colnames(pcaInput)<-rownames(GE)
  
  data3<-as.data.frame(pcaInput)
  ##sanity check
  all(rownames(data3)==ww3$ID)
  data3$Type = factor(ww3$Type2)
  
  calculateFreq(data3$Type)
  
  set.seed(107)
  inTrain <- createDataPartition(
    y = data3$Type,
    ## the outcome data are needed
    p = .75,
    ## The percentage of data in the
    ## training set
    list = FALSE
  )
  
  str(inTrain)
  
  training <- data3[inTrain,]
  testing  <- data3[-inTrain,]
  
  formula = paste(colnames(data3)[ncol(data3)], '~.', sep = '')
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,     # number of folds
                             repeats = 5,     # repeated 5 times
                             classProbs = T,  # remember class probabilities
                             allowParallel = TRUE)    # Explicitly allow caret to use the registered backend
  
  trainingVal = as.matrix(training[,-ncol(training)])
  dim(trainingVal)
  trainingVal2 = convert_0_1matrix(trainingVal, 'col')
  training2 = as.data.frame(trainingVal2)
  training2$Type = training$Type
  
  
  testingVal = as.matrix(testing[,-ncol(testing)])
  dim(testingVal)
  testingVal2 = convert_0_1matrix(testingVal, 'col')
  testing2 = as.data.frame(testingVal2)
  testing2$Type = testing$Type
  
  print(paste0("Training rf2"))
  
  rf2_y = train(as.formula(formula), method = "ranger", data = training2, trControl=fitControl)
  
  testing_pred <- predict(rf2_y, testing2)
  
  zz1<-confusionMatrix(testing_pred, as.factor(testing2$Type))
  
  zz1[[4]][11] ##use (balanced accuracy = sensitivity/2 + specificity/2)
  
  print(paste0("TCGA BalancedAccuracy: ",zz1[[4]][11]))
  
  training_pred <- predict(rf2_y, training2)
  
  zz2<-confusionMatrix(training_pred, as.factor(training2$Type))
  
  zz2[[4]][11] ##balanced accuracy
  
  
  tempIndx = which(colnames(wholePipelineTesting) %in% deGenes)
  wholePipelineTesting1 = wholePipelineTesting[,tempIndx]
  
  pred_spat<-predict(rf2_y,wholePipelineTesting1)
  
  rf2_Prediction = data.frame(Barcode=rownames(wholePipelineTesting), Predicted=pred_spat)
  
  spat_score = spatialValidation(TRUTH,rf2_Prediction)
  
  #spat_score[[1]]: correct_pred_normal/(nrow(Truth_normal)+epsilon)
  #spat_score[[2]]: incorrectly_classified_truthNormal/(nrow(Truth_normal)+epsilon)
  #spat_score[[3]]: correct_pred_tumor/(nrow(Truth_tumor)+epsilon)
  #spat_score[[4]]: incorrectly_classified_truthTumor/(nrow(Truth_tumor)+epsilon)
  
  ans<-list()
  ans[[1]]<-rf2_y
  ans[[2]]<-zz1[[4]][11] ##balanced accuracy
  ans[[3]]<-spat_score[[1]] ##correct_pred_normal %
  ans[[4]]<-spat_score[[2]] ##incorrectly_classified_truthNormal %
  ans[[5]]<-spat_score[[3]] ##correct_pred_tumor %
  ans[[6]]<-spat_score[[4]] ##incorrectly_classified_truthTumor %
  
  w1=w2=1/2
  
  ##### maximize the whole thing. reducthe incorrect stuff. so negative sign in front of it. maximize the correct stuff. Also TCGABA is 
  ## as important as the spatial validation. So w1=w2=1/2
  ans[[7]] = w1*(ans[[2]]) + w2*((spat_score[[1]]) - (spat_score[[2]]) + (spat_score[[3]]) - (spat_score[[4]]))
  
  ans
}

ver = 2.2

date = Sys.Date()

qq<-read.csv(file = paste0("THCA_ge_standardized_ver",ver,".csv"), header = T)
rownames(qq)<-qq$X
count_stand<-qq[,-1]
tt<-read.csv(file = paste0("THCA_ge_standardized_annotation_ver",ver,".csv"), header = T)
ww2<-tt


##sanity check
all(colnames(count_stand)==ww2$ID)

table(ww2$Type)

interest = "Braf.like"

suppl1<-read.csv(file = paste0("THCA_OURClass_2.csv"), header = T)

table(suppl1$OurClass)

ww2$Type2<-ww2$Type
indx<-c()
id2<-gsub("-",".",suppl1$sample)
found<-0
for(i in 1:length(id2)){
  t1<-grep(id2[i],ww2$ID)
  if(length(t1)>0){
    ww2$Type2[t1[1]]<-suppl1$OurClass[i]
    found<-found+1
  }
}
found

table(ww2$Type2)

ww2$Type2<-gsub("-",".",ww2$Type2)
table(ww2$Type2)

nrmlIndx = which(ww2$Type2=="Normal")
indxOfInterest = which(ww2$Type2==interest)

##sanity check
all(colnames(count_stand)==ww2$ID)

## take only normal and tumor-of-interest
count_stand2<-count_stand[,c(nrmlIndx,indxOfInterest)]
ww3<-ww2[c(nrmlIndx,indxOfInterest),]

##sanity check
all(colnames(count_stand2)==ww3$ID)

table(ww3$Type2)

if(interest == "Braf.like") idate = "2026-04-10"
if(interest == "Ras.like") idate = "2026-04-10"
if(interest == "RET.fusion") idate = "2026-04-10"

pSigVal = 0.01

basicStats<-read.csv(file = paste0(interest,"-vs-Normal-basicStats_",idate,".csv"), header = T)
deg2<-basicStats$X[which(basicStats$pvalue_Adj<pSigVal)]
tGenes = deg2 ## unused genes 



st_name<-c("M1","M2","P1","P2","T1","T2","N1","N2")
index = 1
sotti = read.csv(file = paste0(st_name[index],"_TruthForRFGA.csv"), header = T)

##change the barcode to match
bc = gsub("-",".",sotti$X)
TRUTH=data.frame(X=bc,value=sotti$value)

allSpatData<-read.csv(file = paste0(st_name[index],"_rawData.csv"), header = T)
rownames(allSpatData)<-allSpatData$X
count_stand0<-allSpatData[,-1]

######### this piece of code pads the spat data to make it the same size as deg2 data #########

deg3<-intersect(deg2,rownames(count_stand0))

if(length(deg2)==length(deg3)){
  ##
  print("All Same. Nothing to do")
  count_stand1 = count_stand0
  ##
}else {
  #
  foundIndx<-which(deg2 %in% allSpatData$X)
  notFound<-deg2[-foundIndx]
  
  toBeAdded<-matrix(0, nrow = length(notFound), ncol = ncol(count_stand0))
  rownames(toBeAdded)<-notFound
  colnames(toBeAdded)<-colnames(count_stand0)
  count_stand1 = rbind(count_stand0,toBeAdded)
  #
}

deg3<-intersect(deg2,rownames(count_stand1))
##sanity check
length(deg2)==length(deg3)

deg3Indx<-which(rownames(count_stand1) %in% deg3)
GE2<-as.matrix(count_stand1[deg3Indx,])
dim(GE2)
GE3<-convert_0_1matrix(GE2,'row')
dim(t(GE3))
wholePipelineTesting<-as.data.frame(t(GE3))

######### this piece of code pads the spat data to make it the same size as TCGA data #########

noImprovementPerIteration = 0

### initiate population ###

tcga_scores = read.csv(file = paste0(interest,"_vs_Normal_DifferentThresholds_",idate,".csv"), header = T)
tcga_scores$NetScore = 0

for(i in 1:nrow(tcga_scores)){
  
  tcga_scores$NetScore[i] = 0.5*tcga_scores$TCGABalancedAccuracy[i] + 0.5*(tcga_scores$correctly_classified_truthNormal[i] - tcga_scores$incorrectly_classified_truthNormal[i] + tcga_scores$correctly_classified_truthTumor[i] - tcga_scores$incorrectly_classified_truthTumor[i])
  
}



if(interest == "Braf.like") cutoff = 0.9
if(interest == "Ras.like") cutoff = 0.5
if(interest == "RET.fusion") cutoff = 0.5

tcga_scores = tcga_scores[which(tcga_scores$NetScore>cutoff),]

population<-list()
#listOfChromosomes<-list()  ##to be used for looking up scores and hopefully reduce time 

fitnessScores<-c() ## store the scores of the chromosomes in population
#listOfChromosomeScores<-list()  ## same index as listOfChromosomes. to be used for looking up scores and hopefully reduce time 


bestFitness=-1 ##$NetScore
bestSize=0
bestSolnBA=0 ##$TCGABalancedAccuracy
bestSolnccTN=0 ##$correctly_classified_truthNormal
bestSolnicTN=0 ##$incorrectly_classified_truthNormal
bestSolnccTT=0 ##$correctly_classified_truthTumor
bestSolnicTT=0 ##$incorrectly_classified_truthTumor

bestIndx = -1

for(i in 1:nrow(tcga_scores)){
  
  knownGood<-read.csv(file = paste0("deGenes_",interest,"_vs_Normal_thrshld",tcga_scores$TumorThrshld[i],"_",idate,".csv"), header = T)
  temp<-knownGood$deGenes
  
  
  population[[i]] = temp
  #listOfChromosomes[[i]] = temp
  print(paste0("Created population[[",i,"]]"))
  
  fitnessScores[i] = tcga_scores$NetScore[i]
  
  if(tcga_scores$NetScore[i]>bestFitness)  {
    bestIndx = i
    bestFitness = tcga_scores$NetScore[i]
    bestSize = tcga_scores$deGenes[i]
    bestSolnBA=tcga_scores$TCGABalancedAccuracy[i]##$TCGABalancedAccuracy
    bestSolnccTN=tcga_scores$correctly_classified_truthNormal[i] ##$correctly_classified_truthNormal
    bestSolnicTN=tcga_scores$incorrectly_classified_truthNormal[i] ##$incorrectly_classified_truthNormal
    bestSolnccTT=tcga_scores$correctly_classified_truthTumor[i] ##$correctly_classified_truthTumor
    bestSolnicTT=tcga_scores$incorrectly_classified_truthTumor[i] ##$incorrectly_classified_truthTumor
  }
  
}

popSize = length(population)

bestSoln<-population[[bestIndx]]
bestClassifier<-readRDS(file = paste0(interest,"_vs_Normal_RF2_thrshld",tcga_scores$TumorThrshld[bestIndx],"_",idate,".rds"))

for(i in 1:popSize){
  print(paste0("Length of population[[",i,"]]: ",length(population[[i]])))
  print(paste0("Score of population[[",i,"]]: ",fitnessScores[i]))
}

print(paste0("bestFitness: ",bestFitness))
print(paste0("bestSolution Size: ",length(bestSoln)))

  print(paste0("BestFitness: ",bestFitness))
  print(paste0("BestSoln_BA: ",bestSolnBA))
  print(paste0("BestSoln_ccTN: ",bestSolnccTN))
  print(paste0("BestSoln_icTN: ",bestSolnicTN))
  print(paste0("BestSoln_ccTT: ",bestSolnccTT))
  print(paste0("BestSoln_icTT: ",bestSolnicTT)) 
### to make life easier. To reduce time 
discardedChrom<-list()
discardedIndx=1

allGenesUsedIndx<-c()
for(i in 1:popSize){
  rr = which(tGenes %in% population[[i]])
  print(length(rr))
  allGenesUsedIndx<-c(allGenesUsedIndx,rr)
}

allGenesUsedIndx2<-unique(allGenesUsedIndx)
yetToBeUsedGenes<-tGenes[-allGenesUsedIndx2]

print(paste0("yetToBeUsedGenes: ",length(yetToBeUsedGenes)))

iteration<-1
terminationCondition = 1001 #noImprovementInBest for 1000 iterations

print(paste0("Fitness Scores: "))
print(fitnessScores)


noImprovementInBest=0
improvementInBest=F

populationUnchanged = 0
alpha=4 ## maximum number iterations, all chroms in population are same before we forcefully introduce mutation


theta = 0.001 ## minimum difference in fitness scores
# For a 5 in 100 chance, the probability is 0.05.
low_probability <- 0.1
fiftyfifty = 0.5 ##determines single or double point crossover


# Configure parallel processing
# Detect the number of available cores and leave 2 free for the OS
num_cores <- 20
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
cat(paste("Registered parallel backend with", num_cores, "cores\n")) # Optional status message



set.seed(721)
while(noImprovementInBest<terminationCondition){
  
  ## resetting values
  improvementInBest = F ##setting for new iteration
  
  print(paste0("Inside While Loop Iteration: ",iteration))
  
  if(noImprovementInBest>=terminationCondition){
    break
  }
  
  print(paste0("Fitness Scores: "))
  print(fitnessScores)
  
  
  relFitnessScores = (fitnessScores/sum(fitnessScores))
  partialSums=relFitnessScores
  for(i in 2:length(fitnessScores)){
    partialSums[i]<-partialSums[i-1]+relFitnessScores[i]
  }
  
  print(paste0("Partial Sums: "))
  print(partialSums)
  
  set.seed(iteration)
  
  S<-runif(1,min=0,max=1)
  pIndx1<-rouletteSelection(partialSums,S)
  parent1<-population[[pIndx1]]
  S<-runif(1,min=0,max=1)
  pIndx2<-rouletteSelection(partialSums,S)
  while(pIndx2==pIndx1){
    print(paste0("inside little while"))
    S<-runif(1,min=0,max=1)
    pIndx2<-rouletteSelection(partialSums,S)
  }
  parent2<-population[[pIndx2]]
  
  print(paste0("Parent indices: ",pIndx1," & ",pIndx2))
  
  
  if(runif(1,0,1)>fiftyfifty){
    print(paste0("Single Point Crossover"))
    offspring<-singlePntcrossOver(parent1,parent2)
  } else {
    print(paste0("Double Point Crossover"))
    offspring<-doublePntcrossOver(parent1,parent2)
  }
  
  print(paste0("Size_parent1: ",length(parent1)," Size_parent2: ",length(parent2)))
  print(paste0("Size_offspring[[1]]: ",length(offspring[[1]])," Size_offspring[[2]]: ",length(offspring[[2]])))
  
  ##mutation step . make low probability event
  
  mut1<-runif(1,0,1)
  print(paste0("mut1: ",mut1))
  if(populationUnchanged>alpha){
    mut1 = low_probability-1
    print(paste0("Need to introduce mutation. So ","mut1: ",mut1))
  }
  print(paste0("yetToBeUsedGenes: ",length(yetToBeUsedGenes)))
  if(mut1 < low_probability && length(yetToBeUsedGenes)>0)  {
    #whichOffspring<-sample(c(1:2),1)
    print(paste0("Offspring[1] mutates"))
    answers = mutate(yetToBeUsedGenes,offspring[[1]])
    offspring[[1]] = answers[[1]]
    yetToBeUsedGenes = answers[[2]]
  }
  
  
  
  mut2<-runif(1,0,1)
  print(paste0("yetToBeUsedGenes: ",length(yetToBeUsedGenes)))
  if(populationUnchanged>alpha){
    mut2 = low_probability-1
    print(paste0("Need to introduce mutation. So ","mut2: ",mut2))
  }
  print(paste0("mut2: ",mut2))
  if(mut2 < low_probability && length(yetToBeUsedGenes)>0)  {
    print(paste0("Offspring[2] mutates"))
    answers = mutate(yetToBeUsedGenes,offspring[[2]])
    offspring[[2]] = answers[[1]]
    yetToBeUsedGenes = answers[[2]]
  }
  print(paste0("yetToBeUsedGenes: ",length(yetToBeUsedGenes)))
  
  #### check for duplication of genes. if yes, delete the duplicates
  for(j in 1:length(offspring)){
    
    dupYes = areThereDuplicates(offspring[[j]])
    if(dupYes){
      dupIndx = which(duplicated(offspring[[j]]))
      offspring[[j]] = offspring[[j]][-dupIndx]
      print(paste0("duplicate gene deleted!!"))
    }
    
    print(paste0("Size of offspring[",j,"]: ",length(offspring[[j]])))
    
  } ##each offspring
  
  ## score the offsprings
  ## score the offsprings
  offspringFitness<-rep(0,length(offspring))
  discardOffspring<-rep(0,length(offspring))  ## vector to keep track of offspring viability ## 0: okay, 1: to be discarded, 2: already discarded
  alreadySeen<-rep(0,length(offspring))  ## vector to keep track whether generation already present in population ## 0: okay, 1: already in population
  
  print(paste0("# of discarded chromosome: ", length(discardedChrom)))
  
  ###check the discarded pile and stored chromosome list
  for(j in 1:length(offspring)){
    print(paste0("Size of Offspring[[",j,"]]:",length(offspring[[j]])))
    discard = isAlreadyPresent(discardedChrom,offspring[[j]])
    print(paste0("Offspring[[",j,"]] already in discarded pile? ",discard))
    if(discard>-1) {
      discardOffspring[j] = 2 ## already discarded
    }
    seen = isAlreadyPresent(population,offspring[[j]])
    print(paste0("Offspring[[",j,"]] already seen in population? ",seen))
    if(seen>-1){
      offspringFitness[j] = fitnessScores[seen]
      alreadySeen[j] = 1
    }
    
  }
  
  for(i in 1:length(offspring)){
    if(discardOffspring[i]==0 && alreadySeen[i]==0){
      
      tmp<-getChromosomeFitnessScore(offspring[[i]])
      
      # ans[[1]]<-rf2_y
      # ans[[2]]<-zz1[[4]][11] ##balanced accuracy
      # ans[[3]]<-spat_score[[1]] ##correctly_classified_truthNormal %
      # ans[[4]]<-spat_score[[2]] ##incorrectly_classified_truthNormal %
      # ans[[5]]<-spat_score[[3]] ##correctly_classified_truthTumor %
      # ans[[6]]<-spat_score[[4]] ##incorrectly_classified_truthTumor %
      # ans[[7]] = w1*(ans[[2]]) + w2*((spat_score[[1]]) - (spat_score[[2]]) + (spat_score[[3]]) - (spat_score[[4]]))
      
      
      offspringFitness[i]=tmp[[7]]
      
      if((offspringFitness[i]-bestFitness)>theta){
        bestFitness = offspringFitness[i]
        bestSoln = offspring[[i]]
        bestClassifier = tmp[[1]]
        bestSolnBA=tmp[[2]] ##$TCGABalancedAccuracy
        bestSolnccTN=tmp[[3]] ##$correctly_classified_truthNormal
        bestSolnicTN=tmp[[4]] ##$incorrectly_classified_truthNormal
        bestSolnccTT=tmp[[5]] ##$correctly_classified_truthTumor
        bestSolnicTT=tmp[[6]] ##$incorrectly_classified_truthTumor
        improvementInBest=T
        noImprovementInBest=0
      }   
    }
  }
  
  
  print(paste0("Offspring Fitness Scores: "))
  print(offspringFitness)
  
  
  ##code for eliminating Weaker solutions
  
  
  for(i in 1:length(offspringFitness)){
    if(discardOffspring[i]==0){
      ss<-sort.int(fitnessScores, decreasing = F, index.return = T)
      if(offspringFitness[i]<=fitnessScores[ss$ix[1]]){ ## if offspring fitness lower than the lowest pop fitness score, discard it
        
        discardedChrom[[discardedIndx]] = offspring[[i]]
        discardedIndx = discardedIndx + 1
        print(paste0("Offspring[",i,"] added to discard pile"))
        
      } else {
        
        ## discard the lowest population chromosome
        discardedChrom[[discardedIndx]] = population[[ss$ix[1]]]
        discardedIndx = discardedIndx + 1
        print(paste0("population[[",ss$ix[1],"]] added to discard pile"))
        population[[ss$ix[1]]]<-offspring[[i]]
        fitnessScores[ss$ix[1]]<-offspringFitness[i]
      }
    } else if(discardOffspring[i]==1){ 
      discardedChrom[[discardedIndx]] = offspring[[i]]
      discardedIndx = discardedIndx + 1
      print(paste0("Offspring[",i,"] added to discard pile"))
    } else {
      ##already discarded
      print(paste0("Offspring[",i,"] already present in discard pile"))
    }
  }
  
  print(paste0("Fitness Scores: "))
  print(fitnessScores)
  
  print(paste0("bestFitness: ",bestFitness))
  
  if(!improvementInBest){
    noImprovementInBest=noImprovementInBest+1
  }
  
  fitnessScores1 = round(fitnessScores, digits = 5)
  
  ## if populationUnchanged>alpha and if in this iteration the fitness score of the chroms in population are different
  
  if(populationUnchanged>alpha && length(unique(fitnessScores1))>1){
    populationUnchanged = 0 ##reset the counter to just under limit
  }
  
  ## if fitness score of all the chroms in population are same
  
  if(length(unique(fitnessScores1))==1){
    print(paste0("Population unchanged for ",populationUnchanged, " iterations!!"))
    populationUnchanged = populationUnchanged + 1
  }
  
  print(paste0("noImprovementInBest: ",noImprovementInBest))
  print(paste0("End of iteration: ",iteration))
  print(paste0("BestFitness: ",bestFitness))
  print(paste0("BestSoln_BA: ",bestSolnBA))
  print(paste0("BestSoln_ccTN: ",bestSolnccTN))
  print(paste0("BestSoln_icTN: ",bestSolnicTN))
  print(paste0("BestSoln_ccTT: ",bestSolnccTT))
  print(paste0("BestSoln_icTT: ",bestSolnicTT)) 
  iteration = iteration + 1
  
  
  
}## termination
# Stop the parallel cluster when done
stopCluster(cl)
registerDoSEQ() # Register sequential backend again
cat(paste("Registered sequential backend again\n")) # Optional status message

#print(paste0("fitnessScores[ss$ix[1]]: ",fitnessScores[ss$ix[1]]))
print(paste0("BestFitness: ",bestFitness))
print(paste0("BestSoln_BA: ",bestSolnBA))
print(paste0("BestSoln_ccTN: ",bestSolnccTN))
print(paste0("BestSoln_icTN: ",bestSolnicTN))
print(paste0("BestSoln_ccTT: ",bestSolnccTT))
print(paste0("BestSoln_icTT: ",bestSolnicTT))
print(paste0("SolutionSize: ",length(bestSoln)))
df11<-data.frame(Fitness=bestFitness,TCGABA=bestSolnBA,ccTN=bestSolnccTN,icTN=bestSolnicTN,ccTT=bestSolnccTT,icTT=bestSolnicTT)
write.csv(df11, file = paste0("Normal_",interest,"_GA_Performance_SolnSize=",length(bestSoln),"_",date,".csv"), row.names = F)
df12<-data.frame(deGenes=bestSoln)
write.csv(df12, file = paste0("Normal_",interest,"_GA_deGenes_SolnSize=",length(bestSoln),"_",date,".csv"), row.names = F)
saveRDS(bestClassifier, file = paste0("Normal_",interest,"_GA_Classifier_SolnSize=",length(bestSoln),"_",date,".rds"))







