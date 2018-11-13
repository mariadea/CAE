rm(list=ls())

library(reshape)
library(plyr)
library(rpart)
library(caret)

source("../CAE/CAA_utils.R")
source("../CAE/CAA_classify.R")

#load labels
patient_summary<-read.csv('../data/PatientSummary.csv')

setwd("../data/EEG/")

filenames<-list.files()

#subset of patients who survived
patient_summary_lived <- patient_summary[patient_summary$SURV=='Lived',]
lived <-as.vector(patient_summary_lived$SubjectID)

#only take into account filenames of patients who survive
filenames_lived<-rep(NA,length(lived))

i=0
for (f in filenames){
  nam <- substr(f,5,(nchar(f)-4))
  if (nam%in%lived){
    filenames_lived[i]<-f
    i = i+1
  }
}

###Matrix specifying sets for disjoint support - this prevents trivial correlations from ocurring
S = matrix(0,nrow=(66*3),ncol=(66*3))
diag(S)<-1
nrows=66
for ( i in c(24:28,(24+nrows):(28+nrows),(24+2*nrows):(28+2*nrows))) {
  for (j in c(24:28,(24+nrows):(28+nrows),(24+2*nrows):(28+2*nrows))){
    S[i,j]<-1}
}
for ( i in c(29:33,(29+nrows):(33+nrows),(29+2*nrows):(33+2*nrows))) {
  for (j in c(29:33,(29+nrows):(33+nrows),(29+2*nrows):(33+2*nrows))){
    S[i,j]<-1}
}
for ( i in c(34:38,(34+nrows):(38+nrows),(34+2*nrows):(38+2*nrows))) {
  for (j in c(34:38,(34+nrows):(38+nrows),(34+2*nrows):(38+2*nrows))){
    S[i,j]<-1}
}
for ( i in c(39:43,(39+nrows):(43+nrows),(39+2*nrows):(43+2*nrows))) {
  for (j in c(39:43,(39+nrows):(43+nrows),(39+2*nrows):(43+2*nrows))){
    S[i,j]<-1}
}

#choose time window parameters
#how many seconds before the end of recording to start/stop

t_start <- 122400
t_final <- 129600
min_entries <- 3600


#EXTRACT PORTION OF EEG TO BE USED IN PREDICTION

filenames = filenames_lived
p_proj<-c()
label_proj<-c()
All_U<-matrix(nrow=66*3,ncol=0)
All_V<-matrix(nrow=66*3,ncol=0)
All_D<-matrix(nrow=2,ncol=0)
EEG = list()
i = 1
for (f in filenames){
  print(f)
  EEG_i = read.csv(f)
  #transform time to proper format
  EEG_i$ClockDateTime <- as.POSIXlt(EEG_i$ClockDateTime)
  EEG_i <- EEG_i[order(EEG_i$ClockDateTime),]
  #identify window to train CAA on
  EEG_i<- EEG_i[which(t_start < difftime(EEG_i$ClockDateTime,EEG_i[1,4],units="secs") &
                          difftime(EEG_i$ClockDateTime,EEG_i[1,4],units="secs") < t_final),]
  
  if(nrow(EEG_i) < min_entries){next}
  #remove columns that are not neccesary
  EEG_i = EEG_i[,6:71]
  #prevent errors that arise from a column of only zeros
  all_zeros <- c()
  for (c in 1:66){
    if(all( EEG_i[, c]==0)){
      all_zeros<-c(all_zeros,c)
    }
  }
  EEG_i[,all_zeros]<-NULL
    
  n_feat<-ncol(EEG_i)
    
  #matrix containing sets to enforce disjoint support over
  S_i<-S
  if(!is.null(all_zeros)){
    S_i<-S_i[-c(all_zeros,all_zeros+66,all_zeros+(2*66)),-c(all_zeros,all_zeros+66,all_zeros+(2*66))]
  }

#APPLY CAA
      
  EEG_i<-cbind(EEG_i,EEG_i^2,EEG_i^3)
    
  CAA_output<-CAA(EEG_i,0.2,0.2,S_i,ncol(EEG_i))#S)
  if(is.null(CAA_output)){next}
  #store CAA projections found
  U<-CAA_output$U
  V<-CAA_output$V
  if(length(all_zeros)>0){
    for(c in all_zeros){
      U <- rbind(as.matrix(U[1:c,]),rep(0,ncol(U)),as.matrix(U[-(1:c),]))
      V <- rbind(as.matrix(V[1:c,]),rep(0,ncol(V)),as.matrix(V[-(1:c),]))
      EEG_i<-cbind(as.matrix(EEG_i[,1:c]),rep(0,nrow(EEG_i)),as.matrix(EEG_i[,-(1:c)]))
      }
      for(c in all_zeros){
        U <- rbind(as.matrix(U[1:(66+c),]),rep(0,ncol(U)),as.matrix(U[-(1:(66+c)),]))
        V <- rbind(as.matrix(V[1:(66+c),]),rep(0,ncol(V)),as.matrix(V[-(1:(66+c)),]))
        EEG_i<-cbind(as.matrix(EEG_i[,1:c]),rep(0,nrow(EEG_i)),as.matrix(EEG_i[,-(1:c)]))
      }
      for(c in all_zeros){
        U <- rbind(as.matrix(U[1:(2*66+c),]),rep(0,ncol(U)),as.matrix(U[-(1:(2*66+c)),]))
        V <- rbind(as.matrix(V[1:(2*66+c),]),rep(0,ncol(V)),as.matrix(V[-(1:(2*66+c)),]))
        EEG_i<-cbind(as.matrix(EEG_i[,1:c]),rep(0,nrow(EEG_i)),as.matrix(EEG_i[,-(1:c)]))
      }
    }
    
  n_proj<-ncol(U)
  nam <- substr(f,5,(nchar(f)-4))
  p_proj<-append(p_proj,rep(nam,n_proj))
  label_proj<-append(label_proj,rep(as.vector(patient_summary[patient_summary$SubjectID==nam , 'Good_Outcome']),n_proj))
  All_U<-cbind(All_U,U)
  All_V<-cbind(All_V,V)
  D<-CAA_output$D
  All_D<-cbind(All_D,D)
    EEG[[i]] <-EEG_i
    i = i+1
  }
  proj_data<-cbind(p_proj,label_proj)


save(EEG,All_U,All_V,All_D, proj_data,file='../CAA_EEG.RData')

