#Algorithms for classification based on CAA distance


#k-nn algorithm
#Input:  Dist: distance matrix
#         M_proj_info: matrix with columns 'SUBJECT_ID', 'WEIGHT', 'LABEL' 
#                   ('LABEL' should be 1 or -1, if instead it is 'Yes' and 'No' it is internally changed)
#         idx_train, idx_test: boolean vectors indicating if elements in train/test sets
#         k: parameter specifying number of neighbors to consider

CAA_knn<-function(Dist, proj_info, ids_train, ids_test, k, thres){
  eps<-10^-8
  patients_test<-as.vector(unique(proj_info[ids_test,'SUBJECT_ID']))
  pred <- matrix(NA,ncol=2,nrow=length(patients_test))
  i=1
  if(!all(unique(proj_info[,'LABEL'])%in%c(1,0))){
    proj_info[,'LABEL']<-unlist(lapply(proj_info[,'LABEL'],FUN=label_to_num))
  }
  proj_info_train<-proj_info[ids_train,]
  
  for (p in patients_test){
    #identify which projections correspond to the patient
    idx_p<-which(proj_info[,'SUBJECT_ID']==p)
    Dist_p = Dist[idx_p,ids_train]
    if(is.matrix(Dist_p)==FALSE){
      Dist_p<-t(as.matrix(Dist_p))
    }
    pred_p = rep(NA,nrow(Dist_p))
    weights_p<-as.numeric(proj_info[idx_p,'WEIGHT'])
    for(r in 1:nrow(Dist_p)){
      #reorder distances in ascending order and choose k closest
      ordDist_p_r<-order(Dist_p[r,])
      k_nearest_dist<-Dist_p[r,ordDist_p_r][1:k]
      k_nearest<-ordDist_p_r[which(!is.na(Dist_p[r,ordDist_p_r][1:k]))]
      #mean
      pred_r<-mean(as.numeric(proj_info_train[k_nearest,'LABEL']))
      if(!is.na(pred_r)){
      if(abs(pred_r-0.5)>thres){pred_p[r]<-pred_r}
      }
      #weighted mean of labels of k neighbors, where weights correspond to 1-(r2_i-r2_j)
      #pred_p[r]<-sum(as.numeric(proj_info_train[k_nearest,'LABEL'])*(1-(weights_p[r]-as.numeric(proj_info_train[k_nearest,'WEIGHT']))))/sum(1-(weights_p[r]-as.numeric(proj_info_train[k_nearest,'WEIGHT'])))
    }
    #pred[i,]<-c(p,mean(pred_p,na.rm=T))
    pred_p<-na.omit(pred_p)
    
    pred[i,]<-c(p,log(prod((pred_p+eps)/(1-pred_p+eps))))
    
    #pred[i,]<-c(p,mean(pred_p,na.rm=T))
    i = i+1
  }
  return(pred)
}



knn_spar_tune<-function(Dist, proj_info, k_options,spar_options,thres_options, f=10){
  patient_info<-unique(proj_info[,c('SUBJECT_ID','LABEL')])
  if(!all(unique(proj_info[,'LABEL'])%in%c(1,0))){
    proj_info[,'LABEL']<-unlist(lapply(proj_info[,'LABEL'],FUN=label_to_num))
  }
  set.seed(1)
  Folds_tune<-createFolds(patient_info[,'LABEL'],k=f)
  #create table with all combinations of grid search
  param_grid<-matrix(NA,ncol=4,nrow=length(k_options)*length(spar_options)*length(thres_options))
  #<-rev(rep(spar_options,each=length(k_options)))
  #param_grid[,2]<-rep(k_options,times=length(spar_options))
  j=1
  for(spar_opt in rev(spar_options)){
    Dist[Dist > spar_opt] <-NA
  for(thres_opt in thres_options){
    for(k_option in k_options){
    #table with results (true and pred labels)
    Outcomes<-matrix(,ncol=2,nrow=0)
    for(fold_t in Folds_tune){
      fold_t = unlist(fold_t)
      patients_fold_t<-as.vector(patient_info[fold_t,1])
      ids_train<-!(proj_info[,'SUBJECT_ID']%in%patients_fold_t)
      ids_test<-proj_info[,'SUBJECT_ID']%in%patient_info[fold_t,1]
      pred_fold<-CAA_knn(Dist, proj_info, ids_train, ids_test, k_option,thres_opt)
      #add to table of results - true and pred values-
      Outcomes<-rbind(Outcomes,pred_fold)
    }
    colnames(Outcomes)<-c('SUBJECT_ID','score')
    Outcomes<-unique(merge(Outcomes,proj_info[,1:2],by='SUBJECT_ID'))
    Outcomes[,'score']<-as.numeric(as.character(Outcomes[,'score']))
    #calculate roc and record
    roc_test <- roc(Outcomes[,3],Outcomes[,2])
    param_grid[j,1]<-spar_opt
    param_grid[j,2]<-thres_opt
    param_grid[j,3]<-k_option
    param_grid[j,4] = auc(roc_test)
    j = j+1
  }
  }
  }
  #choose best performance using roc
  #return optimal k, and performance at that level
  return(param_grid[which.max(param_grid[,4]),])
  
}


#k-fold cross validation to tune k-nn
# Input: Dist: Distance matrix
#        proj_info: matrix with columns 'SUBJECT_ID', 'WEIGHT', 'LABEL'
#        k_grid: grid over which to find optimal k (through cross-validation)
#        f: number of folds in cross-validation to tune k
k_nn_tune<-function(Dist, proj_info, k_grid, f=10){
  #5-fold cv partition
  patient_info<-unique(proj_info[,c('SUBJECT_ID','LABEL')])
  set.seed(1)
  Folds_tune<-createFolds(patient_info[,'LABEL'],k=f)
  #for each k: estimate performance of k-nn using 5-fold cv
  auc_grid = rep(NA,length(k_grid))
  j=1
  if(!all(unique(proj_info[,'LABEL'])%in%c(1,-1))){
    proj_info[,'LABEL']<-unlist(lapply(proj_info[,'LABEL'],FUN=label_to_num))
  }
  for(k_option in k_grid){
    #table with results (true and pred labels)
    Outcomes<-matrix(,ncol=2,nrow=0)
    for(fold_t in Folds_tune){
      fold_t = unlist(fold_t)
      patients_fold_t<-as.vector(patient_info[fold_t,1])
      ids_train<-!(proj_info[,'SUBJECT_ID']%in%patients_fold_t)
      ids_test<-proj_info[,'SUBJECT_ID']%in%patient_info[fold_t,1]
      pred_fold<-CAA_knn(Dist, proj_info, ids_train, ids_test, k_option)
      #add to table of results - true and pred values-
      Outcomes<-rbind(Outcomes,pred_fold)
    }
    colnames(Outcomes)<-c('SUBJECT_ID','score')
    Outcomes<-unique(merge(Outcomes,proj_info[,1:2],by='SUBJECT_ID'))
    Outcomes[,'score']<-as.numeric(as.character(Outcomes[,'score']))
    #calculate roc and record
    roc_test <- roc(Outcomes[,3],Outcomes[,2])
    auc_grid[j] = auc(roc_test)
    j = j+1
  }
  #choose best performance using roc
  #return optimal k, and performance at that level
  return(list(k = k_grid[which.max(auc_grid)],auc = max(auc_grid)))
  
}


#####Cluster vote
#Input: Dist: Distance matrix
#       proj_info: matrix with columns 'SUBJECT_ID', 'LABEL'
#         idx_train, idx_test: boolean vectors indicating if elements in train/test sets
#         k: parameter specifying number of neighbors to consider

CAA_ClusterVote<-function(Dist,proj_info,ids_train,ids_test,k){
  Dist_train<-Dist[ids_train,ids_train]
  Clus<-pam(Dist_train,k=k,diss=TRUE)
  clus_class<-data.frame(Clus$clustering,proj_info[ids_train,'LABEL'])
  colnames(clus_class)<-c("cluster","label_proj")
  clus_div<-table(clus_class$label_proj,clus_class$cluster)
  #order clusters according to ratio of classes
  ratio<-c()
  if(all(unique(proj_info[,'LABEL'])%in%c('Yes','No'))){
  for(x in colnames(clus_div)){
    ratio<-append(ratio,clus_div[1,x]/(clus_div['Yes',x]+clus_div['No',x]))
  }
  }else{
  if(all(unique(proj_info[,'LABEL'])%in%c('a','c'))){
    for(x in colnames(clus_div)){
      ratio<-append(ratio,clus_div['a',x]/(clus_div['a',x]+clus_div['c',x]))
    }
  }else{print('I dont recognize the labels')}}
  patients_test<-unique(proj_info[ids_test,1])
  pred_label_test<-matrix(NA,ncol=2,nrow=length(patients_test))
  i=1
  for(p in patients_test){
    p_id<-which(proj_info[,1]==p)
    Dist_p_med<-Dist[p_id,Clus$medoids]
    if(length(p_id)==1){
      Dist_p_med<-matrix(Dist_p_med,nrow=1)
    }
    p_clusters<-unique(apply(Dist_p_med,1,FUN=which.min))
    pred_label_test[i,]<-c(p,mean(ratio[p_clusters]))
    i=i+1
  }
  return(pred_label_test)
}

clustervote_tune_k<-function(Dist, proj_info, k_grid, f=10){
  patient_info<-unique(proj_info[,c('SUBJECT_ID','LABEL')])
  set.seed(1)
  Folds_tune<-createFolds(patient_info[,'LABEL'],k=f)
  #for each k: estimate performance of k-nn using 5-fold cv
  auc_grid = rep(NA,length(k_grid))
  j=1
  for(k_option in k_grid){
    #table with results (true and pred labels)
    Outcomes<-matrix(,ncol=2,nrow=0)
    for(fold_t in Folds_tune){
      fold_t = unlist(fold_t)
      ids_train<-!(proj_info[,'SUBJECT_ID']%in%patient_info[fold_t,'SUBJECT_ID'])
      ids_test<-proj_info[,'SUBJECT_ID']%in%patient_info[fold_t,'SUBJECT_ID']
      pred_fold<-CAA_ClusterVote(Dist, proj_info,ids_train,ids_test,k_option)
      #add to table of results - true and pred values-
      Outcomes<-rbind(Outcomes,pred_fold)
    }
    colnames(Outcomes)<-c('SUBJECT_ID','score')
    Outcomes<-unique(merge(Outcomes,proj_info[,c('SUBJECT_ID','LABEL')],by='SUBJECT_ID'))
    Outcomes[,'score']<-as.numeric(as.character(Outcomes[,'score']))
    #calculate roc and record
    roc_test <- roc(Outcomes[,3],Outcomes[,2])
    auc_grid[j] = auc(roc_test)
    j = j+1
  }
  #choose best performance using roc
  #return optimal k, and performance at that level
  return(list(k = k_grid[which.max(auc_grid)],auc = max(auc_grid)))
}

