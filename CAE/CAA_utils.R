
#function to compute L2 norm
L2norm <- function(x){sqrt(sum(x^2))}
L1norm <- function(x){sum(abs(x))}


#distance between two CAA projections
#Calculates distance between a pair of CAA spaces
#inputs: vector c1 = (u1,v1) c2=(u2,v2) -> concatenated u, v
#       scalar k: length of u to know where to break c
CAA_dist<-function(c1,c2,k){
  dif1<-c1-c2
  dif2<-c1-c(c2[(k+1):(2*k)],c2[1:k])
  d<-min(L2norm(dif1[1:k])+L2norm(dif1[(k+1):(2*k)]),L2norm(dif2[1:k])+L2norm(dif2[(k+1):(2*k)]))
  return(d)
}

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}


#Distance to compute distance matrix between all pairs in a set of points
dist_matrix<-function(All_U,All_V){
  Dist<-matrix(ncol=ncol(All_U),nrow=ncol(All_U))
  diag(Dist)<-0
  ###Build distance matrix
  for(i in 1:(ncol(All_U)-1)){
    for(j in (i+1):ncol(All_U)){
      dist<-CAA_dist(append(c(All_U[,i]),c(All_V[,i])),append(c(All_U[,j]),c(All_V[,j])),nrow(All_U))
      Dist[i,j]<-dist
      Dist[j,i]<-dist
    }
  }
  return(Dist)
}

#Distance to compute distance matrix across two sets of points
dist_matrix_cross<-function(All_U_1,All_V_1,All_U_2,All_V_2){
  Dist<-matrix(nrow=ncol(All_U_1),ncol=ncol(All_U_2))
  ###Build distance matrix
  for(i in 1:ncol(All_U_1)){
    for(j in 1:ncol(All_U_2)){
      dist<-CAA_dist(append(c(All_U_1[,i]),c(All_V_1[,i])),append(c(All_U_2[,j]),c(All_V_2[,j])),nrow(All_U_1))
      Dist[i,j]<-dist
    }
  }
  return(Dist)
}


CAA_ray<-function(u,v){
  ray<-(All_U[,i]+All_V[,i])/length(u)
  return(ray)
}

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

angle_matrix<-function(All_rays){
  Dist<-matrix(ncol=ncol(All_rays),nrow=ncol(All_rays))
  diag(Dist)<-0
  ###Build distance matrix
  for(i in 1:(ncol(All_rays)-1)){
    for(j in (i+1):ncol(All_rays)){
      dist<-angle(All_rays[,i],All_rays[,j])
      Dist[i,j]<-dist
      Dist[j,i]<-dist
    }
  }
  return(Dist)
}
  
#Distance between two objects (i.e. patients) based on their CAA characterization
#Input: CAA_obj_1, CAA_obj_2: CAA characterization of two objects/patients, the type of object is the same as output of CAA main function
CAA_object_dist<-function(CAA_char_1,CAA_char_2,metric='max'){
  Dist<-dist_matrix_cross(CAA_char_1$All_U, CAA_char_1$All_V, CAA_char_2$All_U, CAA_char_2$All_V)
  Dist_min_1<-apply(Dist,1,min)
  Dist_min_2<-apply(Dist,2,min)
  if(metric=='max'){
    return(max(c(Dist_min_1,Dist_min_2)))
  }
  if(metric=='mean'){
    return(mean(c(Dist_min_1,Dist_min_2)))
  }
  #calculate 
}

#Distance between two CAA objects doing stepwise pairwise assignment of distances
CAA_pairwise_dist<-function(CAA_char_1,CAA_char_2,metric='mean'){
  Dist<-dist_matrix_cross(CAA_char_1$All_U, CAA_char_1$All_V, CAA_char_2$All_U, CAA_char_2$All_V)
  n_pairs<-min(ncol(CAA_char_1$All_U),ncol(CAA_char_2$All_U))
  pairwise_dist<-rep(NA,n_pairs)
  for( i in 1:(n_pairs-1)){
    pairwise_dist[i]<-min(Dist)
    idx_min<-which(Dist==min(Dist),arr.ind = T)
    Dist<-Dist[-idx_min[1],-idx_min[2]]
  }
  pairwise_dist[n_pairs]<-min(Dist)
  if(metric=='max'){
    return(max(pairwise_dist))
  }
  if(metric=='mean'){
    return(mean(pairwise_dist))
  }
  #calculate 
}

#transform labels so that 'yes'->1 and 'no'->-1
label_to_num<-function(label){
  if(label=='Yes'|label=='a'){
    return(1)
  }
  else{
    return(0)
  }
}


plot_MDS<-function(M_dist,color_label,alph,siz){
  fit_mds <- cmdscale(M_dist,eig=TRUE, k=2) 
  plot_mds <- as.data.frame(fit_mds$points)
  colnames(plot_mds) <-c('x_coord','y_coord')
  plot_mds['label']<-color_label
  qplot(x_coord, y_coord, data=plot_mds, colour = label, alpha=I(alph), size = I(siz), main = 'MDS projection')
}

