#CAA with disjoint support penalty

#input:
#X: matrix containing the data
#c1, c2: L1 constraints for u and v, respectively

library(TunePareto)

#function to compute L2 norm
L2norm <- function(x){sqrt(sum(x^2))}
L1norm <- function(x){sum(abs(x))}



obtain_r2<-function(u,v,x){
  y<-x%*%v
  f<-x%*%u
  y_hat<-mean(y)
  SS_tot<-sum((y-y_hat)^2)
  SS_res<-sum((y-f)^2)
  return(1-SS_res/SS_tot)
}


Generate_partition<-function(n,labels,supervised=TRUE){
    foldList<- generateCVRuns(labels=labels,
                              ntimes = 1,
                              nfold = n,
                              leaveOneOut = FALSE,
                              stratified = supervised)
    partition<-foldList[[1]]  #needs to be modified if doing CV more times
  }


#Softthresholding function. w = u or v
Soft_thres<-function(Pw,lambda){
  w_output<-sign(Pw)*pmax(0, abs(Pw)-lambda)
  #if(all(w_output==0)){return(NULL)}
  w_output<-w_output/L2norm(w_output)
  return(w_output)
}



###function to solve convex problem at each iteration
#input: w=either u or v from previous iteration (whichever is fixed), X, c (c1 or c2, penalty that constrains vector being optimized)
update_w<-function(Co, S, ci, w){
  #initialize output
  w_output <- rep(0,times=ncol(Co))
  Pw <- Co%*%w
  abs_w<- abs(w)
  w_Si<-rep(0,length(abs_w))
  for(r in 1:length(abs_w)){
    w_Si[r]<-sum(abs_w[which(S[r,]!=0)],na.rm=T)
  }
  lambda1<- max(abs(Pw[(abs_w!=0)]/abs_w[(abs_w!=0)]))
  lambda2<- 0
  w_output<- Soft_thres(Pw,lambda2+lambda1*w_Si)
  if(all(is.na(w_output))){
    return(NULL)}
  if(L1norm(w_output) <= ci){
    return(w_output)
    
  }else{
    lambda2_min <- 0
    lambda2_max <- max(Pw[(w_Si==0)])
    while(lambda2_max-lambda2_min > 1e-5){
    lambda2<-(lambda2_min+lambda2_max)/2
    w_output<-Soft_thres(Pw,lambda2+lambda1*w_Si)
    if(L1norm(w_output) < ci || sum(w_output!=0)<=1){
      lambda2_max<-lambda2
      }
    else{
      lambda2_min<-lambda2
    }
  }
  return(w_output)
}
}



###Function to obtain one CAA factor
CAA_it <- function(Co,S,u_init,v_init,c1,c2){
  converged=FALSE
  u_prev <- u_init
  v_prev <- v_init
  numIter = 0
  while(converged==FALSE && numIter < 10000){
    u_update <- update_w(Co,S,c1,v_prev)
    v_update <- update_w(Co,S,c2,u_update)
    if(any(c(is.null(u_update),is.null(v_update)))){
      return(NULL)}
    if(L2norm(u_update-u_prev)<1e-5&&L2norm(v_update-v_prev)<1e-5){
      converged = TRUE
    }
    u_prev <- u_update
    v_prev <- v_update
    numIter = numIter + 1
  }
  if(converged==FALSE){print('Did not converge')}
  return(list(u=u_update,v=v_update))
}


###main CAA function

CAA<-function(X, penalty1, penalty2, S, kproj, scale=TRUE,double_init=TRUE){
  start.time<-Sys.time()
  if(scale==T){
    X<-scale(X,T,T)
  }
  X[is.na(X)]<-mean(X[!is.na(X)])
  Co<-t(X)%*%X
  diag(Co)<-0
  ###Initialize u and v by finding pair of most correlated features
  # (this is the max off-diagonal entry in the covariance matrix)
  N <- ncol(Co)
  U=matrix(,nrow=N,ncol=kproj)
  V=matrix(,nrow=N,ncol=kproj)
  r2coef <-rep(NA,kproj)
  d_rep <-rep(NA,kproj)
  i=1
  r2 = 1
  while(i<=kproj && r2>0){
  ind_max_cor <- which(abs(Co)==max(abs(Co)),arr.ind=T)[1,]

  u_init <- rep(0,times=N)
  v_init <- rep(0,times=N)
  u_init[ind_max_cor[1]]<-1 #c1
  v_init[ind_max_cor[2]]<-1 #c2

  if(double_init==TRUE){
  c1<-0.5*sqrt(N)
  c2<-0.5*sqrt(N)
  CAA_output<-CAA_it(Co,S,u_init,v_init,c1,c2)
  if(!is.null(CAA_output)){
    u_init<-CAA_output$u
    v_init<-CAA_output$v
  }else{
    i=kproj+1
    next}
  }
  c1<-penalty1*sqrt(N)
  c2<-penalty2*sqrt(N)
  CAA_output<-CAA_it(Co,S,u_init,v_init,c1,c2)
  if(!is.null(CAA_output)){
  u_new<-CAA_output$u
  v_new<-CAA_output$v
  d<-as.numeric(t(u_new)%*%Co%*%v_new)
  U[,i]<-u_new
  V[,i]<-v_new
  r2<-obtain_r2(u_new,v_new,X)
  r2coef[i]<-r2
  d_rep[i]<-d
  Co <- Co - d*(u_new%*%t(v_new)+v_new%*%t(u_new))
  i=i+1
  }else{
    i=kproj+1
  }
  }
keep<-(r2coef > 0)&!is.na(r2coef)
if(all(keep==FALSE)){
  print('Error in first convergence')
  return(NULL)
}
U<-as.matrix(U[,keep])
V<-as.matrix(V[,keep])
#Numbering<-c(1:kproj)
D<-t(as.matrix(cbind(r2coef,d_rep)))
D<-as.matrix(D[,keep])
#D<-D[,ordered]
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
return(list(U=U,V=V,D=D))
}



CAA_tune<-function(X,scale=TRUE,steps=10,c_min='predet',c_max='predet',folds=5){   ####INCLUDE WARM STARTS
  start.time<-Sys.time()
  if(scale==T){
    X<-scale(X,T,T)
  }
  N <- ncol(X)
  r <- nrow(X)
  #Co<-t(X)%*%X
  ###Initialize u and v by finding pair of most correlated features
  # (this is the max off-diagonal entry in the covariance matrix)

  partition<-Generate_partition(folds,rep(0,r),FALSE)
  d_rep <-matrix(NA,nrow=folds,ncol=steps)
  if(c_min=='predet'){
    c_min=1/sqrt(N)
  }
  if(c_max=='predet'){
    c_max=1
  }
  tune_grid = seq(c_min,c_max,(c_max-c_min)/steps)
  for (p in 1:folds){
    partition_p<-unlist(partition[p])
    X_p<-X[partition_p,]
    Co<-t(X_p)%*%X_p
    Co_nodiag <- Co
    diag(Co_nodiag) <- rep(0,times=N)
    ind_max_cor <- which(Co_nodiag==min(Co_nodiag),arr.ind=T)[1,]

    u_init <- rep(0,times=N)
    v_init <- rep(0,times=N)
    u_init[ind_max_cor[1]]<-1 #c1
    v_init[ind_max_cor[2]]<-1 #c2

    for(i in rev(1:steps)){
      d_folds<-rep(NA,folds)
      c_param<-tune_grid[i]*sqrt(N)
      CAA_output<-tryCatch(CAA_it(Co,u_init,v_init,c_param,c_param))

    if(!is.null(CAA_output)){
      u_init<-CAA_output$u
      v_init<-CAA_output$v
      d_rep[p,i]<-as.numeric(t(u_init)%*%Co%*%v_init)
    }
    i=i+1
  }

  p=p+1

  }
  d_steps<-colMeans(d_rep)
  idx_d_max<-which.max(d_steps)
  return(c_opt=tune_grid[idx_d_max])
}

CAA_steptune<-function(X,kproj,scale=TRUE,steps=10,c_min='predet',c_max='predet',folds=5){
  start.time<-Sys.time()
  if(scale==T){
    X<-scale(X,T,T)
  }
  #Co_k<-t(X)%*%X
  N <- ncol(X)
  r <- nrow(X)
  #Co<-t(X)%*%X
  ###Initialize u and v by finding pair of most correlated features
  # (this is the max off-diagonal entry in the covariance matrix)
  #U=matrix(,nrow=N,ncol=kproj)
  #V=matrix(,nrow=N,ncol=kproj)
  #r2coef <-rep(NA,kproj)
  #d_final <-rep(NA,kproj)
  partition<-Generate_partition(folds,rep(0,r),FALSE)
  d_rep <-matrix(NA,nrow=folds,ncol=steps)
  if(c_min=='predet'){
    c_min=1/sqrt(N)
  }
  if(c_max=='predet'){
    c_max=1
  }
  tune_grid = seq(c_min,c_max,(c_max-c_min)/steps)
  for (k in 1:kproj){

    for (p in 1:folds){
      partition_p<-unlist(partition[p])
      X_p<-X[partition_p,]
      Co<-t(X_p)%*%X_p
      Co_nodiag <- Co
      diag(Co_nodiag) <- rep(0,times=N)
      ind_max_cor <- which(Co_nodiag==min(Co_nodiag),arr.ind=T)[1,]

      u_init <- rep(0,times=N)
      v_init <- rep(0,times=N)
      u_init[ind_max_cor[1]]<-1 #c1
      v_init[ind_max_cor[2]]<-1 #c2
      for(i in rev(1:steps)){
        d_folds<-rep(NA,folds)
        c_param<-tune_grid[i]*sqrt(N)
        CAA_output<-tryCatch(CAA_it(Co,u_init,v_init,c_param,c_param))

        if(!is.null(CAA_output)){
          u_init<-CAA_output$u
          v_init<-CAA_output$v
          d_rep[p,i]<-as.numeric(t(u_init)%*%Co%*%v_init)
        }
        i=i+1
      }

      p=p+1

    }
    d_steps<-colMeans(d_rep)
    c_opt=tune_grid[which.max(d_steps)]
    c1<-c_opt*sqrt(ncol(X))
    c2<-c_opt*sqrt(ncol(X))
    Co_k_nodiag <- Co_k
    diag(Co_k_nodiag) <- rep(0,times=N)
    ind_max_cor <- which(Co_k_nodiag==min(Co_k_nodiag),arr.ind=T)[1,]

    u_0 <- rep(0,times=N)
    v_0 <- rep(0,times=N)
    u_0[ind_max_cor[1]]<-1 #c1
    v_0[ind_max_cor[2]]<-1 #c2
    CAA_output<-tryCatch(CAA_it(Co,u_0,v_0,c1,c2))
    if(!is.null(CAA_output)){
      u_new<-CAA_output$u
      v_new<-CAA_output$v
      d<-as.numeric(t(u_new)%*%Co%*%v_new)
      U[,k]<-u_new
      V[,k]<-v_new
      r2<-obtain_r2(u_new,v_new,X)
      r2coef[k]<-r2
      d_final[k]<-d
      Co_k <- Co_k - d*(u_new%*%t(v_new)+v_new%*%t(u_new))
      k=k+1
    }
    else{
      k=kproj+1
    }
  }
    keep<-(colSums(is.na(U)) == 0)
    U<-as.matrix(U[,keep])
    V<-as.matrix(V[,keep])
    Numbering<-c(1:kproj)
    D<-t(as.matrix(cbind(r2coef,Numbering,d_final)))
    #D<-D[,ordered]
    D<-as.matrix(D[,keep])
    keep_posR2<-which(D[1,]>0)
    D<-D[,keep_posR2]
    U<-U[,keep_posR2]
    V<-V[,keep_posR2]
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    return(list(U=U,V=V,D=D))
}
#
# CAA_output<-CAA_steptune(X,5)
# # c_opt$c_opt<-CAA_tune(X)
# # c_opt<-c_opt$c_opt
# # CAA_output<-CAA(X,c_opt,c_opt,5,F,T)
# # # # #
# U<-CAA_output$U
# # # #
# V<-CAA_output$V
# D<-CAA_output$D
