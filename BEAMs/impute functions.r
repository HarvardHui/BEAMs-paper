#### Uncovering the Consequences of BEAMs in Omics Data Analysis ####
# Author: Harvard Wai Hann Hui
# Date: 09/February/2023
# Script: Functions

##################################################################
####### Create correlation matrix (Shah et al. 2017) 
##################################################################
CorrMatrixNegFixed <- function(blocks, a, corr, off){
  Corr <- NULL
  for (i in 1:blocks) {
    Corr[[i]] <- matrix(NA, ncol = a, nrow = a)
    diag= rep(1,a)
    offdiag = corr
    Corr[[i]][lower.tri(Corr[[i]])] <- offdiag
    Corr[[i]][upper.tri(Corr[[i]])] <- t(Corr[[i]])[upper.tri(t(Corr[[i]]))]
    diag(Corr[[i]]) <- diag
    for ( k in ((ncol(Corr[[i]])/2)+1):(ncol(Corr[[i]])) ) {
      for( j in 1:(nrow(Corr[[i]])/2) ) {
        Corr[[i]][j,k] <- -1*Corr[[i]][j,k]
      }
    }
    
    for ( k in ((nrow(Corr[[i]])/2)+1):(nrow(Corr[[i]])) ) {
      for( j in 1:(ncol(Corr[[i]])/2) ) {
        Corr[[i]][k,j] <- -1*Corr[[i]][k,j]
      }
    }
    
  }
  res <- as.matrix(bdiag(Corr))
  res[which(res == 0, arr.ind = TRUE)] <- off
  return(res)
}


##################################################################
####### Create dataframe from corr. matrix (Shah et al. 2017)
# n = no. of samples
# p = no. of features
# covar = correlation matrix
##################################################################
SimulatedData <- function(n, p, covar, low, high){
  Means <- runif(p, min = low, max = high)
  data <- rmvnorm(n, mean = Means, sigma = covar)
  return(data)
}


##################################################################
####### Class effect simulation (Affects the first 60 rows)
##################################################################
class_effect <- function(df, class_factor, halfway=TRUE){
  if(halfway==TRUE)
  {
    df[1:60,c(which(class_factor==0))]=df[1:60,c(which(class_factor==0))]*1.5
  }else{
    df[,c(which(class_factor==0))]=df[,c(which(class_factor==0))]*1.5
  }
  return(df)
}


##################################################################
####### Batch effect simulation (Affects the first 10 columns)
# Z = multiplicative factor
# Y = additive factor
##################################################################
batch_effect<-function(df, Z, Y, multiplicative=FALSE, additive=FALSE){
  if(multiplicative==TRUE && additive==FALSE){
    batch_df <- cbind(df[,1:10]*Z, df[,11:20])#create batch effect(first 10 cols vs last 10 cols)
  } else if (multiplicative==FALSE && additive==TRUE){
    batch_df <- cbind(df[,1:10]+Y, df[,11:20])
  } else if (multiplicative==TRUE && additive==TRUE){
    batch_df <- cbind(Z*(df[,1:10]+Y), df[,11:20])#z(x+Y), where y is the additive factor, and z is the multiplicative factor.
  } else{
    batch_df<-df
  }
  return(batch_df)
  
}


##################################################################
####### Simulation of Missing values
# total = total MV proportion
# mcar = MCAR MVs proportion
##################################################################
mv.sim<-function(df,total,mcar){
  df2<-2^df
  subdfhole=c()
  iterr=0
  for(number in 1:ncol(df2)){
    #here generating MNAR
    mis_prop=total-mcar
    data_res = df2[,number]
    a=data_res
    if(mis_prop!=0){
      cutoff = quantile(a, mis_prop,na.rm=T)
      a[a< cutoff] = NA
      data_res=a
    }
    #here generating MCAR
    mis_prop=mcar
    q=round(length(data_res)*mis_prop)
    tmplis=which(!is.na(data_res))
    mi = sample(tmplis, q)
    data_res[mi] = NA
    
    subdfhole=cbind(subdfhole,data_res)
  }
  colnames(subdfhole)=colnames(df2)
  jid=c()
  
  ### Which rows are completely NA
  for(j in 1:nrow(subdfhole)){
    tmp=subdfhole[j,]
    ge=tmp[is.na(tmp)]
    if(length(ge)==ncol(subdfhole))
      jid=append(jid,j)
  }
  
  rownames(subdfhole)<-rownames(df2)
  
  df2<-log(df2,base=2)
  subdfhole<-log(subdfhole,base=2)
  protein<-rownames(subdfhole)
  return(list(df2,subdfhole))
}


##################################################################
####### M1 imputation by KNN
##################################################################
impute.M1 <- function (df, k, mv.threshold = 0.8){
  rownames(df) <- seq(nrow(df))
  drop_inds=c()
  for (i in 1:nrow(df)){
    if ((length(which(is.na(df[i,])))/ncol(df))>=mv.threshold){
      drop_inds<-c(drop_inds,i)
    }
  }
  
  df2<-df
  if (length(drop_inds)>0){
    df2<-df[-drop_inds,]
  }
  
  # impute by KNN
  imp_df <- impute.knn(df2, k=k, colmax=0.99, rowmax=0.99)
  out_df <- imp_df$data
  
  rownames(out_df) <- rownames(df2)
  return(out_df)
}


##################################################################
####### M2 imputation by KNN
##################################################################
impute.M2 <- function (df, k, batch_factor, mv.threshold = 0.8){
  rownames(df) <- seq(nrow(df))
  batches <- unique(batch_factor)
  batch_list <- list()
  for (h in 1:length(batches)){
    batch_list[[h]]<- df[,which(batch_factor==h)]
  }

  drop_inds=c()
  for (p in 1:length(batches)){
    cur_batch <- batch_list[[p]]
    for (i in 1:nrow(df)){
      if ((length(which(is.na(cur_batch[i,])))/ncol(cur_batch))>=mv.threshold){
        drop_inds<-c(drop_inds,i)
      }
    }
  }
  drop_inds <- unique(drop_inds)
  df2 <- df
  if (length(drop_inds)>0){
    df2 <- df[-drop_inds,]
  }
  
  for (g in 1:length(batches)){
    batch_list[[g]]<- df2[,which(batch_factor==g)]
  }
  # impute by KNN
  
  out_df<-c()
  for (j in 1:length(batches)){
    cur_imp_batch <- impute.knn(batch_list[[j]], k=k, colmax=0.99, rowmax=0.99)
    cur_imp_batch <- cur_imp_batch$data
    out_df <- cbind(out_df, cur_imp_batch)
  }

  rownames(out_df) <- rownames(df2)
  return(out_df)
}


##################################################################
####### Batch effect correction by ComBat
# pdat should include class and batch factors
##################################################################
do.combat.sim <- function(df, pdat){
  df = as.data.frame(df)
  batch <- pdat$batch
  class <- pdat$class
  
  mod = model.matrix(~as.factor(class), data=df)
  
  combat_edata = ComBat(dat=df, batch=batch, mod=mod, par.prior = T, prior.plots = F)
  
  return(combat_edata)
}


##################################################################
####### Perform t-test based on class factor
# output either "features" for significant features or "pvals"
# for p-values of all features
##################################################################
t.test.sim <- function(df, class_factor, output="features"){
  class1 = df[,which(class_factor==0)]
  class2 = df[,which(class_factor==1)]
    result=c()
    for(i in 1:nrow(df))
    {
      if(sd(class1[i,])<1e-6 && sd(class2[i,])<1e-6)
      {
        class1[i,1]=jitter(class1[i,1])
      }
      a=t.test(class1[i,],class2[i,])
      result=rbind(result,a$p.value)
    }
    result=p.adjust(result,method = "BH")
    result=as.data.frame(result)
    rownames(result)=rownames(class1)
    
    sig.feats = which(result<0.05)

    if (output == "features"){
      return(sig.feats)
    }
    if (output == "pvals"){
      return(result)
    }
}


##################################################################
####### Reduce dataframe to same features
# if both == FALSE, only df1 is reduced to df2's feature space
# if both == TRUE, both df1 and df2 are reduced to same space
##################################################################
same.features <- function(df1, df2, both = FALSE){ # df1 = data to be reduced, df2 = reference
  df1.features = rownames(df1)
  df2.features = rownames(df2)
  
  same.features <- Reduce(intersect, list(df1.features,df2.features))
  feature_ind1 <- match(same.features, df1.features)
  out_df1 <- df1[feature_ind1,]
  
  feature_ind2 <- match(same.features, df2.features)
  out_df2 <- df2[feature_ind2,]
  if (both == TRUE){
    return(list(out_df1,out_df2))
  }
  if (both == FALSE){
    return(out_df1)
  }
}


##################################################################
####### Get TPR/FPR values of imputed data
# true and false cases are determined using the supplied true dataset
##################################################################
tprfpr <- function(true, imp, class_factor, output="tpr"){
  tps <- t.test.sim(true, class_factor, output="features")
  imp.true.features <- tps[which(tps %in% rownames(imp))]
  
  imp.ps <- t.test.sim(imp, class_factor, output="features")
  imp.tps <- imp.ps[which(imp.ps %in% imp.true.features)]
  tpr <- (length(imp.tps)/(length(imp.tps) + (length(imp.true.features)-length(imp.tps)))) * 100
  
  imp.fps <- imp.ps[-which(imp.ps %in% imp.true.features)]
  imp.negatives <- rownames(imp)[-which(rownames(imp) %in% imp.ps)]
  imp.fn <- which(imp.negatives %in% tps)
  if (length(imp.fn) > 1){
    imp.tn <- imp.negatives[-imp.fn]
  } else{
    imp.tn <- imp.negatives
  }
  fpr <- (length(imp.fps)/(length(imp.fps) + length(imp.tn))) * 100
  if (output == "tpr"){
    return(tpr)
  }
  if (output == "fpr"){
    return(fpr)
  }
}


##################################################################
####### Check if results are significantly different
##################################################################
check.p <- function(res_df){
  output=c()
  iter=seq(from=1,to=ncol(res_df), by=2)
  for (i in iter){
    v = i+1
    t = t.test(res_df[,i], res_df[,v])
    output = c(output, t$p.val)
  }
  return(output)
}
