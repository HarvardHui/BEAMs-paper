# Data simulation functions (Shah et al. 2017)
CorrMatrixNegFixed <- function(blocks, a, corr, off) {
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

SimulatedData <- function(n, p, covar, low, high)  {
  ## No of Samples == N
  ## No of Features = p
  Means <- runif(p, min = low, max = high)
  data <- rmvnorm(n, mean = Means, sigma = covar)
  return(data)
}

##### Class effect simulation for the first 60 rows
class_effect <- function(df, class_factor, halfway=TRUE)
{
  if(halfway==TRUE)
  {
    df[1:60,c(which(class_factor==0))]=df[1:60,c(which(class_factor==0))]*1.5
  }else{
    df[,c(which(class_factor==0))]=df[,c(which(class_factor==0))]*1.5
  }
  return(df)
}

##### Batch effect simulation (Affects the first 10 columns)
batch_effect<-function(df, Z, Y, multiplicative=FALSE, additive=FALSE)
{
  if(multiplicative==TRUE && additive==FALSE){
    batch_df <- cbind(df[,1:10]*Z, df[,11:20])#create batch effect(first 10 cols vs last 10 cols)
  } else if (multiplicative==FALSE && additive==TRUE){
    batch_df <- cbind(df[,1:10]+Y, df[,11:20])#Maybe we can calculate global average, take the square root and then add
  } else if (multiplicative==TRUE && additive==TRUE){
    batch_df <- cbind(Z*(df[,1:10]+Y), df[,11:20])#z(x+Y), where y is the additive factor, and z is the multiplicative factor.
  } else{
    batch_df<-df
  }
  return(batch_df)
  
}

# Simulate MV function
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

impute.M1 <- function (df, mv.threshold = 0.8){
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
  imp_df <- impute.knn(df2, k=3, colmax=0.99, rowmax=0.99)
  out_df <- imp_df$data
  
  rownames(out_df) <- rownames(df2)
  return(out_df)
}

impute.M2 <- function (df, batch_factor, mv.threshold = 0.8){
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
    cur_imp_batch <- impute.knn(batch_list[[j]], k=3, colmax=0.99, rowmax=0.99)
    cur_imp_batch <- cur_imp_batch$data
    out_df <- cbind(out_df, cur_imp_batch)
  }

  rownames(out_df) <- rownames(df2)
  return(out_df)
}
  
do.combat.sim <- function(df, pdat){
  df = as.data.frame(df)
  batch <- pdat$batch
  class <- pdat$class
  
  mod = model.matrix(~as.factor(class), data=df)
  
  combat_edata = ComBat(dat=df, batch=batch, mod=mod, par.prior = T, prior.plots = F)
  
  return(combat_edata)
}
  
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
} # [[1]] if only df1 is reduced, [[2]] also if df2 is reduced

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


# KNN optimization
do.knn<-function(df, k){
  imp.knn <- df
  imp.knn[is.finite(df) == FALSE] <- NA
  mv_ind<-which(is.na(df), arr.ind = TRUE)
  mv_rows<-unique(mv_ind[,1])
  
  klist<-list()
  
  for (i in mv_rows){
    klist[[i]]<-list()
    feature_mvs<-which(is.na(df[i,]))
    feature_obs<-which(!is.na(df[i,]))
    cand_vectors<-df[,feature_obs]
    for (j in 1:length(feature_mvs)){
      feature_mvs_ind<-feature_mvs[j]
      target_vector<-df[,feature_mvs_ind]
      
      if (length(feature_obs) < k){
        imp.knn[i,feature_mvs_ind]<-mean(df[,feature_mvs_ind], na.rm = TRUE)
      } else {
        calc.dist<-data.frame((cand_vectors-target_vector)^2)
        dist<-sqrt(colMeans(calc.dist, na.rm = TRUE))
        dist[is.nan(dist) | is.na(dist)] <- Inf
        dist[dist==0] <- ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
        
        if (sum(is.finite(dist)) < k) {
          #stop(message = "Fewer than K finite distances found")
          imp.knn[i,feature_mvs_ind]<-mean(df[,feature_mvs_ind], na.rm = TRUE)
        } else {
          k_sample_ind <- order(dist)[1:k]
          k_samples <- feature_obs[k_sample_ind]
          wghts <- 1/dist[k_sample_ind]/sum(1/dist[k_sample_ind])
          imp.knn[i, feature_mvs_ind] <- wghts %*% df[i, k_samples]
        }
      }
      klist[[i]][[feature_mvs_ind]]<-k_samples
    }
  }
  return(list(original_data=df, imputed_data=imp.knn, klist=klist))
}