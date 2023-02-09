#### Uncovering the Consequences of BEAMs in Omics Data Analysis ####
# Author: Harvard Wai Hann Hui
# Date: 09/February/2023
# Script: Simulation of BEAMs

############################################################################################################################
#### This Script creates datasets with control and BEAMs MV proportions.
#### It then performs KNN MVI by M1 and M2 strategies, and analyzes the
#### imputation accuracy, TPR, and FPR.
#### 
#### Please check the "Set Parameters" section to ensure that parameters are
#### set correctly, or as desired
############################################################################################################################

###### Load all necessary libraries and functions
library(here)
library(corrplot)
library(gtools)
library(ggfortify)
library(ggplot2)
library(mice)
library(sva)
library(vegan)
library(imputeLCMD)
library(limma)
library(reshape2)
library(missForest)

wd <- here()
source(paste0(wd, "/BEAMs/impute functions.r"))

##### Set parameters
set.seed(2023)
class_factor <- c(rep(c(rep(0,5), rep(1,5)),2))
batch_factor <- c(rep(1,10),rep(2,10))
Z = 1.2 # Batch multiplication factor
Y = sqrt(5) # Batch additive factor
N = 100 # N number of iterations
k = 3 # k nearest neighbours
mv.threshold=0.8 # % of MVs accepted per feature

# Use these for severe BEAMs
#mv_prop1=c(0.1,0.025)
#mv_prop2=c(0.5,0.125)

# Use these for moderate BEAMs
mv_prop1=c(0.2,0.05)
mv_prop2=c(0.4,0.1)

# Control MV level
mv_prop3=c(0.3,0.075)

###### Save results
rmse_res <- data.frame(matrix(ncol=0,nrow=0))
feature_res <- data.frame(matrix(ncol=0,nrow=0))
tpr_res <- data.frame(matrix(ncol=0,nrow=0))
fpr_res <- data.frame(matrix(ncol=0,nrow=0))

############## Data simulation (Functions from Shah et al. 2017)
# Create correlation matrix
CorrelationMatrix0.7 <- CorrMatrixNegFixed(30, 30, 0.7, 0.2)
# Create dataframe with 20 samples and 900 features
DataSimul <- SimulatedData(20, 900, as.matrix(nearPD(CorrelationMatrix0.7)[[1]]), low = 5, high =15)
DataSimul <- t(DataSimul)

# Create phenotype data to store class and batch factors
pdat = data.frame(class=class_factor, batch=batch_factor)

############# Class effect simulation
DataSimul_class <- class_effect(DataSimul, class_factor, halfway=TRUE)

############# Batch effect simulation
batch_df <- batch_effect(DataSimul_class, Z, Y, multiplicative = TRUE, additive = TRUE)
rownames(batch_df)<-seq(nrow(batch_df))

# Split dataframe by batch for BEAMs simulation
batch.1 <- batch_df[,which(batch_factor==1)]
batch.2 <- batch_df[,which(batch_factor==2)]

# PCA to check batch effect
pc<-t(batch_df)
pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
pca_res<-prcomp(pc, scale. = FALSE)
pca_df <-data.frame(pca_res$x)
pca_df$batch <- as.factor(batch_factor)
pca_df$class <- as.factor(class_factor)
ggplot(pca_df, aes(x = PC1, y = PC2, color = batch, shape = class)) +
  geom_point(size=3) +
  xlab("PC1") +
  ylab("PC2") + theme_bw()

# PCA after batch effect correction by ComBat
pc.bc<-do.combat.sim(as.data.frame(batch_df),pdat)
pc<-t(pc.bc)
pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
pca_res<-prcomp(pc, scale. = FALSE)
pca_df <-data.frame(pca_res$x)
pca_df$batch <- as.factor(batch_factor)
pca_df$class <- as.factor(class_factor)
ggplot(pca_df, aes(x = PC1, y = PC2, color = batch, shape = class)) +
  geom_point(size=3) +
  xlab("PC1") +
  ylab("PC2") + theme_bw()

############################################ Start of simulation iterations #########################

for (i in 1:N){
  set.seed(i)
  # Simulate BEAMs
  b1.beam.df <- mv.sim(batch.1, mv_prop1[1], mv_prop1[2])
  b2.beam.df <- mv.sim(batch.2, mv_prop2[1], mv_prop2[2])
  
  beam.df <- cbind(b1.beam.df[[2]], b2.beam.df[[2]])
  rownames(beam.df)<-seq(nrow(beam.df))
  
  # Simulate Control MVs
  control.df <- mv.sim(batch_df, mv_prop3[1], mv_prop3[2])
  control.df <- control.df[[2]]
  rownames(control.df)<-seq(nrow(control.df))
  
  # M1 (Global) imputation
  beam.m1 <- impute.M1(beam.df, k=k, mv.threshold=mv.threshold)
  control.m1 <- impute.M1(control.df, k=k, mv.threshold=mv.threshold)
  
  # M2 (Batch-sensitized) imputation
  beam.m2 <- impute.M2(beam.df, k=k, batch_factor, mv.threshold=mv.threshold)
  control.m2 <- impute.M2(control.df, k=k, batch_factor, mv.threshold=mv.threshold)
  
  # Apply ComBat batch effect correction
  beam.m1.bc <- do.combat.sim(beam.m1, pdat) # imp
  control.m1.bc <- do.combat.sim(control.m1, pdat) # imp
  
  beam.m2.bc <- do.combat.sim(beam.m2, pdat) # imp
  control.m2.bc <- do.combat.sim(control.m2, pdat) # imp
  
  batch_df.bc <- do.combat.sim(batch_df, pdat) # true
  
  
  ##### get same feature space for true and missing datasets
  mis.beam.m1 <- same.features(beam.df, beam.m1.bc) # mis
  mis.control.m1 <- same.features(control.df, control.m1.bc) # mis
  true.beam.m1 <- same.features(batch_df.bc, beam.m1.bc) # true
  true.control.m1 <- same.features(batch_df.bc, control.m1.bc) #true
  
  mis.beam.m2 <- same.features(beam.df, beam.m2.bc) # mis
  mis.control.m2 <- same.features(control.df, control.m2.bc) # mis
  true.beam.m2 <- same.features(batch_df.bc, beam.m2.bc) # true
  true.control.m2 <- same.features(batch_df.bc, control.m2.bc) # true
  
  ##### same feature space within each method
  beam.m1.bc.s <- same.features(beam.m1.bc, control.m1.bc, both=TRUE)
  control.m1.bc.s <- beam.m1.bc.s[[2]] #imp
  beam.m1.bc.s <- beam.m1.bc.s[[1]] #imp
  true.m1.s <- same.features(batch_df.bc, beam.m1.bc.s) # true
  mis.beam.m1.s <- same.features(mis.beam.m1, beam.m1.bc.s) # mis
  mis.control.m1.s <- same.features(mis.control.m1, beam.m1.bc.s) #mis
  
  control.m2.bc.s <- same.features(control.m2.bc, beam.m2.bc, both=TRUE)
  beam.m2.bc.s <- control.m2.bc.s[[2]] #imp
  control.m2.bc.s <- control.m2.bc.s[[1]] #imp
  true.m2.s <- same.features(batch_df.bc, beam.m2.bc.s) # true
  mis.beam.m2.s <- same.features(mis.beam.m2, beam.m2.bc.s) # mis
  mis.control.m2.s <- same.features(mis.control.m2, beam.m2.bc.s) #mis
  
  # get NRMSE
  rmse.beam.m1 <- nrmse(beam.m1.bc, mis.beam.m1, true.beam.m1)
  rmse.control.m1 <- nrmse(control.m1.bc, mis.control.m1, true.control.m1)
  
  rmse.beam.m2 <- nrmse(beam.m2.bc, mis.beam.m2, true.beam.m2)
  rmse.control.m2 <- nrmse(control.m2.bc, mis.control.m2, true.control.m2)
  
  rmse.beam.m1.s <- nrmse(beam.m1.bc.s, mis.beam.m1.s, true.m1.s)
  rmse.control.m1.s <- nrmse(control.m1.bc.s, mis.control.m1.s, true.m1.s)
  
  rmse.beam.m2.s <- nrmse(beam.m2.bc.s, mis.beam.m2.s, true.m2.s)
  rmse.control.m2.s <- nrmse(control.m2.bc.s, mis.control.m2.s, true.m2.s)
  
  rmse_vals <- c(rmse.beam.m1, rmse.control.m1, rmse.beam.m2, rmse.control.m2, 
                 rmse.beam.m1.s, rmse.control.m1.s, rmse.beam.m2.s, rmse.control.m2.s)
  rmse_res <- rbind(rmse_res, rmse_vals)
  
  
  # Get number of features retained in each imputation
  f.b.m1<-nrow(beam.m1)
  f.c.m1<-nrow(control.m1)
  
  f.b.m2<-nrow(beam.m2)
  f.c.m2<-nrow(control.m2)
  
  same.m1 <- nrow(control.m1.bc.s)
  same.m2 <- nrow(control.m2.bc.s)
  
  feature_vals <- c(f.b.m1, f.c.m1, f.b.m2, f.c.m2, same.m1, same.m2)
  feature_res <- rbind(feature_res, feature_vals)
  
  
  # get TPR
  tpr.b.m1 <- tprfpr(batch_df.bc, beam.m1.bc, class_factor, output="tpr")
  tpr.c.m1 <- tprfpr(batch_df.bc, control.m1.bc, class_factor, output="tpr")
  
  tpr.b.m2 <- tprfpr(batch_df.bc, beam.m2.bc, class_factor, output="tpr")
  tpr.c.m2 <- tprfpr(batch_df.bc, control.m2.bc, class_factor, output="tpr")
  
  tpr.b.m1.s <-tprfpr(true.m1.s, beam.m1.bc.s, class_factor, output="tpr")
  tpr.c.m1.s <-tprfpr(true.m1.s, control.m1.bc.s, class_factor, output="tpr")
  
  tpr.b.m2.s <-tprfpr(true.m2.s, beam.m2.bc.s, class_factor, output="tpr")
  tpr.c.m2.s <-tprfpr(true.m2.s, control.m2.bc.s, class_factor, output="tpr")
  
  tpr_vals <- c(tpr.b.m1, tpr.c.m1, tpr.b.m2, tpr.c.m2, tpr.b.m1.s, tpr.c.m1.s,
                tpr.b.m2.s, tpr.c.m2.s)
  tpr_res <- rbind(tpr_res, tpr_vals)
  
  
  # get FPR
  fpr.b.m1 <- tprfpr(batch_df.bc, beam.m1.bc, class_factor, output="fpr")
  fpr.c.m1 <- tprfpr(batch_df.bc, control.m1.bc, class_factor, output="fpr")
  
  fpr.b.m2 <- tprfpr(batch_df.bc, beam.m2.bc, class_factor, output="fpr")
  fpr.c.m2 <- tprfpr(batch_df.bc, control.m2.bc, class_factor, output="fpr")
  
  fpr.b.m1.s <-tprfpr(true.m1.s, beam.m1.bc.s, class_factor, output="fpr")
  fpr.c.m1.s <-tprfpr(true.m1.s, control.m1.bc.s, class_factor, output="fpr")
  
  fpr.b.m2.s <-tprfpr(true.m2.s, beam.m2.bc.s, class_factor, output="fpr")
  fpr.c.m2.s <-tprfpr(true.m2.s, control.m2.bc.s, class_factor, output="fpr")
  
  fpr_vals <- c(fpr.b.m1, fpr.c.m1, fpr.b.m2, fpr.c.m2, fpr.b.m1.s, fpr.c.m1.s,
                fpr.b.m2.s, fpr.c.m2.s)
  fpr_res <- rbind(fpr_res, fpr_vals)
  
  # Check current iteration
  print(paste0("Current iteration:",i))
}

# Set column names and plot colours
colnames(rmse_res) <- c('BEAM.M1', 'Control.M1', 'BEAM.M2', 'Control.M2','BEAM.M1.r','Control.M1.r','BEAM.M2.r','Control.M2.r')
colnames(feature_res) <- c('BEAM.M1', 'Control.M1', 'BEAM.M2', 'Control.M2','M1.r','M2.r')
colnames(tpr_res) <- c('BEAM.M1', 'Control.M1', 'BEAM.M2', 'Control.M2','BEAM.M1.r','Control.M1.r','BEAM.M2.r','Control.M2.r')
colnames(fpr_res) <- c('BEAM.M1', 'Control.M1', 'BEAM.M2', 'Control.M2','BEAM.M1.r','Control.M1.r','BEAM.M2.r','Control.M2.r')
mycolor<-c("#CDF0EA","#6E85B7","#A8E890","#42855B","yellow3","yellow4","pink","pink4")

# NRMSE boxplot
rmse_plot<-stack(as.data.frame(rmse_res))
colnames(rmse_plot)<-c('values','legend')
ggplot(rmse_plot, aes(x = legend, y = values, fill = legend)) + theme_classic() + scale_fill_manual(values=mycolor) +
  geom_boxplot() + labs(x='Imputation method', y='NRMSE') + theme(axis.text.x=element_blank())

# Features boxplot
feature_plot<-stack(as.data.frame(feature_res))
colnames(feature_plot)<-c('values','legend')
ggplot(feature_plot, aes(x = legend, y = values, fill = legend)) + theme_classic() + scale_fill_manual(values=mycolor) +
  geom_boxplot() + labs(x='Imputation method', y='no. of features') + theme(axis.text.x=element_blank())

# TPR boxplot
tpr_plot<-stack(as.data.frame(tpr_res))
colnames(tpr_plot)<-c('values','legend')
ggplot(tpr_plot, aes(x = legend, y = values, fill = legend)) + theme_classic() + scale_fill_manual(values=mycolor) +
  geom_boxplot() + labs(x='Imputation method', y='TPR') + theme(axis.text.x=element_blank())

# FPR boxplot
fpr_plot<-stack(as.data.frame(fpr_res))
colnames(fpr_plot)<-c('values','legend')
ggplot(fpr_plot, aes(x = legend, y = values, fill = legend)) + theme_classic() + scale_fill_manual(values=mycolor) +
  geom_boxplot() + labs(x='Imputation method', y='FPR') + theme(axis.text.x=element_blank())


# Boxplot of intra-sample variance (using dataset from final iteration)
require(grid)
require(gridExtra)
var.tm1.1 <- var(true.m1.s)
var.tm1 <- melt(var.tm1.1)

var.tm2.1 <- var(true.m2.s)
var.tm2 <- melt(var.tm2.1)

var.bm1.1 <- var(beam.m1.bc.s)
var.bm1 <- melt(var.bm1.1)

var.cm1.1 <- var(control.m1.bc.s)
var.cm1 <- melt(var.cm1.1)

var.bm2.1 <- var(beam.m2.bc.s)
var.bm2 <- melt(var.bm2.1)

var.cm2.1 <- var(control.m2.bc.s)
var.cm2 <- melt(var.cm2.1)

v1=ggplot() + geom_boxplot(data=var.tm1, aes(x=Var1, y=value),colour='black', fill='lightgreen') + theme_bw() + theme(axis.text.x=element_blank()) + labs(title="true M1", x='Samples', y='Variance') + ylim(5,13)
v2=ggplot() + geom_boxplot(data=var.tm2, aes(x=Var1, y=value),colour='black', fill='lightgreen') + theme_bw() + theme(axis.text.x=element_blank()) + labs(title="true M2", x='Samples', y='Variance') + ylim(3,12)
v3=ggplot() + geom_boxplot(data=var.bm1, aes(x=Var1, y=value),colour='black', fill='lightgreen') + theme_bw() + theme(axis.text.x=element_blank()) + labs(title="BEAMs M1", x='Samples', y='Variance') + ylim(5,13)
v4=ggplot() + geom_boxplot(data=var.cm1, aes(x=Var1, y=value),colour='black', fill='lightgreen') + theme_bw() + theme(axis.text.x=element_blank()) + labs(title="Control M1", x='Samples', y='Variance') + ylim(5,13)
v5=ggplot() + geom_boxplot(data=var.bm2, aes(x=Var1, y=value),colour='black', fill='lightgreen') + theme_bw() + theme(axis.text.x=element_blank()) + labs(title="BEAMs M2", x='Samples', y='Variance') + ylim(3,12)
v6=ggplot() + geom_boxplot(data=var.cm2, aes(x=Var1, y=value),colour='black', fill='lightgreen') + theme_bw() + theme(axis.text.x=element_blank()) + labs(title="Control M2", x='Samples', y='Variance') + ylim(3,12)

combined_var = grid.arrange(arrangeGrob(v1,v2,v3,v4,v5,v6, ncol=3, layout_matrix = rbind(c(1,3,4),c(2,5,6))))
print(combined_var)


# Check statistical significance
p.rmse <- check.p(rmse_res)
p.features <- check.p(feature_res)
p.tpr <- check.p(tpr_res)
p.fpr <- check.p(fpr_res)

# Check values of extra features in BEAM M1
test1 = rownames(beam.m1)
test2 = rownames(beam.m1.bc.s)
test3 = match(test2, test1)
test4 = test1[-test3]
for (i in 1:length(test4)){
  print(test4[i] %in% test2)
}

mis.1 = beam.df[as.numeric(test4),]
true.1 = batch_df[as.numeric(test4),]
imp.1 = beam.m1[test4,]

imp.t = as.numeric(imp.1)
true.t = as.numeric(true.1)

dat = cbind(imp=imp.t, true=true.t, value=seq(length(imp.t)))

ggplot(data=dat) + geom_point(aes(x=value, y=imp), color = "red2", alpha = 0.5) +
  geom_point(aes(x=value, y=true), color = "blue2", alpha = 0.5) +labs(y="Value", x="Data point") +
  geom_vline(xintercept = (length(imp.t)/2), linetype = "dashed") + theme_bw()

###### Additional PCA check

# M1 PCA
#pc<-t(beam.m1.bc.s)
#pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
#pca_res<-prcomp(pc, scale. = FALSE)
#autoplot(pca_res, data = pc1, colour = 'batch') + theme_bw()

#pc<-t(control.m1.bc.s)
#pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
#pca_res<-prcomp(pc, scale. = FALSE)
#autoplot(pca_res, data = pc1, colour = 'batch') + theme_bw()

#pc<-t(true.m1.s)
#pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
#pca_res<-prcomp(pc, scale. = FALSE)
#autoplot(pca_res, data = pc1, colour = 'batch') + theme_bw()


# M2 PCA
#pc<-t(beam.m2.bc.s)
#pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
#pca_res<-prcomp(pc, scale. = FALSE)
#autoplot(pca_res, data = pc1, colour = 'batch') + theme_bw()

#pc<-t(control.m2.bc.s)
#pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
#pca_res<-prcomp(pc, scale. = FALSE)
#autoplot(pca_res, data = pc1, colour = 'batch') + theme_bw()

#pc<-t(true.m2.s)
#pc1<-cbind(pc, batch=c(rep('batch 1',10),rep('batch 2',10)))
#pca_res<-prcomp(pc, scale. = FALSE)
#autoplot(pca_res, data = pc1, colour = 'batch') + theme_bw()
