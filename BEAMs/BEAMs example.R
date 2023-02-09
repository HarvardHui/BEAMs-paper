#### Uncovering the Consequences of BEAMs in Omics Data Analysis
# Author: Harvard Wai Hann Hui
# Date: 09/February/2023
# Script: BEAMs example

############################################################################################################################
#### This Script creates example datasets with and without BEAMs using three batches
#### As this is purely for visualisation purposes, batch effects are not simulated
####
#### Please check the "Set Parameters" section to ensure that parameters are
#### set correctly, or as desired
############################################################################################################################


###### Load all necessary libraries and functions
library(here)
library(imputeLCMD)

wd <- here()
source(paste0(wd, "/BEAMs/impute functions.r"))

###### Set parameters
set.seed(2023)
batch_factor <- c(rep(1,9),rep(2,9), rep(3,9))

# BEAMs proportion
mv_prop1=c(0.15,0.0475)
mv_prop2=c(0.3,0.075)
mv_prop3=c(0.5,0.135)

# Control MVs proportion
mv_prop4=c(0.3,0.075)

############## Data simulation (Functions from Shah et al. 2017)
CorrelationMatrix0.7 <- CorrMatrixNegFixed(30, 30, 0.7, 0.2)
DataSimul <- SimulatedData(27, 900, as.matrix(nearPD(CorrelationMatrix0.7)[[1]]), low = 5, high =15)
DataSimul <- t(DataSimul)

##### Non-BEAMs simulation
mis1 = mv.sim(DataSimul, mv_prop4[1],mv_prop4[2])
mis1 = mis1[[2]]

##### BEAMS simulation
mis2.1 = mv.sim(DataSimul[,which(batch_factor==1)], mv_prop1[1], mv_prop1[2])
mis2.2 = mv.sim(DataSimul[,which(batch_factor==2)], mv_prop2[1], mv_prop2[2])
mis2.3 = mv.sim(DataSimul[,which(batch_factor==3)], mv_prop3[1], mv_prop3[2])

mis2 = cbind(mis2.1[[2]], mis2.2[[2]], mis2.3[[2]])

# Visualization by heatmap
visdat::vis_miss(as.data.frame(mis1))
visdat::vis_miss(as.data.frame(mis2))
