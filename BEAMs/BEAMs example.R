#### Uncovering the Consequences of BEAMs in Omics Data Analysis

# BEAMs example

library(here)
library(imputeLCMD)



wd <- here()
source(paste0(wd, "/BEAMs/impute functions.r"))

set.seed(2023)
class_factor <- c(rep(c(rep(0,3), rep(1,3), rep(2,3)),3))
batch_factor <- c(rep(1,9),rep(2,9), rep(3,9))
Z = 1.2
Y = sqrt(5)
mv_prop1=c(0.15,0.0475)
mv_prop2=c(0.3,0.075)
mv_prop3=c(0.5,0.135)
mv_prop4=c(0.3,0.075)

############## Data simulation
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

# Heatmap visualization
visdat::vis_miss(as.data.frame(mis1))
visdat::vis_miss(as.data.frame(mis2))
