#### Uncovering the Consequences of BEAMs in Omics Data Analysis
# Author: Harvard Wai Hann Hui
# Date: 09/February/2023
# Script: CPTAC study 6 BEAMs example


############################################################################################################################
#### This Script visualizes the BEAMs present in the CPTAC study 6 dataset,
#### and the feature changes after M1 and M2 MVI.
#### Processing steps were adapted from Ludger et al. 2020, which was based on the
#### introduction to DEP vignette.
####
#### Please check the "Set Parameters" section to ensure that parameters are
#### set correctly, or as desired
############################################################################################################################


###### Load all necessary libraries and functions
library(here)
library(stringr)
library(dplyr)
library(visdat)
library(imputeLCMD)
library(preprocessCore)

wd <- here()
source(paste0(wd, "/BEAMS/impute functions.r"))

###### Set Parameters
batch_factor = rep(c(rep(1,3),rep(2,3),rep(3,3)),3)
k = 3
mv.threshold = 0.8

################## CPTAC data processing adapted from Ludger et al. 2020

######## peptides CPTAC data
load(paste0(wd, "/BEAMs/CPTAC/peptides_CPTAC.RData"))

######## proteingroups file with manually removed ups proteins labelled as contaminants
dataset.CPTAC <- read.table(paste0(wd,"/BEAMs/CPTAC/proteinGroups_cont_curated.txt"),
                            sep="\t", quote="", comment.char = "", header=TRUE)

######## Remove reverse sequences
dataset.CPTAC <- dplyr::filter(dataset.CPTAC, Reverse != "+", Potential.contaminant != "+", Only.identified.by.site != "+")

######## Get unique proteins
dataset.CPTAC$Gene.names %>% duplicated() %>% any()
## TRUE

######## No duplicated protein IDs, as expected
dataset.CPTAC$Protein.IDs %>% duplicated() %>% any()
## FALSE

dataset.CPTAC <- DEP::make_unique(dataset.CPTAC, "Protein.IDs", "Protein.IDs", delim = ";")

######## Get expression matrix with protein names
protein.names = dataset.CPTAC$Protein.IDs
protein.names = str_split(protein.names, pattern = ";")
prot_names = sapply(protein.names, "[[", 1)

LFQ.columns <- grep("LFQ.", colnames(dataset.CPTAC)) # get LFQ column numbers
batch_df <- dataset.CPTAC[,LFQ.columns]
rownames(batch_df) <- prot_names
batch_df <- log2(batch_df)
batch_df <- as.matrix(batch_df)
batch_df[which(is.infinite(batch_df))]<-NA

######## Check missingness
vis_miss(as.data.frame(batch_df))

# Arrange by batch
batch1<- batch_df[,which(batch_factor==1)]
batch2<- batch_df[,which(batch_factor==2)]
batch3<- batch_df[,which(batch_factor==3)]
batch.all <- cbind(batch1,batch2,batch3)
vis_miss(as.data.frame(batch.all))

######## Imputation by M1 and M2
M1.imp <- impute.M1(batch_df, k=k, mv.threshold = mv.threshold)
rownames(M1.imp) <- prot_names[as.numeric(rownames(M1.imp))]
M2.imp <- impute.M2(batch_df, k=k, batch_factor, mv.threshold = mv.threshold)
rownames(M2.imp) <- prot_names[as.numeric(rownames(M2.imp))]

####### Check number of UPS proteins retained
# Total UPS proteins
length(which(grepl('ups', prot_names)))
## 40

# UPS proteins after M1 MVI
length(which(grepl('ups', rownames(M1.imp))))
## 30

# UPS proteins after M2 MVI
length(which(grepl('ups', rownames(M2.imp))))
## 9
