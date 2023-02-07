library(here)
library(devtools)
library(MSqRob)
library(Biobase)
library(MSnbase)
library(limma)
library(ggplot2)
library(dplyr)
library(pROC)
library(sva)
#library(proBatch)
library(vegan)
library(grid)
library(gridExtra)
library(missForest)
#BiocManager::install("imputeLCMD")
library(imputeLCMD)

#library(MSqRobData)
wd <- here()
source(paste0(wd, "/BEAMs/CPTAC/functions.r"))
source(paste0(wd, "/BEAMs/CPTAC/DEP.r"))
source(paste0(wd, "/BEAMS/impute functions.r"))

########
correct_batch = TRUE
# Set results
pAUC_res<-data.frame(matrix(nrow=0,ncol=6))
k_plot<-c()
#####################################

# peptides CPTAC data
load(paste0(wd, "/BEAMs/CPTAC/peptides_CPTAC.RData"))

# proteingroups file with manually removed ups proteins labelled as contaminants
dataset.CPTAC <- read.table(paste0(wd,"/BEAMs/CPTAC/proteinGroups_cont_curated.txt"),
                            sep="\t", quote="", comment.char = "", header=TRUE)

# Remove reverse sequences
dataset.CPTAC <- dplyr::filter(dataset.CPTAC, Reverse != "+", Potential.contaminant != "+", Only.identified.by.site != "+")

# Get unique proteins
dataset.CPTAC$Gene.names %>% duplicated() %>% any()
## TRUE

# No duplicated protein IDs, as expected
dataset.CPTAC$Protein.IDs %>% duplicated() %>% any()
## FALSE

dataset.CPTAC <- DEP::make_unique(dataset.CPTAC, "Protein.IDs", "Protein.IDs", delim = ";")

# Make experimental design
experimental.design <- data.frame(label = peptides.CPTAC %>% sampleNames %>% as.character, condition = peptides.CPTAC %>% sampleNames %>% substr(2,3) %>% as.character, replicate = c("LTQ-Orbitrap_86" %>% rep(3), "LTQ-OrbitrapO_65" %>% rep(3), "LTQ-OrbitrapW_56" %>% rep(3)) %>% rep(3) %>% as.character, stringsAsFactors = FALSE)

class_factor = c(rep(0,9),rep(1,9),rep(2,9))
batch_factor = rep(c(rep(1,3),rep(2,3),rep(3,3)),3)

protein.names = dataset.CPTAC$Protein.IDs
protein.names = str_split(protein.names, pattern = ";")
prot_names = sapply(protein.names, "[[", 1)

LFQ.columns <- grep("LFQ.", colnames(dataset.CPTAC)) # get LFQ column numbers
batch_df <- dataset.CPTAC[,LFQ.columns]
rownames(batch_df) <- prot_names
batch_df <- log2(batch_df)
batch_df <- as.matrix(batch_df)
batch_df[which(is.infinite(batch_df))]<-NA

# Check missingness
require(visdat)
vis_miss(as.data.frame(batch_df))

# Normalize the data
batch_df.norm <- normalize.quantiles(batch_df, keep.names = TRUE)

# Imputation by M1 and M2
M1.imp <- impute.M1(batch_df, 0.8)
rownames(M1.imp) <- prot_names[as.numeric(rownames(M1.imp))]
M2.imp <- impute.M2(batch_df, batch_factor, 0.8)
rownames(M2.imp) <- prot_names[as.numeric(rownames(M2.imp))]

batch1<- batch_df[,which(batch_factor==1)]
batch2<- batch_df[,which(batch_factor==2)]
batch3<- batch_df[,which(batch_factor==3)]
batch.all <- cbind(batch1,batch2,batch3)
vis_miss(as.data.frame(batch.all))

# Check number of UPS proteins retained
length(which(grepl('ups', prot_names)))
length(which(grepl('ups', rownames(M1.imp))))
length(which(grepl('ups', rownames(M2.imp))))
