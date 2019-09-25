################
# Author : Thi Ngoc Ha Nguyen
# Date   : 08/20/2019
# Email  : thi-ngoc-ha.nguyen@i2bc.paris-saclay.fr
################

rm(list=ls())
set.seed(12678)

################
# Load libraries
################

library(knitr)

################
# Load list of functions 
################
source("useful_functions.R")
#######################################################################
## scripts
#######################################################################
## ------------------------------------------------------------------------
# loading small clinical cohort
dataNanoString <- as.data.frame(read.delim("Data_for_genes/DESeq_UP_pvalue05_normal_vs_tumor.csv",sep=",",header=T,check.names=F)) 
dataNanoString <- data.frame(ID=dataNanoString$gene_name, dataNanoString[-2], check.names=F)[,-c(2:7, 32:36,38,39)]

# keeping genes with log2FC >2
dataNanoString <- dataNanoString[which(dataNanoString$log2FoldChange>2),]
# folder to save result
dirResult = "Data_for_genes/Result_genes/"
dir.create(file.path(dirResult), showWarnings = FALSE, recursive = TRUE)
## ------------------------------------------------------------------------
# transpose
geneID <- dataNanoString$ID

dataNanoString <- as.data.frame(t(dataNanoString[,-1]))

colnames(dataNanoString) <- geneID

for (i in 1:nrow(dataNanoString)){
  
  if(as.character(substr(rownames(dataNanoString)[i],start = 1,stop = 6)) == "normal"){
    
    dataNanoString$condition[i] <- "normal"
    
  }else{
    
    dataNanoString$condition[i] <- "tumoral"
    
  }
  
}

dataNanoString$condition <- factor(dataNanoString$condition)

# loading TCGA dataset
dataTCGA <-as.data.frame(read.delim("Data_for_genes/genes_expression_tcga.tsv",sep="\t",header=T,check.names=F))
#save(dataTCGA,file = paste0(dirResult,"gene_expression_tcga.RData"))
dataTCGA <- data.frame(ID=dataTCGA$gene_name, dataTCGA[-2], check.names=F)[,-c(2:7)]

# transpose
geneID <- dataTCGA$ID
dataTCGA <- as.data.frame(t(dataTCGA[,-1]))

colnames(dataTCGA) <- geneID


# loading sample and corresponding condition
dataAnnotation <-as.data.frame(read.delim("Data_for_genes/sample_condition.tsv",sep="\t",header=T,check.names=F))
# map sample count and corresponding condition
dataTCGA$ID_kallisto <- rownames(dataTCGA)

dataTCGA<-merge(dataAnnotation,dataTCGA,
                 by.x="ID_kallisto",
                 by.y="ID_kallisto",
                 all.x=FALSE,
                 all.y=FALSE)

dataTCGA <- dataTCGA[,-1]

## ------------------------------------------------------------------------

# Signature selected on NanoString using stable lasso regression

signature.stb <- extractSignatureStb(dataNanoString, 0.5)
signature.stb

write.csv2(signature.stb,file=paste0(dirResult,"signature-gene.txt"))


NUM_RUNS=100

# dataframe of signature in TCGA dataset
dataFrameSigTCGA_NT <- dataTCGA[c("condition", signature.stb)]

dataFrameSigTCGA_NH <- extract.risk(dataTCGA, risk.level = "HR", signature.stb)

dataFrameSigTCGA_NI <- extract.risk(dataTCGA, risk.level = "IR", signature.stb)

dataFrameSigTCGA_NL <- extract.risk(dataTCGA, risk.level = "LR", signature.stb)


# Calculate mean and sd AUC of gene signature in each condition

resSign1 <- takeDataReturnAUCvalues(dataFrameSigTCGA_NT, NUM_RUNS)

resSign2 <- takeDataReturnAUCvalues(dataFrameSigTCGA_NH, NUM_RUNS)

resSign3 <- takeDataReturnAUCvalues(dataFrameSigTCGA_NI, NUM_RUNS)

resSign4 <- takeDataReturnAUCvalues(dataFrameSigTCGA_NL, NUM_RUNS)

auc1 <- do.call("rbind", lapply(resSign1, "[[", 1))
auc2 <- do.call("rbind", lapply(resSign2, "[[", 1))
auc3 <- do.call("rbind", lapply(resSign3, "[[", 1))
auc4 <- do.call("rbind", lapply(resSign4, "[[", 1))

row1 <- rbind(round(mean(auc1),2),round(sd(auc1),2)) 
row2 <- rbind(round(mean(auc2),2),round(sd(auc2),2)) 
row3 <- rbind(round(mean(auc3),2),round(sd(auc3),2)) 
row4 <- rbind(round(mean(auc4),2),round(sd(auc4),2)) 

# Concatenate all rows to the final table
finalTab <- rbind(row1, row2, row3, row4)

rownames(finalTab) <- c("Mean Normal vs Tumor", "SD Normal vs Tumor",
                        "Mean Normal vs  HR", "SD Normal vs  HR",
                        "Mean Normal vs  IR","SD Normal vs  IR",
                        "Mean Normal vs  LR","SD Normal vs  LR")
 
finalTab

save.image(paste0(dirResult,"gene_selNanoString_valTCGA.RData"))

## ------------------------------------------------------------------------

write.csv(finalTab, paste0(dirResult,"performance_geneSignature_selNanoString_valTCGA.csv"), row.names = TRUE)
write.table(finalTab, paste0(dirResult,"performance_geneSignature_selNanoString_valTCGA.txt"))