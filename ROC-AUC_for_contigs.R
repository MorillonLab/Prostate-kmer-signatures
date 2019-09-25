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
# list of prostate specific probes 
frame.probes <-as.data.frame(read.delim("Data_for_DEkupl/list_mixed.csv",sep=",",header=T,check.names=T))
listProbesMixed <-frame.probes$Probes #list of selection probes

# list of lncRNA probes
frame.probes.lncRNA <-as.data.frame(read.delim("Data_for_DEkupl/list_new_lnc.csv",sep=",",header=T,check.names=T))
listProbesLncRNA <-frame.probes.lncRNA$Probes #list of selection probes

# loading nanoString dataset
dataNanoString <-as.data.frame(read.delim("Data_for_DEkupl/nanoString_count.csv",sep=",",header=T,check.names=T)) 
dataNanoString <- dataNanoString[,-1]

# loading TCGA dataset
dataTCGA <-as.data.frame(read.delim("Data_for_DEkupl/TCGA_count.csv",sep=",",header=T,check.names=T)) 
dataTCGA <- dataTCGA[,-1]

# folder to save result
dirResult = "Data_for_DEkupl/Result_contigs/"
dir.create(file.path(dirResult), showWarnings = FALSE, recursive = TRUE)
## ------------------------------------------------------------------------


# Restriction of the dataset to prostate specific probes
dataNanoStringMixed <- dataNanoString[which(colnames(dataNanoString)%in%listProbesMixed
                                            | colnames(dataNanoString)%in%c("condition","risk"))]

# Signature selected on NanoString for prostate specific probes including PCA3
signature_mixed <- extractSignature(dataNanoStringMixed)
signature_mixed

# Restriction of the dataset to prostate specific probes
dataNanoStringLncRNAProbes <- dataNanoString[
  which(colnames(dataNanoString)%in%listProbesLncRNA
        | colnames(dataNanoString)%in%c("condition","risk"))]

# Signature selected on NanoString using lncRNA contigs only 
signature_new_lnc <- extractSignature(dataNanoStringLncRNAProbes)
signature_new_lnc


write.csv2(signature_mixed,file=paste0(dirResult,"signature_mixed.txt"))
write.csv2(signature_new_lnc,file=paste0(dirResult,"signature_new_lnc.txt"))


# selelection NanoString validation NanoString minimal
# 
listSignature_NanoString <- list("PCA3",signature_mixed,signature_new_lnc)


NUM_RUNS=100
# Concatenate mean AUC and prepare data for drawing ROC curve
resSign1 <- concatMeanAUCReturnFrameROC(dataNanoString, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "normal-tumor")

resSign2 <- concatMeanAUCReturnFrameROC(dataNanoString, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "HR")

resSign3 <- concatMeanAUCReturnFrameROC(dataNanoString, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "IR")

resSign4 <- concatMeanAUCReturnFrameROC(dataNanoString, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "LR")

resSign5 <- concatMeanAUCReturnFrameROC(dataNanoString, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "combine")

resSign6 <- concatMeanAUCReturnFrameROC(dataNanoString, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "LH")


# Concatenate all rows to the final table
finalTab <- rbind(resSign1[[1]],resSign2[[1]],resSign3[[1]],
                  resSign4[[1]],resSign5[[1]],resSign6[[1]])

rownames(finalTab) <- c("Mean Normal vs Tumor", "SD Normal vs Tumor",
                        "Mean Normal vs  HR", "SD Normal vs  HR",
                        "Mean Normal vs  IR","SD Normal vs  IR",
                        "Mean Normal vs  LR","SD Normal vs  LR",
                        "Mean Normal+LR vs HR+IR", "SD Normal+LR vs HR+IR",
                        "Mean LR vs HR", "SD LR vs HR")

finalTab

save.image(paste0(dirResult, "selNanoString_valNanoString.RData"))

## ------------------------------------------------------------------------

write.csv(finalTab, paste0(dirResult, "performance_selNanoString_valNanoString.csv"), row.names = TRUE)
write.table(finalTab, paste0(dirResult, "performance_selNanoString_valNanoString.txt"))

# draw ROC curve
png(filename = paste0(dirResult, "ROC_signature_selNanoString_valNanoString.png"),width=800,height=800)
figure <- drawROCcurve(resSign1[[3]],"Normal vs Tumor", mode="NanoString")
plot(figure)
dev.off()

## selelection NanoString validation TCGA

# If ctg_81545_37852 is in a signature (from NanoString dataset)
# replace it with two contigs : ctg_81545 and ctg_37852 (in TCGA dataset)
for (i in 1:length(listSignature_NanoString)){for (j in 1:length(listSignature_NanoString[[i]])){
  if (listSignature_NanoString[[i]][[j]] == "ctg_81545_37852"){
    listSignature_NanoString[[i]][[j]] <- "ctg_81545"
    listSignature_NanoString[[i]]<- append(listSignature_NanoString[[i]], "ctg_37852")
  }
}
}

NUM_RUNS=100
# Concatenate mean AUC and prepare data for drawing ROC curve
resSign1 <- concatMeanAUCReturnFrameROC(dataTCGA, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "normal-tumor")

resSign2 <- concatMeanAUCReturnFrameROC(dataTCGA, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "HR")

resSign3 <- concatMeanAUCReturnFrameROC(dataTCGA, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "IR")

resSign4 <- concatMeanAUCReturnFrameROC(dataTCGA, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "LR")

resSign5 <- concatMeanAUCReturnFrameROC(dataTCGA, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "combine")

resSign6 <- concatMeanAUCReturnFrameROC(dataTCGA, NUM_RUNS, listSignature_NanoString,
                                        oneCondition = "LH")

# Concatenate all rows to the final table
finalTab <- rbind(resSign1[[1]],resSign2[[1]],resSign3[[1]],
                  resSign4[[1]],resSign5[[1]],resSign6[[1]])

rownames(finalTab) <- c("Mean Normal vs Tumor", "SD Normal vs Tumor",
                        "Mean Normal vs  HR", "SD Normal vs  HR",
                        "Mean Normal vs  IR","SD Normal vs  IR",
                        "Mean Normal vs  LR","SD Normal vs  LR",
                        "Mean Normal+LR vs HR+IR", "SD Normal+LR vs HR+IR",
                        "Mean LR vs HR", "SD LR vs HR")

finalTab

save.image(paste0(dirResult, "selNanoString_valTCGA.RData"))

## ------------------------------------------------------------------------

write.csv(finalTab, paste0(dirResult, "performance_signature_selNanoString_valTCGA.csv"), row.names = TRUE)
write.table(finalTab, paste0(dirResult, "performance_signature_selNanoString_valTCGA.txt"))

# draw ROC curve 
png(filename = paste0(dirResult, "ROC_signature_selNanoString_valTCGA.png"),width=800,height=800)
figure <- drawROCcurve(resSign1[[3]],"Normal vs Tumor", mode="TCGA")
plot(figure)
dev.off()

## infoSessions
sessionInfo()