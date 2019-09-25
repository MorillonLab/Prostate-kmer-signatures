################
# Author : Thi Ngoc Ha Nguyen
# Date   : 04/10/2019
# Email  : thi-ngoc-ha.nguyen@i2bc.paris-saclay.fr
################

rm(list=ls())
set.seed(12678)

################
# Load libraries
################

library(glmnet)
library(caTools)
library(caret)
library(precrec)
library(zoo)
library(parallel)

################
# Load list of functions 
################

################
# Functions for probe selection
################
## A) Functions to extract subgroup from the dataset 

extract.risk <- function(frame, risk.level, top){
  # Extract data for specific risk condition
  #     
  # Args:
  #   frame     : Data frame containes all the data
  #   risk.level: A string specify the level or risk condition (ex: "HR" for high risk ...)
  #   top       : List of probes are selected by Lasso logistic regression
  #
  # Returns:
  #   Data frame of selected probes in normal vs risk.level
  
  top_data <- subset(frame, ((frame$risk == risk.level & frame$condition == "tumoral") | frame$condition == "normal"), select = c("condition",top))
  
  top_data$condition <- as.factor(as.character(top_data$condition))
  
  return (top_data)
  
}

extract.combine <- function(frame, top){
  # Extract data for combination condition.
  #
  # Args:
  #   frame: Data frame containes all the data.
  #   top  : List of probes are selected by Lasso logistic regression.
  # 
  # Returns:
  #   Data frame of selected probes in combination condition.
  
  # Remove all rows with risk level equals not_avaiable
  frame <- subset(frame, (frame$risk=="HR"|frame$risk=="IR"|frame$risk=="LR"|frame$condition=="normal"))
  
  frame <- transform(frame, condition= ifelse(frame$condition=="normal"|frame$risk=="LR", "normal", "tumoral"))
  
  frame <- frame[,c("condition",top)]
  
  return (frame)
  
}

extract.condition <- function(frame, cond1, cond2, top){
  # Extract data for two conditions
  #     
  # Args:
  #   frame : Data frame containes all the data
  #   cond1 : A string specify the first condition
  #   cond2: A string specify the other condition
  #   top       : List of probes are selected by Lasso logistic regression
  #
  # Returns:
  #   Data frame of selected probes in cond1 vs cond2
  
  frame <- subset(frame, ((frame$risk==cond1&frame$condition=="tumoral")|(frame$risk==cond2&frame$condition=="tumoral")), select = c("risk",top))
  
  frame <- transform(frame, condition= ifelse(frame$risk==cond1, "normal", "tumoral"))
  
  frame <- frame[which(!colnames(frame)%in%c("risk"))]
  
  frame$condition <- as.factor(as.character(frame$condition))
  
  frame <-frame[c("condition",setdiff(names(frame),"condition"))]
  
  return (frame)
  
}

################
## B) Function to extract signature
################
extractSignature <- function(frame, incPCA3=FALSE){
  # Extract signature
  #
  # Args:
  #   frame : Data frame containes selected probes
  #   incPCA3 : whether including PCA3 (TRUE: including PCA3 and vise versa) 
  #
  # Returns:
  #   The signature of the selected probes
  
  if(incPCA3 == FALSE){
    
    frame <- frame[which(!colnames(frame)%in%c("PCA3"))]
    
  }
  
  # Extract dataframe for using Lasso logistic regression
  frame.lasso <-frame[which(!colnames(frame)%in%c("risk"))] 
  
  # Perform upsampling for correcting unbalanced design
  up_data <- upSample(x = frame.lasso[,-1],
                      y = frame.lasso$condition)
  frame.lasso.up <- cbind(condition=up_data$Class,up_data[,-ncol(up_data)])
  
  # Convert training data to matrix format
  x <- model.matrix(condition~.,frame.lasso.up)
  
  # Perform k-fold cross validation for glmnet to find optimal value of lambda
  cvGlm <- cv.glmnet(x = x, y = frame.lasso.up[,1], family = "binomial", alpha =1, type.measure = "mse")
  
  sel <- glmnet(x = x, y = frame.lasso.up[,1],lambda=cvGlm$lambda.1se, family = "binomial", alpha = 1)
  
  id <- which( abs( as.vector(sel$beta) ) > 0)
  
  # List of signature selected by Lasso logistic regression
  signature <- names(frame.lasso.up)[id]
  
  return(signature)
  
}
################
## C) Function to extract stable signature
################
extractSignatureStb <- function(frame, thres){
  # Extract signature
  #
  # Args:
  #   frame : Data frame containes selected probes
  #
  # Returns:
  #   The signature of the stable selected probes
  
  # Perform upsampling for correcting unbalanced design
  up_data <- upSample(x = frame[,-length(frame)],
                      y = frame$condition)
  frame.up <- cbind(condition=up_data$Class,up_data[,-ncol(up_data)])
  
  # Convert training data to matrix format
  X <- model.matrix(condition~.,frame.up)
  
  Y <- frame.up$condition
  
  n <- nrow(X) 
  
  p <- ncol(X)
  
  Gene <- colnames(X) 
  
  # Perform k-fold cross validation for glmnet to find optimal value of lambda
  cvGlm <- cv.glmnet(x = X, y = Y, family = "binomial", alpha =1, type.measure = "mse")
  
  cores <- detectCores()
  
  cl <- makeCluster(cores-2)
  
  NUM_RUNS <- 2e3
  
  stabsel <- function(i){
    cat("+")
    b_sort <- sort(sample(1:n,round(3*n/4)))
    out <- glmnet(X[b_sort,],Y[b_sort], family = "binomial",
                  lambda=cvGlm$lambda.1se, alpha=1)
    return(tabulate(which(out$beta[,ncol(out$beta)]!=0),p))
  }
  
  clusterEvalQ(cl, expr = c(library(glmnet)))
  clusterExport(cl,c('stabsel','frame', 'NUM_RUNS','n','glmnet','cvGlm','X','Y','p'),envir=environment())
  
  res.cum <- Reduce("+", parLapply(cl, 1:NUM_RUNS, stabsel))
  
  stopCluster(cl)
  
  prob.sel <- res.cum/NUM_RUNS
  plot(sort(prob.sel))
  
  gene.stabsel <- Gene[prob.sel >= thres]
  
  return (gene.stabsel)
}

################
## D) Function to calculate performance of signature in specific condition
################
takeDataReturnAUCvalues<- function(frame, NUM_RUNS, mode){
  # Calculate mean and sd for signature in specific condition
  #
  # Args:
  #   frame: Data frame containes selected probes and specific conditon
  #   NUM_RUNS : number of time we sample from the data frame
  #   method : string indicates method for estimating model accuracy (K-fold cross validation or leave one out cross validation -LOOCV)
  #
  # Results:
  #   AUC ROC curve in specific condition
  #   List predicted scores and observed labels by the signature
  
  n <- detectCores()
  
  cl <- makeCluster(n-2)
  
  calculateAUC <- function(i){
    
    # Sampling the data into training and testing set
    # Extract data labels
    Y = frame$condition
    
    # Create logical vector of the same length as Y 
    # with random SplitRatio*length(Y) elements set to TRUE
    msk = sample.split(Y, SplitRatio=0.7)
    
    imbal.train <- as.data.frame(frame[ msk,])
    
    imbal.test <- as.data.frame(frame[!msk,])
    
    # Train model 
    ctrl <- trainControl(method = "LOOCV",
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary,
                         sampling = "up")
    
    # Use Boosted Logistic Regression method in caret package
    up_inside <- train(condition ~ ., data = imbal.train,
                       method = "LogitBoost",
                       nIter = 100,
                       metric = "ROC",
                       trControl = ctrl)  
    
    # Predition for test data
    pred = predict(up_inside, imbal.test[,-1, drop = FALSE], type = "prob")[, "tumoral"]
    
    # Calculate auc for ROC curve
    resAUC <- auc(evalmod(scores = pred, labels = imbal.test$condition))
    
    rocAUC <- resAUC$aucs[1]
    
    if(rocAUC<0.5){ # Reverse roc value in case it lower than 0.5
      
      rocAUC = 1 - rocAUC
      
    }
    
    # Save predicted scores 
    scores <- pred
    
    # Save observed labels
    label <- imbal.test$condition
    
    return(list(rocAUC, scores, label))
  }
  
  clusterEvalQ(cl, expr = c(library(caTools),library(caret), library(precrec)))
  clusterExport(cl,c('calculateAUC','frame', 'NUM_RUNS','sample.split','trainControl','train','auc','evalmod'),envir=environment())
  
  res <- parLapply(cl, 1:NUM_RUNS, calculateAUC)
  
  stopCluster(cl)
  
  return(res)
  
}  

################
## E) Function to concatenate mean and prepare data to draw ROC curve for signatures in one condition
################
concatMeanAUCReturnFrameROC <- function(frame, NUM_RUNS, signature, oneCondition){
  # Concatenate mean and sd for all signatures in one conditions
  # and calculate data for drawing ROC curve
  #
  # Args:
  #   frame        : Data frame containes all the data.
  #   NUM_RUNS     : Number of times running
  #   signature    : List of signature
  #   oneCondition : Specific condition
  #
  # Results:
  #   List of mean and standard deviation AUC ROC curve for all signatures in one conditions
  #   Object for drawing ROC curve

  listOfFrame <- list()
  
  conditionName <- ""
  
  title <- ""
  

  sig <- c("PCA3","mixed", "new_lnc")
  
  color <- c("#FF8900","#FFBD4F","#77ABD6")
  
  # Extract data frame for all signature in one condition
  if(oneCondition == "normal-tumor"){
    
    conditionName <- "Normal vs Tumor"
    
    for (i in 1:length(signature)){
      
      listOfFrame[[i]] <- frame[c("condition",signature[[i]])]
      
    }
    
  }else if (oneCondition == "LH"){
    
    conditionName <- "LR vs HR"
    
    for (i in 1:length(signature)){
      
      listOfFrame[[i]] <- extract.condition(frame,"LR","HR",signature[[i]])
      
    }
    
  }else if (oneCondition == "combine"){
    
    conditionName <- "Normal+LR vs HR+IR"
    
    for (i in 1:length(signature)){
      
      listOfFrame[[i]] <- extract.combine(frame,signature[[i]])
      
    }
    
  }else{
    
    conditionName <- paste("Normal vs ", oneCondition)
    
    for (i in 1:length(signature)){
      
      listOfFrame[[i]] <- extract.risk(frame, oneCondition, signature[[i]])
      
    }
      
  }
    
  # Calculate mean and sd for specific condition
  res1 <- takeDataReturnAUCvalues(listOfFrame[[1]], NUM_RUNS)
  res2 <- takeDataReturnAUCvalues(listOfFrame[[2]], NUM_RUNS)
  res3 <- takeDataReturnAUCvalues(listOfFrame[[3]], NUM_RUNS)
  
  core1 <- do.call("join_scores", lapply(res1, "[[", 2))
  core2 <- do.call("join_scores", lapply(res2, "[[", 2))
  core3 <- do.call("join_scores", lapply(res3, "[[", 2))
  
  
  label1 <- do.call("join_labels", lapply(res1, "[[", 3))
  label2 <- do.call("join_labels", lapply(res2, "[[", 3))
  label3 <- do.call("join_labels", lapply(res3, "[[", 3))
  
  auc1 <- do.call("rbind", lapply(res1, "[[", 1))
  auc2 <- do.call("rbind", lapply(res2, "[[", 1))
  auc3 <- do.call("rbind", lapply(res3, "[[", 1))
  
  col1 <- rbind(round(mean(auc1),2),round(sd(auc1),2)) 
  col2 <- rbind(round(mean(auc2),2),round(sd(auc2),2)) 
  col3 <- rbind(round(mean(auc3),2),round(sd(auc3),2)) 
  
  # Takes predicted scores from 3 signatures vs common signature and coverts them to a list
  listScores = join_scores(core1, core2, core3)
  
  # Takes observed labels from 3 signatures vs common signature and converts them to a list.
  listLabels = join_labels(label1, label2, label3)

  # Format data for performance evaluation calculation
  mdat <- mmdata(listScores,listLabels, modnames = sig, dsids = 1:NUM_RUNS, expd_first ="dsids")
  
  # Calculate performance evaluation measures
  mmcurves <- evalmod(mdat, raw_curves = TRUE)
  
  # Convert a curves and points object to a data frame for ggplot2
  mmdf <- fortify(mmcurves, raw_curves = FALSE)
  
  # Extract data frame for drawing average ROC curve
  frameROC <- subset(mmdf, curvetype == "ROC")
  
  # Concatenate mean and sd of all signature in one conditions
  myrow <- cbind(col1, col2, col3)
  
  colnames(myrow) <- c("PCA3", "Signature mixed", "Signature new_lnc")
  
  rownames(myrow) <- c("Mean", "SD")  
 
  return (list(myrow,conditionName, frameROC))
  
}
################
## F) Function to draw smooth ROC curve for signatures in one condition
################
drawROCcurve <- function(frame, conditionName, mode){
  # Draw smooth ROC curve
  #
  # Args:
  #   frame        : Data frame containes all the data for drawing ROC curve.
  #   conditionName: Specific condition
  #   mode         : String value, TCGA or Nanostring  
  #
  # Results:
  #   Display smooth ROC curve
  
  # Curves' color
  color <- c("#FF8900","#FFBD4F","#77ABD6")
  
  color.gene <- c("#FF8900")
  
  # Transfer specificity & sensitivity in scale percentage
  frame$y<-frame$y*100
  
  frame$x<-(1-frame$x)*100
  
  # Remove NA values
  subData<-frame[which(is.na(frame$modname)==F),]
  
  # Data frame with processed values for all probes
  allProbesData<-data.frame()
  
  # Take the specificity & sensitivity for each probe, create a sliding window (shift by 1) for the sensitivity,
  # in which the mean will be computed, and assigned for each original sensitivity value
  for(i in unique(subData$modname)){
    
    # Take data for one probe
    oneProbeData<-subData[which(subData$modname==i),]
    
    # Take n extracted value (=window for rollapply), and make the mean on them
    extracted_values<-round(length(oneProbeData$y)/10)
    
    # For each value of y, make the mean in the window of n extracted values
    result<-rollapply(data=oneProbeData$y, width=extracted_values,FUN=mean)
    
    # The length may differ between the result and the original
    the_rest<-length(oneProbeData$y)-length(result)
    
    if(the_rest!=0){
      
      oneProbeData$y[1:(length(oneProbeData$y)-the_rest)]<-result
      
    }
    
    # Round the percentage, in order to remove them (because they will create a spot on the plot)
    oneProbeData$y<-round(oneProbeData$y)
    
    # Remove
    oneProbeData<-oneProbeData[!duplicated(oneProbeData$y), ]
    
    # Arrange the extremity of the roc curve
    oneProbeData$y[1]<-0
    
    oneProbeData$y[length(oneProbeData$y)]<-100
    
    oneProbeData$x[1]<-100
    
    oneProbeData$x[length(oneProbeData$x)]<-0
    
    oneProbeData<-oneProbeData[order(-oneProbeData$x),]
    
    allProbesData<-rbind(allProbesData,oneProbeData)
    
  }
  
  listOfROC<-allProbesData
  
  whiteBackground<-theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank())
  
  plotROC <- ggplot(listOfROC,aes(x=x,y=y,group=modname,color=modname))+theme_bw()+whiteBackground+
    geom_path(aes(group=modname),size=2)+
    ylab("sensitivity (%)")+
    xlab("specificity (%)")+
    scale_x_reverse(breaks = seq(100,0,-10),limits=c(100,0))+scale_y_continuous(breaks = seq(0,100,10),limits=c(0,100))+geom_abline(intercept=100)+
    theme(axis.text.x=element_text(color="black",size=12,angle=45,hjust=1),axis.text.y = element_text(color="black",size=14))+
    theme(plot.title = element_text(hjust = 0.5,size=16),axis.title=element_text(size=16,face="bold"),legend.text=element_text(size=14),legend.title=element_text(size=16))
  
  if(mode == "TCGA"){
    
    plotROC <- plotROC + scale_color_manual(values = color, name = "Signature")+
      ggtitle(paste("ROC", conditionName, "in TCGA", sep = " "))
    
  }else if(mode == "NanoString"){
    
    plotROC <- plotROC + scale_color_manual(values = color, name = "Signature")+
      ggtitle(paste("ROC", conditionName, "in Nanostring", sep = " "))
  
  }else {
    
    plotROC <- plotROC + scale_color_manual(values = color.gene, name = "Signature")+
      ggtitle(paste("ROC", conditionName, "in TCGA", sep = " "))
    
  }
  
  return (plotROC)
  
}