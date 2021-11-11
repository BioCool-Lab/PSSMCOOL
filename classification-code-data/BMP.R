# install.packages('caret', dependencies = TRUE)
library(caret)
#in this script I'm going to do classification using the data set prepared by Alireza
setwd("F:\\article400\\javad2")
bmp.R2.submission.data.df <- read.csv("DataSet2.csv")
dim(bmp.R2.submission.data.df)#1730  102
View(bmp.R2.submission.data.df)
#Assigning the Uniprot IDs for each protein pairs to the row name
rownames(bmp.R2.submission.data.df) <- 
  bmp.R2.submission.data.df$interactions
#Removing the Uniprot IDs 
bmp.R2.submission.data.df <- 
  bmp.R2.submission.data.df[,-1]
dim(bmp.R2.submission.data.df)#1730  101
table(bmp.R2.submission.data.df$class)
#Interaction Non-Interaction 
#865             865 
bmp.R2.submission.data.df$class <- 
  as.factor(bmp.R2.submission.data.df$class)

#setting.the.trainControl===========
setting.the.trainControl.3 <- function()
{
  #setting the trainControl function parameter: repeated CV; downsampling; 
  set.seed(100)
  fitControl <- trainControl(## 10-fold CV
    method = "cv",
    returnData = TRUE,
    classProbs = TRUE,
  )
  return(fitControl)
  
}
#setting cross validation parameters
trainControl.for.PSSM <- setting.the.trainControl.3()

#10-fold cross-validation using "Bagged CART (treebag)" classifier=======
cross.validation.bulit.model.treebag <- 
  train(class ~ ., data = bmp.R2.submission.data.df, 
      method = "treebag", 
      trControl = trainControl.for.PSSM, 
      verbose = FALSE)

print(cross.validation.bulit.model.treebag$results)
# parameter Accuracy    Kappa  AccuracySD    KappaSD
# 1      none 0.995947 0.991893 0.005486098 0.01097446


#10-fold cross-validation using "Single C5.0 Tree (C5.0Tree)" classifier=======
cross.validation.bulit.model.C5.0Tree <- 
  train(class ~ ., data = bmp.R2.submission.data.df, 
        method = "C5.0Tree", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE)


print(cross.validation.bulit.model.C5.0Tree$results)
# parameter  Accuracy     Kappa  AccuracySD    KappaSD
# 1      none 0.9965351 0.9930693 0.005582827 0.01116793

#10-fold cross-validation using "Partial Least Squares (pls)" classifier=======
cross.validation.bulit.model.pls <- 
  train(class ~ ., data = bmp.R2.submission.data.df, 
        method = "pls", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE)

print(cross.validation.bulit.model.pls$results)
# ncomp  Accuracy       Kappa AccuracySD    KappaSD
# 1     1 0.5034885  0.01032276 0.01831448 0.03438894
# 2     2 0.4861705 -0.02718683 0.05446108 0.10915104
# 3     3 0.5427787  0.08574924 0.04465618 0.08864498





