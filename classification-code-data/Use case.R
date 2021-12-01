# Installing PSSMCOOL and loading it
install.packages("PSSMCOOL")
library(PSSMCOOL)

# setting up working environment and downloading necessary files from GitHub
current_directory <- "D:/Paper submission/Alireza Mohammadi/second revise/R code/"
setwd(current_directory)
# Downloading the required PSSM files
pssm_url <- 'https://github.com/BioCool-Lab/PSSMCOOL/raw/main/classification-code-data/all_needed_pssms90.zip'
download.file(pssm_url, './all_needed_pssm90.zip', method = 'auto', quiet = FALSE)
unzip('all_needed_pssm90.zip', exdir = 'all_needed_pssm90')
PSSM_directory <- 'all_needed_pssm90/all_needed_pssms90/'

# Downloading positive data and loading it to R
url <- "https://raw.githubusercontent.com/BioCool-Lab/PSSMCOOL/main/classification-code-data/positive.csv"
download.file(url, './PositiveData.csv')
positive_data <- read.csv("./PositiveData.csv", header = TRUE)

# Downloading negative data and loading it to R
url <- "https://raw.githubusercontent.com/BioCool-Lab/PSSMCOOL/main/classification-code-data/negative.csv"
download.file(url, './NegativeData.csv')
negative_data <- read.csv("./NegativeData.csv", header = TRUE)


############################### Positive feature extraction ############################################
# Feature extraction
positiveFeatures<- c()
for(i in 1:dim(positive_data)[1]) {
  ff<-FPSSM2(paste0(PSSM_directory, positive_data[i,1],'.fasta.pssm'), 
             paste0(PSSM_directory, positive_data[i,2],'.fasta.pssm'), 20)
  positiveFeatures<-rbind(positiveFeatures, ff)
}

# Adding row names and class
positiveFirstColumn <- c()
for(i in 1:dim(positive_data)[1]) {
  dd <- paste(positive_data[i,1], '-' ,positive_data[i,2])
  positiveFirstColumn <- rbind(positiveFirstColumn, dd)
}

pos_class <- rep("Interaction", dim(positiveFeatures)[1])
positiveFeatures2 <- cbind(positiveFirstColumn, positiveFeatures, pos_class)
############################### Negative feature extraction ############################################
# Feature extraction
negativeFeatures <- c()
for(i in 1:dim(negative_data)[1]) {
  ff2<-FPSSM2(paste0(PSSM_directory, negative_data[i,1],'.fasta.pssm'), 
              paste0(PSSM_directory, negative_data[i,2],'.fasta.pssm'), 20)
  negativeFeatures<-rbind(negativeFeatures, ff2)
}
# Adding row names and class
negativeFirstColumn <- c()
for(i in 1:dim(negative_data)[1]) {
  dd2 <- paste(negative_data[i,1], '-' ,negative_data[i,2])
  negativeFirstColumn <- rbind(negativeFirstColumn, dd2)
}
neg_class <- rep("Non.Interaction", dim(negativeFeatures)[1])
negativeFeatures2 <- cbind(negativeFirstColumn, negativeFeatures, neg_class)

# Merging two feature vectors
mainDataSet <- rbind(positiveFeatures2, negativeFeatures2)


############################### Preparing data set for model training ############################################
# In the following we are going to carry out classification on the data we have prepared so far (mainDataSet)
# First we need to install and load caret package and its dependencies
install.packages('caret', dependencies = TRUE)
library(caret)
bmp.R2.submission.data.df <- as.data.frame(mainDataSet)
colnames(bmp.R2.submission.data.df)[1] <- "interactions"
dim(bmp.R2.submission.data.df)#1730  102
#Assigning the Uniprot IDs for each protein pairs to the row name
rownames(bmp.R2.submission.data.df) <- bmp.R2.submission.data.df$interactions
#Removing the Uniprot IDs
bmp.R2.submission.data.df <-bmp.R2.submission.data.df[,-1]
View(bmp.R2.submission.data.df)
colnames(bmp.R2.submission.data.df) <- c(paste0('Frt', 1: dim(positiveFeatures)[2]), 'Class')
dim(bmp.R2.submission.data.df)#1730  101
table(bmp.R2.submission.data.df$Class)
#Interaction Non-Interaction 
#865             865 
bmp.R2.submission.data.df$Class <- 
  as.factor(bmp.R2.submission.data.df$Class)
write.csv(bmp.R2.submission.data.df, 'DataSet.csv')


############################### Training model with three classifier ############################################

#setting.the.trainControl===========
bmp.R2.submission.data.df <- read.csv("DataSet.csv")
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
  train(Class ~ ., data = bmp.R2.submission.data.df, 
        method = "treebag", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE)
print(cross.validation.bulit.model.treebag$results)
#parameter  Accuracy     Kappa  AccuracySD    KappaSD
#1      none 0.9965351 0.9930707 0.005582867 0.01116413


#10-fold cross-validation using "Single C5.0 Tree (C5.0Tree)" classifier=======
cross.validation.bulit.model.C5.0Tree <- 
  train(Class ~ ., data = bmp.R2.submission.data.df, 
        method = "C5.0Tree", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE)
print(cross.validation.bulit.model.C5.0Tree$results)
#parameter  Accuracy     Kappa  AccuracySD     KappaSD
#1      none 0.9976911 0.9953822 0.004028016 0.008056142


#10-fold cross-validation using "Partial Least Squares (pls)" classifier=======
cross.validation.bulit.model.pls <- 
  train(Class ~ ., data = bmp.R2.submission.data.df, 
        method = "pls", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE)
print(cross.validation.bulit.model.pls$results)
#  ncomp  Accuracy      Kappa AccuracySD    KappaSD
#1     1 0.5005948 0.00231924 0.01330889 0.02666192
#2     2 0.5070671 0.01413974 0.03944775 0.07848283
#3     3 0.5324142 0.06473979 0.02790055 0.05655902


