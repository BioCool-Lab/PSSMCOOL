setwd("F://article400/javad2")
library(PSSMCOOL)
interactFile <- read.delim('annot__dat (3).txt', skip = 1, sep = "", header = FALSE)
print(dim(interactFile))
print(typeof(interactFile))
interactFile2 <- interactFile[, 1:2]
interactFile2 <- as.data.frame(interactFile2)
interactFile2 <- unique(interactFile2)
interactFile3 <- interactFile2[!grepl('intact', interactFile2$V1), ]
interactFile3 <- interactFile3[!grepl('intact', interactFile3$V2), ]
interactFile3 <- interactFile3[!grepl('PRO', interactFile3$V1), ]
interactFile3 <- interactFile3[!grepl('PRO', interactFile3$V2), ]
interactFile3 <- interactFile3[!grepl('chebi', interactFile3$V1), ]
interactFile3 <- interactFile3[!grepl('chebi', interactFile3$V2), ]
interactFile3 <- interactFile3[!grepl('-', interactFile3$V1), ]
interactFile3 <- interactFile3[!grepl('-', interactFile3$V2), ]
interactFile3 <- interactFile3[substr(interactFile3$V1, 1, 9) == 'uniprotkb', ]
interactFile3
g2<- lapply(interactFile3, function(x) substr(x, 11, nchar(x)))
positiveData <- as.data.frame(g2)

setwd('C:\\Users\\LapTop\\Downloads\\Compressed\\name\\human\\human_pssms90')
allFiles <- list.files()
allFiles2 <- gsub('.{0,11}$', '', allFiles)
s2 <- union(positiveData$V1, positiveData$V2)
length(intersect(s2, allFiles2))
length(s2)
positiveData2 <- positiveData[positiveData$V1 %in% intersect(s2, allFiles2) & positiveData$V2 %in% intersect(s2, allFiles2),]

dim(positiveData2)
s3 <- union(positiveData2$V1, positiveData2$V2)
length(s3)
length(intersect(s3, allFiles2))

write.csv(positiveData2, "F://article400/javad2/positive.csv",row.names = FALSE)


################################# Negative Generation ################################
library(combinat)

neg1 <- union(positiveData2$V1, positiveData2$V2)

negWhole <- combn2(neg1)

negWhole <- as.data.frame(negWhole)
negWhole2 <- negWhole


for(i in 1:dim(positiveData2)[1]) {
  negWhole2 <- negWhole2[!(negWhole2$V1 == positiveData2[i,1] & negWhole2$V2 == positiveData2[i,2]),] 
}

negativeData <- negWhole2[sample(1:dim(positiveData2)[1], dim(positiveData2)[1]),]

ggg<- positiveData2 == negativeData

write.csv(negativeData, "F://article400/javad2/negative.csv",row.names = FALSE)

getwd()
positiveFeatures<- c()
for(i in 1:dim(positiveData2)[1]) {
  ff<-FPSSM2(paste0(positiveData2[i,1],'.fasta.pssm'), paste0(positiveData2[i,2],'.fasta.pssm'), 20)
  positiveFeatures<-rbind(positiveFeatures, ff)
}

colnames(positiveFeatures) <- 1:100
rownames(positiveFeatures) <- 1:dim(positiveFeatures)[1]

positiveFirstColumn <- c()

for(i in 1:dim(positiveData2)[1]) {
  dd <- paste(positiveData2[i,1], '-' ,positiveData2[i,2])
  positiveFirstColumn <- rbind(positiveFirstColumn, dd)
}

class <- rep(1, dim(positiveFeatures)[1])
positiveFeatures2 <- cbind(positiveFirstColumn, class, positiveFeatures)

############################### Negative ############################################

negativeFirstColumn <- c()

for(i in 1:dim(negativeData)[1]) {
  dd2 <- paste(negativeData[i,1], '-' ,negativeData[i,2])
  negativeFirstColumn <- rbind(negativeFirstColumn, dd2)
}

#########################################################################################


negativeFeatures<- c()

for(i in 1:dim(negativeData)[1]) {
  ff2<-FPSSM2(paste0(negativeData[i,1],'.fasta.pssm'), paste0(negativeData[i,2],'.fasta.pssm'), 20)
  negativeFeatures<-rbind(negativeFeatures, ff2)
}

colnames(negativeFeatures) <- 1:100
rownames(negativeFeatures) <- 1:dim(negativeFeatures)[1]
class <- rep(0, dim(positiveFeatures)[1])
negativeFeatures2 <- cbind(negativeFirstColumn,class, negativeFeatures)

mainDataSet <- rbind(positiveFeatures2, negativeFeatures2)

colnames(mainDataSet)[1] <- 'interactions'

rownames(mainDataSet) <- 1:dim(mainDataSet)[1]
View(mainDataSet)

write.csv(mainDataSet, "F://article400/javad2/DataSet.csv",row.names = FALSE)
getwd()
