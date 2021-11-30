setwd("F://article400/javad2")

positiveData <- read.delim("positive.csv", skip = 1, sep = ",", header = FALSE)

negativeData <- read.delim("negative.csv", skip = 1, sep = ",", header = FALSE)

pos1 <- union(positiveData$V1, positiveData$V2)
neg1 <- union(negativeData$V1, negativeData$V2)

allData <- union(as.vector(pos1), as.vector(neg1))

length(allData)

setwd("C:\\Users\\LapTop\\Downloads\\Compressed\\name\\human\\human_pssms90")

library(PSSMCOOL)


#################################### PSSMAC ########################################################
pssmacTimes <- c()
for(i in 1:length(allData)) {
  pssmacTimes <- c(pssmacTimes, as.matrix(system.time({PSSMAC(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### DPC ########################################################
DPCTimes <- c()
for(i in 1:length(allData)) {
  DPCTimes <- c(DPCTimes, as.matrix(system.time({DPC_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### trigram ########################################################
trigramTimes <- c()
for(i in 1:length(allData)) {
  trigramTimes <- c(trigramTimes, as.matrix(system.time({trigrame_pssm(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### psepssm ########################################################
psepssmTimes <- c()
for(i in 1:length(allData)) {
  psepssmTimes <- c(psepssmTimes, as.matrix(system.time({pse_pssm(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### k_separated ########################################################
k_separatedTimes <- c()
for(i in 1:length(allData)) {
  k_separatedTimes <- c(k_separatedTimes, as.matrix(system.time({k_seperated_bigrame(paste0(allData[i], ".fasta.pssm"),5)}))[1])
}

#################################### EDP ########################################################
EDPTimes <- c()
for(i in 1:length(allData)) {
  EDPTimes <- c(EDPTimes, as.matrix(system.time({EDP_MEDP(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### AB-PSSM ########################################################
ABPSSMTimes <- c()
for(i in 1:length(allData)) {
  ABPSSMTimes <- c(ABPSSMTimes, as.matrix(system.time({AB_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### AATP_TPC ########################################################
aatpTimes <- c()
for(i in 1:length(allData)) {
  aatpTimes <- c(aatpTimes, as.matrix(system.time({AATP_TPCC(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### cs_pse ########################################################
cspseTimes <- c()
for(i in 1:length(allData)) {
  cspseTimes <- c(cspseTimes, as.matrix(system.time({CS_PSe_PSSM(paste0(allData[i], ".fasta.pssm"), "total")}))[1])
}

#################################### FPSSM ########################################################
FPSSMTimes <- c()
for(i in 1:length(allData)) {
  FPSSMTimes <- c(FPSSMTimes, as.matrix(system.time({FPSSM(paste0(allData[i], ".fasta.pssm"), 20)}))[1])
}

#################################### SCSH2 ########################################################
SCSHTimes <- c()
for(i in 1:length(allData)) {
  SCSHTimes <- c(SCSHTimes, as.matrix(system.time({scsh2(paste0(allData[i], ".fasta.pssm"),2)}))[1])
}

#################################### rpssm ########################################################
rpssmTimes <- c()
for(i in 1:length(allData)) {
  rpssmTimes <- c(rpssmTimes, as.matrix(system.time({rpssm(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### pssm_cc ########################################################
pssm_ccTimes <- c()
for(i in 1:length(allData)) {
  pssm_ccTimes <- c(pssm_ccTimes, as.matrix(system.time({pssm_cc(paste0(allData[i], ".fasta.pssm"), 18)}))[1])
}

#################################### discreteCosin ########################################################
discreteCosinTimes <- c()
for(i in 1:length(allData)) {
  discreteCosinTimes <- c(discreteCosinTimes, as.matrix(system.time({Discrete_Cosine_Transform(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### dwt_pssm ########################################################
dwtTimes <- c()
for(i in 1:length(allData)) {
  dwtTimes <- c(dwtTimes, as.matrix(system.time({dwt_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### disulfid ########################################################
disulfidTimes <- c()
for(i in 1:length(allData)) {
  disulfidTimes <- c(disulfidTimes, as.matrix(system.time({disulfid(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### dp_pssm ########################################################
dppssmTimes <- c()
for(i in 1:length(allData)) {
  dppssmTimes <- c(dppssmTimes, as.matrix(system.time({DP_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### DFMCA ########################################################
DFMCATimes <- c()
for(i in 1:length(allData)) {
  DFMCATimes <- c(DFMCATimes, as.matrix(system.time({DFMCA_PSSM(paste0(allData[i], ".fasta.pssm"),7)}))[1])
}

#################################### GREYPSSM ########################################################
GREYTimes <- c()
for(i in 1:length(allData)) {
  GREYTimes <- c(GREYTimes, as.matrix(system.time({grey_pssm_pseAAC(paste0(allData[i], ".fasta.pssm"))}))[1])
}


#################################### smoothed ########################################################
smoothedTimes <- c()
for(i in 1:length(allData)) {
  smoothedTimes <- c(smoothedTimes, as.matrix(system.time({smoothed_PSSM(paste0(allData[i], ".fasta.pssm"),7,11,c(2,3,8,9))}))[1])
}


#################################### kiderafactor ########################################################
kideraTimes <- c()
for(i in 1:length(allData)) {
  kideraTimes <- c(kideraTimes, as.matrix(system.time({kiderafactor(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### MBMGAC ########################################################
MBMGACTimes <- c()
for(i in 1:length(allData)) {
  MBMGACTimes <- c(MBMGACTimes, as.matrix(system.time({MBMGACPSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### LPC_PSSM ########################################################
LPCTimes <- c()
for(i in 1:length(allData)) {
  LPCTimes <- c(LPCTimes, as.matrix(system.time({LPC_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### PSSM400 ########################################################
PSSM400Times <- c()
for(i in 1:length(allData)) {
  PSSM400Times <- c(PSSM400Times, as.matrix(system.time({pssm400(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### PSSMBLOCK ########################################################
PSSMBLOCKTimes <- c()
for(i in 1:length(allData)) {
  PSSMBLOCKTimes <- c(PSSMBLOCKTimes, as.matrix(system.time({PSSMBLOCK(paste0(allData[i], ".fasta.pssm"),5)}))[1])
}


#################################### PSSMSD ########################################################
PSSMSDTimes <- c()
for(i in 1:length(allData)) {
  PSSMSDTimes <- c(PSSMSDTimes, as.matrix(system.time({PSSM_SD(paste0(allData[i], ".fasta.pssm"))}))[1])
}


#################################### PSSMSEG ########################################################
PSSMSEGTimes <- c()
for(i in 1:length(allData)) {
  PSSMSEGTimes <- c(PSSMSEGTimes, as.matrix(system.time({pssm_seg(paste0(allData[i], ".fasta.pssm"))}))[1])
}


#################################### somapssm ########################################################
somapssmTimes <- c()
for(i in 1:length(allData)) {
  somapssmTimes <- c(somapssmTimes, as.matrix(system.time({SOMA_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}

#################################### svdpssm ########################################################
svdpssmTimes <- c()
for(i in 1:length(allData)) {
  svdpssmTimes <- c(svdpssmTimes, as.matrix(system.time({SVD_PSSM(paste0(allData[i], ".fasta.pssm"))}))[1])
}


AllFeatureTimes <- cbind(pssmacTimes, DPCTimes, trigramTimes, psepssmTimes, k_separatedTimes, ABPSSMTimes, cspseTimes,
                         EDPTimes, aatpTimes, FPSSMTimes, SCSHTimes, rpssmTimes, pssm_ccTimes, discreteCosinTimes,
                         dwtTimes, disulfidTimes, dppssmTimes, DFMCATimes, GREYTimes, smoothedTimes, kideraTimes,
                         MBMGACTimes, LPCTimes, PSSM400Times, PSSMBLOCKTimes, PSSMSDTimes, PSSMSEGTimes, somapssmTimes,
                         svdpssmTimes)

rownames(AllFeatureTimes) <- allData


write.csv(AllFeatureTimes, "F://article400/javad2/AllFeatureTimes.csv")


#################################### GET PROTEINS LENGSES #######################################

getwd()

proteinLengthes <- c()
for(i in 1:length(allData)) {
  x <- read.delim(paste0(allData[i], ".fasta.pssm"), skip = 2 , sep = "", header = FALSE)
  proteinLengthes <- c(proteinLengthes, dim(x)[1])
}

Protein_Lengthes <- cbind(allData, proteinLengthes)

colnames(Protein_Lengthes) <- c('Protein', 'Length')
View(Protein_Lengthes)

write.csv(Protein_Lengthes, "F://article400/javad2/ProteinLengthes.csv", row.names = FALSE)
