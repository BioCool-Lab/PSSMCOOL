## ----vali,echo=FALSE,fig.cap="Figure 1: List of implemented features in PSSMCOOL package",out.width = '70%'----
knitr::include_graphics("figures/feature_table.jpg")

## -----------------------------------------------------------------------------
library(PSSMCOOL)

## ----pssm-ac,echo=FALSE,fig.cap="Figure 2: process of extracting PSSM-AC feature vector from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/pssm_ac.jpg")

## -----------------------------------------------------------------------------
 w<-PSSMAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
 head(w, n = 50)

## ----dpc-pssm,echo=FALSE,fig.cap="Figure 3: process of extracting DPC-PSSM feature vector from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/dpc-pssm.jpg")

## -----------------------------------------------------------------------------
 ss<-DPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ss, n = 50)

## ----trigram,echo=FALSE,fig.cap="Figure 4: process of extracting trigram-PSSM feature vector from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/trigram.jpg")

## -----------------------------------------------------------------------------
as<-trigrame_pssm(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GSI6.txt.pssm"))
head(as, n = 50)

## ----pse-pssm,echo=FALSE,fig.cap="Figure 5: process of extracting Pse-PSSM feature vector from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/pse-pssm.jpg")

## -----------------------------------------------------------------------------
 v<-pse_pssm(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(v, n = 50)

## ----k-separated,echo=FALSE,fig.cap="Figure 6: process of extracting K-separated-bigam-PSSM feature vector from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/k-separated.jpg")

## -----------------------------------------------------------------------------
 w<-k_seperated_bigrame(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),5)
head(w, n = 50)

## ----eedp,echo=FALSE,fig.cap="Figure 7: process of extracting EDP-EEDP-MEDP feature vectors from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/EEDP.jpg")

## -----------------------------------------------------------------------------
 as<-EDP_MEDP(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GS61.txt.pssm"))
head(as, n = 50)

## ----ab-pssm,echo=FALSE,fig.cap="Figure 8: process of extracting AB-PSSM feature vectors from PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/AB-PSSM.jpg")

## -----------------------------------------------------------------------------
  zz<- AB_PSSM(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"))
head(zz[1], n = 50)

## -----------------------------------------------------------------------------
 as<-AATP_TPCC(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GQS7.txt.pssm"))
head(as, n = 50)

## -----------------------------------------------------------------------------
 A<-CS_PSe_PSSM(system.file("extdata", "C7GSI6.txt.pssm", package="PSSMCOOL"),"total")
head(A, n = 50)

## ----fpssm,echo=FALSE,fig.cap="Figure 9: process of making FPSSM Matrix and extracting corresponding feature vectors",out.width = '70%'----
knitr::include_graphics("figures/s-fpssm.jpg")

## -----------------------------------------------------------------------------
 q<-FPSSM(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),20)
head(q, n = 50)

## ----scsh2,echo=FALSE,fig.cap="Figure 10: process of extracting scsh2 feature vector",out.width = '70%'----
knitr::include_graphics("figures/SCSH2.jpg")

## ----scshtable,echo=FALSE,fig.cap="Figure 11: tables of all 2-mers and all 3-mers",out.width = '70%'----
knitr::include_graphics("D:/Edit_thesis/thesis_pics/scshtable.jpg")

## -----------------------------------------------------------------------------
 zz<- scsh2(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"),2)
head(zz, n = 200)

## -----------------------------------------------------------------------------
 w<-rpssm(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)

## ----ccpssm,echo=FALSE,fig.cap="Figure 12: process of extracting CC-PSSM feature vector",out.width = '70%'----
knitr::include_graphics("figures/cc-pssm.jpg")

## -----------------------------------------------------------------------------
 aa<-pssm_cc(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),18)
head(aa, n = 50)

## -----------------------------------------------------------------------------
 as<-Discrete_Cosine_Transform(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)

## ----dwt,echo=FALSE,fig.cap="Figure 13: Schematic diagram of a DWT with 4 levels",out.width = '70%'----
knitr::include_graphics("figures/dwt.jpg")

## -----------------------------------------------------------------------------
 as<-dwt_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)

## ----disulfid,echo=FALSE,fig.cap="Figure 14: The process of extracting disulfide-PSSM feature from the PSSM matrix",out.width = '70%'----
knitr::include_graphics("figures/disulfid.jpg")

## -----------------------------------------------------------------------------
  aq<-disulfid(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(aq[,1:50])

## -----------------------------------------------------------------------------
 ss<-DP_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ss, n = 50)

## -----------------------------------------------------------------------------
  as<-DFMCA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7)
head(as, n = 50)

## -----------------------------------------------------------------------------
 as<-grey_pssm_pseAAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)

## ----smooth,echo=FALSE,fig.cap="Figure 15: process of smoothed-pssm matrix generation, (A) represents the PSSM matrix and (B) represents the smoothed-pssm matrix",out.width = '70%'----
knitr::include_graphics("figures/smoothed.jpg")

## -----------------------------------------------------------------------------
 w<-smoothed_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7,11,c(2,3,8,9))
head(w[,1:50], n = 50)

## -----------------------------------------------------------------------------
  w<-kiderafactor(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),c(2,3,8,9))
head(w[,1:50], n = 50)

## -----------------------------------------------------------------------------
 w<-MBMGACPSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)

## -----------------------------------------------------------------------------
 w<-LPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)

## ----pssm400,echo=FALSE,fig.cap="Figure 16: process of extracting PSSM400 feature vector, which for amino acid S, represents the corresponding rows in PSSM Matrix",out.width = '70%'----
knitr::include_graphics("figures/pssm400.jpg")

## -----------------------------------------------------------------------------
 q<-pssm400(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
head(q, n = 50)

## -----------------------------------------------------------------------------
 as<-PSSMBLOCK(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),5)
head(as, n = 50)

## ----pssmsd,echo=FALSE,fig.cap="Figure 17: process of extracting PSSM-SD feature vector values for column j",out.width = '70%'----
knitr::include_graphics("figures/pssmsd.jpg")

## -----------------------------------------------------------------------------
 ww<-PSSM_SD(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ww, n = 50)

## -----------------------------------------------------------------------------
 q<-pssm_seg(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),3)
head(q, n = 50)

## -----------------------------------------------------------------------------
 w<-SOMA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)

## -----------------------------------------------------------------------------
 w<-SVD_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 20)

