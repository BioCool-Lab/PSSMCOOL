# PSSMCOOL

This package includs all features that has been extracted from Position-Specific Scoring Matrix (PSSM) in POSSUM website: http://possum.erc.monash.edu/ and some of other features
from several articles. each function in this package corresponds to one feature that has extracted from PSSM Matrix.

###### **Install from GitHub :  devtools::install_github("BioCool-Lab/PSSMCOOL")**


## 1 Preface
<br></br>
<font size="4"> Feature extraction or feature encoding is a fundamental step in the construction of high-quality machine learning-based models. Specifically, this is a key step for determining the effectiveness of the trained models in bioinformatics applications. In the last two decades, a variety of feature encoding schemes have been proposed in order to exploit useful patterns from protein sequences. Such schemes are often based on sequence information or physicochemical properties of amino acids. Although direct features derived from sequences themselves (such as amino acid compositions, dipeptide compositions, and counting of k-mers) are regarded as essential for training models, an increasing number of studies have shown that evolutionary information in the form of PSSM profiles is much more informative than sequence information alone. Accordingly, PSSM-based feature descriptors have been commonly used as indispensable primary features to construct models, filling a major gap in the current bioinformatics research. For example, PSSM-based feature descriptors have successfully improved the prediction performance of structural and functional properties of proteins across a wide spectrum of bioinformatics applications. These predictions can be applied for protein fold recognition and the prediction of protein structural classes, protein-protein interactions, protein subcellular localization, RNA-binding sites, and protein functions. But at the same time, there is no comprehensive, simple tool in the R programming language for extracting all of these features from the PSSM and displaying it in the output. PSSMCOOL package here is developed in R for these purposes. First, in Figure 1 which is a table that lists all of the features implemented in this package with their feature-lengths is brought. Then each one of these features will be explained in full detail.</font>

<br></br>

<font size="4"> PSSMCOOL Package is currently available on CRAN website:</font>

[ https://CRAN.R-project.org/package=PSSMCOOL]( https://CRAN.R-project.org/package=PSSMCOOL)

<br></br>

<font size="4"> for issues about this package :</font>

[https://github.com/BioCool-Lab/PSSMCOOL/issues](https://github.com/BioCool-Lab/PSSMCOOL/issues)

<br></br>
![](vignettes/figures/feature_table.jpg)
<pre>                  Figure 1: List of implemented features in PSSMCOOL package </pre>

<br>*feature vector length depends on the choice of parameter</br>  **these features produce Matrix of features which its dimension depends on choice of parameter

<br></br>

```{r}
library(PSSMCOOL)
```

## 2 PSSM-AC
<font size="4">This feature, which stands for auto-covariance transformation, for column j, calculates the average of this column as shown in Figure 2, Then subtracts the resulting number from the elements on the rows i and i + g of this column, and finally multiplies them and calculates the sum by changing the variable i from 1 to L-g. Because the variable j changes between 1 and 20 and the variable g changes between 1 and 10, eventually a feature vector of length 200 will be obtained.</font>

<br></br>
![](vignettes/figures/pssm_ac.jpg)
<br></br>
<pre>                   Figure 2: process of extracting PSSM-AC feature vector from PSSM Matrix </pre>

![](vignettes/figures/screens/pssm_ac.JPG)

## 3 DPC-PSSM
<font size="4">This feature stands for dipeptide composition, which multiplies the values that are located in two consecutive rows and two different columns. Having calculated these values for different rows and columns, they are summed. Next, for both columns, the sum is divided by L-1. Since the result depends on two different columns, eventually a feature vector of length 400, according to Figure 3 and following equation ,will be obtained.</font>
<br></br>

![](vignettes/figures/dpc-pssm.jpg)
<pre>                   Figure 3: process of extracting DPC-PSSM feature vector from PSSM Matrix </pre>
<br></br>

![](vignettes/figures/screens/dpc_pssm.JPG)


## 4 Trigram-PSSM
<font size="4">This feature vector is of length 8000, which is extracted from the PSSM. If we multiply elements available in three consecutive rows and three different columns of the PSSM by each other, and apply this to all rows (all three consecutive rows) and then sum these numbers, eventually one of the elements of feature vector with length 8000 corresponding to the three selected columns will be obtained. Because we have 20 different columns, The final feature vector will be of length 8000 = 20 * 20 * 20. Figure 4: shows these steps. For example, in this figure for three marked rows and columns, the numbers obtained from the intersection of these rows and columns marked with a blue dotted circle around them, are multiplied to each other. </font>

![](vignettes/figures/trigram.jpg)

<pre>                   Figure 4: process of extracting trigram-PSSM feature vector from PSSM </pre>
<br></br>

![](vignettes/figures/screens/trigram.JPG)
<br></br>

## 5 Pse-PSSM
<font size="4"> The length of this feature vector is 320. The first 20 numbers of this feature vector are the mean of 20 columns in PSSM, and the next values for each column are the mean squares of the difference between the elements of row i and i + lag in this column. Because the lag value varies between 1 and 15, the final feature vector will have a length of 320. Figure 5: and following equation shows the process of this function and the corresponding mathematical equation, respectively.</font>

![](vignettes/figures/pse-pssm.jpg)
<pre>                  Figure 5: process of extracting Pse-PSSM feature vector from PSSM </pre>
<br></br>

![](vignettes/figures/screens/pse_pssm.JPG)

## 6 k-Separated-bigram-PSSM
<font size="4"> This feature is almost identical to the DPC feature, and in fact, the DPC feature is part of this feature (for k = 1) and for two different columns, it considers rows that have distance k.</font>

![](vignettes/figures/k-separated.jpg)
<pre>                  Figure 6: process of extracting K-separated-bigam-PSSM feature vector from PSSM </pre>

![](vignettes/figures/screens/k-separated.JPG)


## 7 EDP-EEDP-MEDP
<font size="4"> In this group of features, in order to use uniform dimensions to show proteins of different lengths, in the first step, the average evolutionary score between adjacent residues is calculated using the following equations:</font>

![](vignettes/figures/screens/edp.JPG)
![](vignettes/figures/EEDP.jpg)
<pre>                  Figure 7: process of extracting EDP-EEDP-MEDP feature vectors from PSSM </pre>
<br></br>

![](vignettes/figures/screens/edp2.JPG)

## 8 AB-PSSM
<font size="4"> This feature consists of two types of feature vectors. At first, each protein sequence is divided into 20 equal parts, each of which is called a block, and in each block, the row vectors of the PSSM
related to that block are added together. The resulting final vector is divided by the length of that block, which is equal to 5% of protein length. Finally, by placing these 20 vectors side by side, the first feature vector of length 400 is obtained. The second feature for each amino acid in each column is the average of the positive numbers in that column and for each block, and these 20 values, corresponding to 20 blocks, are placed next to each other, and therefore for each of the 20 types of amino acids, a vector of length 20 is obtained, and by placing these together,the second feature vector of length 400, is obtained. Figure 8 represents this process.</font>

![](vignettes/figures/AB-PSSM.jpg)
<pre>                  Figure 8: process of extracting AB-PSSM feature vectors from PSSM </pre>

![](vignettes/figures/screens/ab-pssm.JPG)

## 9 AATP-TPC
<font size="4"> In this feature, at first, a TPM matrix is constructed from the PSSM, which has represented by a vector corresponding to the following equation:</font>
![](vignettes/figures/screens/aatp1.JPG)

<font size="4"> In the above equation, the numerator is the same as the equation related to DPC-PSSM feature without considering its coefficient. By placing these components together, a TPC feature vector of length 400 is obtained, and if we add the AAC feature vector of length 20 which is the average of columns of the PSSM to the beginning of this vector, AATP feature vector of length 420 is obtained.</font>

![](vignettes/figures/screens/aatp2.JPG)

![](vignettes/figures/screens/cs-pse1.JPG)
![](vignettes/figures/screens/cs-pse2.JPG)
![](vignettes/figures/screens/cs-pse3.JPG)

![](vignettes/figures/screens/df-pssm1.JPG)
![](vignettes/figures/screens/df-pssm2.JPG)

## 12 SCSH2

<font size="4"> To generate this feature vector, the consensus sequence corresponding to the protein sequence is extracted using the PSSM. Then, by placing these two sequences next to each other, a matrix with dimensions of 2 * L will be created. In the next step, each component in the upper row of this matrix is connected to two components in the lower row of this matrix, and thus a graph similar to a bipartite graph could be created.
Now in this graph, each path of length 2 specifies a 3-mer and each path of length 1 denotes a 2-mer corresponding to these two sequences. Now if we consider a table consisting of two rows and 8000 columns so that the first row contains all possible 3-mers of 20 amino acids, then for every 3-mer obtained from this graph, we put number 1 below the corresponding cell With that 3-mer in the aforementioned table and 0 in other cells. so This gives us a vector of length 8000. For the 2-mers obtained from this graph, a vector of length 400 is obtained in a similar way. figures 10, 11 show these processes. </font>
<br></br>

![](vignettes/figures/screens/scsh1.JPG)
![](vignettes/figures/screens/scsh2.JPG)

![](vignettes/figures/screens/rpssm1.JPG)
![](vignettes/figures/screens/rpssm2.JPG)

![](vignettes/figures/screens/cc-pssm1.JPG)
![](vignettes/figures/screens/cc-pssm2.JPG)

![](vignettes/figures/screens/discrete-cos.JPG)

![](vignettes/figures/screens/discrete-wave.JPG)
![](vignettes/figures/screens/discrete-wave2.JPG)

## 17 Disulfide_PSSM
<font size="4">For the purpose of predicting disulfide bond in protein at first, the total number of cysteine amino
acids in the protein sequence is counted and their position in the protein sequence is identified.
Then, using a sliding window with a length of 13, moved on the PSSM from top to bottom so
that the middle of the window is on the amino acid cysteine, then the rows below the matrix obtained
from the PSSM with the dimension of 13 x 20 are placed next to each other to get a feature vector
with a length of 260 = 20 * 13 per cysteine. If the position of the first and last cysteine in the
protein sequence is such that the middle of sliding window is not on cysteine residue while moving
on PSSM, then the required number of zero rows from top and bottom is added to the PSSM
matrix to achieve this goal.Thus, for every cysteine amino acid presented in protein sequence, a
feature vector with a length of 260 is formed.Then all the pairwise combinations of these cysteines
is wrote in the first column of a table. In front of each of these pairwise combinations, the
corresponding feature vectors are stuck together to get a feature vector of length 520 for each of
these compounds.Finally, the table obtained in this way will have the number of rows equal to the
number of all pairwise combinations of these cysteines and the number of columns will be equal to
521 (the first column includes the name of these pair combinations). It is easy to divide this
table into training and testing data and predict the desired disulfide bonds between cysteines.Figure 14 shows a schematic of this process:</font>

![](vignettes/figures/screens/disulfid1.JPG)
![](vignettes/figures/screens/disulfid2.JPG)

## 18 DP-PSSM
<p>The extraction of this feature is obtained using the following equations from the PSSM:</p>

![\begin{aligned}
P_{DP-PSSM}^{\alpha}&=[T',G']=[p_1,p_2,...,p_{40+40\times{\alpha}}]\\
T'&=[\bar{T}_1^P,\bar{T}_1^N,\bar{T}_2^P,\bar{T}_2^N,...,\bar{T}_{20}^P,\bar{T}_{20}^N]
\end{aligned}](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cbegin%7Baligned%7D%20P_%7BDP-PSSM%7D%5E%7B%5Calpha%7D%26%3D%5BT%27%2CG%27%5D%3D%5Bp_1%2Cp_2%2C...%2Cp_%7B40&plus;40%5Ctimes%7B%5Calpha%7D%7D%5D%5C%5C%20T%27%26%3D%5B%5Cbar%7BT%7D_1%5EP%2C%5Cbar%7BT%7D_1%5EN%2C%5Cbar%7BT%7D_2%5EP%2C%5Cbar%7BT%7D_2%5EN%2C...%2C%5Cbar%7BT%7D_%7B20%7D%5EP%2C%5Cbar%7BT%7D_%7B20%7D%5EN%5D%20%5Cend%7Baligned%7D)
<br></br>
![\left\{\begin{array}{ll}\bar{T}_j^P=\frac{1}{NP_j}\sum T_{i,j} &  ,if \ T_{i,j}\geq 0\\ 
\bar{T}_j^N=\frac{1}{NN_j}\sum T_{i,j} &  ,if \  T_{i,j}< 0
\end{array}\right.](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bll%7D%5Cbar%7BT%7D_j%5EP%3D%5Cfrac%7B1%7D%7BNP_j%7D%5Csum%20T_%7Bi%2Cj%7D%20%26%20%2Cif%20%5C%20T_%7Bi%2Cj%7D%5Cgeq%200%5C%5C%20%5Cbar%7BT%7D_j%5EN%3D%5Cfrac%7B1%7D%7BNN_j%7D%5Csum%20T_%7Bi%2Cj%7D%20%26%20%2Cif%20%5C%20T_%7Bi%2Cj%7D%3C%200%20%5Cend%7Barray%7D%5Cright.)
<br></br>
![\begin{aligned}
G'&=[G_1,G_2,...,G_{20}]\\
G_j&=[\bar{\Delta}_{1,j}^P,\bar{\Delta}_{1,j}^N,\bar{\Delta}_{2,j}^P,\bar{\Delta}_{2,j}^N,...,\bar{\Delta}_{\alpha,j}^P,\bar{\Delta}_{\alpha,j}^N]
\end{aligned}](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cbegin%7Baligned%7D%20G%27%26%3D%5BG_1%2CG_2%2C...%2CG_%7B20%7D%5D%5C%5C%20G_j%26%3D%5B%5Cbar%7B%5CDelta%7D_%7B1%2Cj%7D%5EP%2C%5Cbar%7B%5CDelta%7D_%7B1%2Cj%7D%5EN%2C%5Cbar%7B%5CDelta%7D_%7B2%2Cj%7D%5EP%2C%5Cbar%7B%5CDelta%7D_%7B2%2Cj%7D%5EN%2C...%2C%5Cbar%7B%5CDelta%7D_%7B%5Calpha%2Cj%7D%5EP%2C%5Cbar%7B%5CDelta%7D_%7B%5Calpha%2Cj%7D%5EN%5D%20%5Cend%7Baligned%7D)
<br></br>
![\left\{\begin{array}{ll}\bar{\Delta}_{k,j}^P=\frac{1}{NDP_j}\sum [T_{i,j}-T_{i+k,j}]^2 &  ,if \ T_{i,j}-T_{i+k,j}\geq 0\\ 
\bar{\Delta}_{k,j}^N=\frac{-1}{NDN_j}\sum [T_{i,j}-T_{i+k,j}]^2 &  ,if \  T_{i,j}-T_{i+k,j}< 0
\end{array}\right.\\
0<k\leq{\alpha}](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bll%7D%5Cbar%7B%5CDelta%7D_%7Bk%2Cj%7D%5EP%3D%5Cfrac%7B1%7D%7BNDP_j%7D%5Csum%20%5BT_%7Bi%2Cj%7D-T_%7Bi&plus;k%2Cj%7D%5D%5E2%20%26%20%2Cif%20%5C%20T_%7Bi%2Cj%7D-T_%7Bi&plus;k%2Cj%7D%5Cgeq%200%5C%5C%20%5Cbar%7B%5CDelta%7D_%7Bk%2Cj%7D%5EN%3D%5Cfrac%7B-1%7D%7BNDN_j%7D%5Csum%20%5BT_%7Bi%2Cj%7D-T_%7Bi&plus;k%2Cj%7D%5D%5E2%20%26%20%2Cif%20%5C%20T_%7Bi%2Cj%7D-T_%7Bi&plus;k%2Cj%7D%3C%200%20%5Cend%7Barray%7D%5Cright.%5C%5C%200%3Ck%5Cleq%7B%5Calpha%7D)
<br></br>
<p>In the above equations, <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT_%7Bi%2Cj%7D%7D"> represents the value on the i-th row and j-th column of the normalized PSSM, which is denoted by <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BM_T%7D">. This matrix is constructed from the PSSM using the following equations:</p>
<br></br>

![mean_i=\frac{1}{20}\sum_{i=1}^{20}E_{i,k}\\
STD_i=\sqrt{\frac{\sum_{u=1}^{20}[E_{i,u}-mean_i]^2}{20}}\\
T_{i,j}=\frac{E_{i,j}-mean_i}{STD_i}](https://latex.codecogs.com/svg.latex?%5Clarge%20mean_i%3D%5Cfrac%7B1%7D%7B20%7D%5Csum_%7Bi%3D1%7D%5E%7B20%7DE_%7Bi%2Ck%7D%5C%5C%20STD_i%3D%5Csqrt%7B%5Cfrac%7B%5Csum_%7Bu%3D1%7D%5E%7B20%7D%5BE_%7Bi%2Cu%7D-mean_i%5D%5E2%7D%7B20%7D%7D%5C%5C%20T_%7Bi%2Cj%7D%3D%5Cfrac%7BE_%7Bi%2Cj%7D-mean_i%7D%7BSTD_i%7D)

<p>In above equations <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5Coverline%7BT%7D_j%5EP%7D"> represents the mean of the positive values of <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5C%7BT_%7Bi%2Cj%7D%7Ci%3D1%2C2%2C...%2CL%5C%7D%7D"> and <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5Coverline%20T_j%5EN%7D"> represents the mean of the negative values of the above set, which in fact the above set represents the j-th column of the matrix <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BM_T%7D"> the expression <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BNP_j%7D"> indicates the number of positive values of set <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5C%7BT_%7Bi%2Cj%7D%7Ci%3D1%2C2%2C...%2CL%5C%7D%7D"> and <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BNN_j%7D"> is related to the number of negative values of the mentioned set. It is clear that this feature vector arises from the connection of two vectors <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BG%27%7D">, <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT%27%7D">. According to the equations, it is clear that the length of the first feature vector is 40 and the length of the second feature vector is <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5Calpha%5Ctimes%2040%7D">, which by selecting 2 in the used article, a feature vector of length 120 is created from the PSSM.</p>

#### Usage of this feature in PSSMCOOL package:
```
 ss<-DP_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ss, n = 50)
```
```
##  [1]  1.0739 -0.4023  1.1132 -0.5915  0.9459 -0.5134  1.3519 -0.6148  2.1035
## [10] -0.7794  1.1112 -0.5070  0.9216 -0.5843  1.3479 -0.5950  1.2737 -0.5028
## [19]  1.1349 -0.6704  1.3017 -0.6891  1.1585 -0.5086  0.9554 -0.5457  1.4547
## [28] -0.7166  1.5252 -0.6998  1.0084 -0.3219  0.6664 -0.4227  1.3648 -0.7855
## [37]  1.0723 -0.6254  1.1957 -0.6171  1.7783 -1.7591  1.5965 -1.6100  1.1498
## [46] -1.5052  1.4737 -1.5524  1.5066 -1.6937
```

## 19 DFMCA-PSSM
<p>In this feature, each of the columns of the PSSM is considered as a non-static time series. Assuming that <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5C%7Bx_i%5C%7D">, <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5C%7By_i%5C%7D"> for i=1,2,...,L represent two different columns of the PSSM, then two cumulative time series X,Y of these two columns are obtained according to the following equations:</p>

![\left\{\begin{array}{ll}X_k=\sum_{i=1}^K x_i & k=1,2,...,L \\ 
Y_k=\sum_{i=1}^k y_i & k=1,2,...,L 
\end{array}\right.](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bll%7DX_k%3D%5Csum_%7Bi%3D1%7D%5EK%20x_i%20%26%20k%3D1%2C2%2C...%2CL%20%5C%5C%20Y_k%3D%5Csum_%7Bi%3D1%7D%5Ek%20y_i%20%26%20k%3D1%2C2%2C...%2CL%20%5Cend%7Barray%7D%5Cright.)

<p>Now, using these two series, two backward moving average are obtained according to the following equations:</p>

![\left\{\begin{array}{ll}\tilde{X}_{k,s}=\frac{1}{s}\sum_{i=-(s-1)}^0 X_{(k-i)} \\ 
\tilde{Y}_{k,s}=\frac{1}{s}\sum_{i=-(s-1)}^0 Y_{(k-i)} 
\end{array}\right.\\
1<s\leq L](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bll%7D%5Ctilde%7BX%7D_%7Bk%2Cs%7D%3D%5Cfrac%7B1%7D%7Bs%7D%5Csum_%7Bi%3D-%28s-1%29%7D%5E0%20X_%7B%28k-i%29%7D%20%5C%5C%20%5Ctilde%7BY%7D_%7Bk%2Cs%7D%3D%5Cfrac%7B1%7D%7Bs%7D%5Csum_%7Bi%3D-%28s-1%29%7D%5E0%20Y_%7B%28k-i%29%7D%20%5Cend%7Barray%7D%5Cright.%5C%5C%201%3Cs%5Cleq%20L)

<p>Finally, each element of the DFMCA feature vector is obtained using the above equations and the following formula:</p>

![f_{DFMCA}^2(s)=\frac{1}{L-s+1}\sum_{k=1}^{L-s+1}(X_k-\tilde{X}_{k,s})(Y_k-\tilde{Y}_{k,s})](https://latex.codecogs.com/svg.latex?%5Clarge%20f_%7BDFMCA%7D%5E2%28s%29%3D%5Cfrac%7B1%7D%7BL-s&plus;1%7D%5Csum_%7Bk%3D1%7D%5E%7BL-s&plus;1%7D%28X_k-%5Ctilde%7BX%7D_%7Bk%2Cs%7D%29%28Y_k-%5Ctilde%7BY%7D_%7Bk%2Cs%7D%29)

<p>According to the above equation, it is clear that each element of this feature vector is obtained by using two different columns of the PSSM, and since we have 20 different columns and the order of the columns does not matter, the length of the obtained feature vector will be equal to <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cbinom%7B20%7D%7B2%7D%3D%5Cfrac%7B20%5Ctimes%2019%7D%7B2%7D%3D190"> </p>

#### Usage of this feature in PSSMCOOL package:
```
 as<-DFMCA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7)
head(as, n = 50)
```
```
##  [1] 1.5671 1.0027 1.3320 1.1677 0.5471 2.3738 1.3395 1.1099 1.2298 1.0010
## [11] 1.1077 1.0803 0.9034 1.5549 0.8586 1.5949 0.8651 0.4748 0.5141 1.0397
## [21] 1.0206 1.1935 1.0423 0.7185 1.5056 1.1396 1.0414 1.0555 0.9881 0.9916
## [31] 1.0804 0.9890 1.1909 0.9159 1.4816 1.0687 0.6386 0.7624 1.0960 1.0048
## [41] 0.8881 0.4023 1.2723 1.0321 0.7996 0.8973 0.5799 0.5536 0.9635 0.6525
```

## 20 grey_pssm_pseAAC

<p>This function produces a feature vector of length 100. The first 20 components of this vector are the same as the normalized frequency of 20 standard amino acids in the protein. The second 20 components of this vector are the average of the 20 columns of the PSSM corresponding to the protein, and the grey system model method is used to define the next 60 components. If we show this feature vector with length 100 as follows:</p>

![V=(\psi_1,\psi_2,...,\psi_{100})](https://latex.codecogs.com/svg.latex?%5Clarge%20V%3D%28%5Cpsi_1%2C%5Cpsi_2%2C...%2C%5Cpsi_%7B100%7D%29)

<p> then the first 20 components are as follows:</p>

![\psi_i=f_i \quad (i=1,2,...,20)](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cpsi_i%3Df_i%20%5Cquad%20%28i%3D1%2C2%2C...%2C20%29)

<p>Where <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7Bf_i%7D"> is the normalized frequency of type i amino acids of the 20 standard amino acids in the protein chain. If we denote the entries of the PSSM by <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7Bp_%7Bi%2Cj%7D%7D">, then the next 20 components of this feature vector are obtained according to the following equation:</p>

![\psi_{j+20}=\alpha_j \quad (j=1,2,...,20),\\
\alpha_j=\frac{1}{L}\sum_{i=1}^L p_{i,j}](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cpsi_%7Bj&plus;20%7D%3D%5Calpha_j%20%5Cquad%20%28j%3D1%2C2%2C...%2C20%29%2C%5C%5C%20%5Calpha_j%3D%5Cfrac%7B1%7D%7BL%7D%5Csum_%7Bi%3D1%7D%5EL%20p_%7Bi%2Cj%7D)

<p> the next 60 components are obtained by following equations:</p>

![\psi_{j+40}=\delta_j \quad (j=1,2,...,60)](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cpsi_%7Bj&plus;40%7D%3D%5Cdelta_j%20%5Cquad%20%28j%3D1%2C2%2C...%2C60%29)

<p> in this equation <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5Cdelta_j%7D">'s are obtained as follows:</p>

![\left\{\begin{array}{ll} \delta_{3j-2}=f_ja_1^j\\
\delta_{3j-1}=f_ja_2^j \quad (j=1,2,...,20) \\
\delta_{3j}=f_jb^j
\end{array}\right.](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bll%7D%20%5Cdelta_%7B3j-2%7D%3Df_ja_1%5Ej%5C%5C%20%5Cdelta_%7B3j-1%7D%3Df_ja_2%5Ej%20%5Cquad%20%28j%3D1%2C2%2C...%2C20%29%20%5C%5C%20%5Cdelta_%7B3j%7D%3Df_jb%5Ej%20%5Cend%7Barray%7D%5Cright.)

<p> where <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7Bf_j%7D">'s are as in above and <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7Ba_1%5Ej%7D">, <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7Ba_2%5Ej%7D">, <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7Bb%5Ej%7D"> are obtained as follows:</p>

![\begin{bmatrix}a_1^j\\a_2^j\\b^j\end{bmatrix}=(B_j^TB_j)^{-1}B_j^TU_j \quad (j=1,2,...,20)](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cbegin%7Bbmatrix%7Da_1%5Ej%5C%5Ca_2%5Ej%5C%5Cb%5Ej%5Cend%7Bbmatrix%7D%3D%28B_j%5ETB_j%29%5E%7B-1%7DB_j%5ETU_j%20%5Cquad%20%28j%3D1%2C2%2C...%2C20%29)

<p>where: </p>

![B_j=\begin{bmatrix}
-p_{2,j}&-(p_{1,j}+0.5p_{2,j})&1\\
-p_{3,j}&-(\sum_{i=1}^2 p_{i,j}+0.5p_{3,j})&1\\
\vdots&\vdots&\vdots\\
-p_{L,j}&-(\sum_{i=1}^{L-1} p_{i,j}+0.5p_{L,j})&1
\end{bmatrix}](https://latex.codecogs.com/svg.latex?%5Clarge%20B_j%3D%5Cbegin%7Bbmatrix%7D%20-p_%7B2%2Cj%7D%26-%28p_%7B1%2Cj%7D&plus;0.5p_%7B2%2Cj%7D%29%261%5C%5C%20-p_%7B3%2Cj%7D%26-%28%5Csum_%7Bi%3D1%7D%5E2%20p_%7Bi%2Cj%7D&plus;0.5p_%7B3%2Cj%7D%29%261%5C%5C%20%5Cvdots%26%5Cvdots%26%5Cvdots%5C%5C%20-p_%7BL%2Cj%7D%26-%28%5Csum_%7Bi%3D1%7D%5E%7BL-1%7D%20p_%7Bi%2Cj%7D&plus;0.5p_%7BL%2Cj%7D%29%261%20%5Cend%7Bbmatrix%7D)

<p> and:</p>

![U_j=\begin{bmatrix}
p_{2,j}-p_{1,j}\\
p_{3,j}-p_{2,j}\\
\vdots\\
p_{L,j}-p_{L-1,j}
\end{bmatrix}](https://latex.codecogs.com/svg.latex?%5Clarge%20U_j%3D%5Cbegin%7Bbmatrix%7D%20p_%7B2%2Cj%7D-p_%7B1%2Cj%7D%5C%5C%20p_%7B3%2Cj%7D-p_%7B2%2Cj%7D%5C%5C%20%5Cvdots%5C%5C%20p_%7BL%2Cj%7D-p_%7BL-1%2Cj%7D%20%5Cend%7Bbmatrix%7D)

#### Usage of this feature in PSSMCOOL package:
```
 as<-grey_pssm_pseAAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)
```
```
##  [1]  0.0534  0.0229  0.0611  0.0687  0.0611  0.1069  0.0153  0.0534  0.0305
## [10]  0.0382  0.0534  0.0840  0.0076  0.0763  0.0687  0.0611  0.0305  0.0153
## [19]  0.0305  0.0611  0.3773  0.2965  0.3371  0.3049  0.1816  0.4356  0.3395
## [28]  0.2866  0.3109  0.2576  0.2622  0.3150  0.2776  0.3187  0.2383  0.3846
## [37]  0.2755  0.1791  0.2180  0.2808 -0.0485 -0.0001 -0.0207 -0.0220  0.0000
## [46] -0.0051 -0.0432  0.0000 -0.0131 -0.0618
```

## 21 Smoothed_pssm 

<p>This feature has been used to predict RNA binding sites in proteins, and therefore a specific feature vector is generated for each residue. To generate this feature vector, a matrix called smoothed_pssm is first created from the PSSM using a parameter ws called the smooth window size, which usually has a default value of 7 and can Change between 3, 5, 7, 9 and 11. The i-th row in the smoothed_PSSM is obtained by summing ws row vectors around the i-th row, which of course corresponds to the i-th row in the protein. If we show this problem with a mathematical equation, then we will have:</p>

![V_{smoothed_i}=V_{i-\frac{(ws-1)}{2}}+...+V_i+...+V_{i+\frac{(ws-1)}{2}}](https://latex.codecogs.com/svg.latex?%5Clarge%20V_%7Bsmoothed_i%7D%3DV_%7Bi-%5Cfrac%7B%28ws-1%29%7D%7B2%7D%7D&plus;...&plus;V_i&plus;...&plus;V_%7Bi&plus;%5Cfrac%7B%28ws-1%29%7D%7B2%7D%7D)

<p>In this regard, <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BV_i%7D">'s represent the row vectors of the PSSM. To obtain the first and last rows of the smoothed_PSSM corresponding to the N-terminal and C-terminal of the protein, zero vectors are added to the beginning and end of the PSSM. The Figure 15 represents these processes schematically.</p>
![](vignettes/figures/screens/smoothed.JPG)

<p>Now, using another parameter w called the slider window size, which its default value is 11 and can change in the interval (3,41)(ste=2), for the resid <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%5Calpha_i%7D">, the feature vector obtained from the smoothed-PSSM will be in the form as follows:</p>

![(V_{smoothed_{i-\frac{(w-1)}{2}}},...,V_{smoothed_i},...,V_{smoothed_{i+\frac{(w-1)}{2}}})](https://latex.codecogs.com/svg.latex?%5Clarge%20%28V_%7Bsmoothed_%7Bi-%5Cfrac%7B%28w-1%29%7D%7B2%7D%7D%7D%2C...%2CV_%7Bsmoothed_i%7D%2C...%2CV_%7Bsmoothed_%7Bi&plus;%5Cfrac%7B%28w-1%29%7D%7B2%7D%7D%7D%29)

<p>Here, as in the previous case, if the residue in question is the first or last residue of the protein, number of <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7B%28%5Cfrac%7Bw-1%7D%7B2%7D%29%7D"> zero vectors are added to the beginning or end of the smoothed matrix to obtain the feature vector. Therefore, the parameter w will determine the length of the feature vector per residue. If its value is 11, the length of the obtained feature vector will be equal to 220 = 11 * 20. eventually, the feature vector values are normalized between -1 and 1.</p>

#### Usage of this feature in PSSMCOOL package:

```
w<-smoothed_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7,11,c(2,3,8,9))
head(w[,1:50], n = 50)
```
```
##        1      2      3      4      5      6      7      8      9     10     11
## 2 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 3 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 8 1.8638 2.8079 0.7699 0.3285 1.0823 2.6192 0.7699 0.5128 2.2446 3.1497 4.3768
## 9 2.5949 3.0769 0.8891 0.3759 1.3512 2.8881 0.8891 0.7817 2.5135 3.6497 4.8768
##       12     13     14     15     16     17     18     19    20     21     22
## 2 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.000 0.0000 0.0000
## 3 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.000 0.0000 0.0000
## 8 2.1762 3.5374 2.4894 0.5613 2.9282 2.0073 2.7474 3.0000 1.926 2.5949 3.0769
## 9 2.2954 3.8064 3.2204 1.0613 3.8808 2.7384 3.2474 3.1192 2.195 2.8164 3.0769
##       23     24     25     26     27     28     29     30     31     32     33
## 2 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 3 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 8 0.8891 0.3759 1.3512 2.8881 0.8891 0.7817 2.5135 3.6497 4.8768 2.2954 3.8064
## 9 0.8891 0.3759 1.3512 2.8164 0.8891 0.7704 2.4841 4.1023 4.8768 2.2659 3.7591
##       34     35     36     37     38     39     40     41     42     43     44
## 2 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 3 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 8 3.2204 1.0613 3.8808 2.7384 3.2474 3.1192 2.1950 2.8164 3.0769 0.8891 0.3759
## 9 3.2204 1.0907 3.9526 2.8881 3.9311 3.2689 2.8068 2.8164 2.8458 1.3512 0.3041
##       45     46     47     48     49     50
## 2 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 3 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
## 8 1.3512 2.8164 0.8891 0.7704 2.4841 4.1023
## 9 1.2795 2.3638 0.6676 0.6987 1.6507 4.4831
```

## 22 Kiderafactor
<font size="30">For producing this feature vector similar to the smoothed-PSSM feature, firstly PSSM is smoothed by appending zero vectors to its head and tail and a sliding window with odd size is utilized. Then this smoothed PSSM is condensed by the Kidera factors to produce feature vector for each residue.</font>
#### Usage of this feature in PSSMCOOL package:
```
 w<-kiderafactor(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),c(2,3,8,9))
head(w[,1:50], n = 50)
```
```
##        1    2     3      4    5      6     7      8      9     10     11    12
## 2  0.000 0.00 0.000  0.000 0.00  0.000  0.00  0.000  0.000  0.000  0.000 0.000
## 3  0.000 0.00 0.000  0.000 0.00  0.000  0.00  0.000  0.000  0.000  0.000 0.000
## 8  0.000 0.00 0.000  0.000 0.00  0.000  0.00  0.000  0.000  0.000 -1.072 0.650
## 9 -1.072 0.65 0.667 -0.199 0.44 -0.484 -0.03 -0.575 -0.617 -0.217 -1.120 0.953
##      13     14    15     16     17     18     19     20     21    22    23
## 2 0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000 0.000 0.000
## 3 0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000 0.000 0.000
## 8 0.667 -0.199 0.440 -0.484 -0.030 -0.575 -0.617 -0.217 -1.120 0.953 0.924
## 9 0.924 -0.099 0.508 -0.433  0.103 -0.637 -0.305  0.152 -1.298 1.250 1.203
##       24    25     26     27     28     29     30     31   32    33     34
## 2  0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000 0.00 0.000  0.000
## 3  0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000 0.00 0.000  0.000
## 8 -0.099 0.508 -0.433  0.103 -0.637 -0.305  0.152 -1.298 1.25 1.203 -0.792
## 9 -0.792 0.346 -0.760 -0.068 -0.546 -0.492 -0.017 -1.334 1.06 1.274 -1.052
##      35     36     37     38     39     40     41    42    43     44    45
## 2 0.000  0.000  0.000  0.000  0.000  0.000  0.000 0.000 0.000  0.000 0.000
## 3 0.000  0.000  0.000  0.000  0.000  0.000  0.000 0.000 0.000  0.000 0.000
## 8 0.346 -0.760 -0.068 -0.546 -0.492 -0.017 -1.334 1.060 1.274 -1.052 0.247
## 9 0.247 -0.929 -0.177 -0.833 -0.541 -0.152 -1.399 1.138 1.512 -1.311 0.074
##       46     47     48     49     50
## 2  0.000  0.000  0.000  0.000  0.000
## 3  0.000  0.000  0.000  0.000  0.000
## 8 -0.929 -0.177 -0.833 -0.541 -0.152
## 9 -0.997 -0.406 -0.859 -0.767 -0.260
```
## 23 MBMGACPSSM
<font size="4">In this feature three different autocorrelation descriptors based on PSSM are adopted, which include: normalized Moreau-Broto autocorrelation, Moran autocorrelation and Geary autocorrelation
descriptors.Autocorrelation descriptor is a powerful statistical tool and defined based on the distribution of amino acid properties along the sequence, which measures the correlation between two
residues separated by a distance of d in terms of their evolution scores.</font>
#### Usage of this feature in PSSMCOOL package:
```
 w<-MBMGACPSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
```
##  [1] 0.377 0.296 0.337 0.305 0.182 0.436 0.340 0.287 0.311 0.258 0.262 0.315
## [13] 0.278 0.319 0.238 0.385 0.275 0.179 0.218 0.281 0.153 0.093 0.142 0.110
## [25] 0.043 0.253 0.126 0.110 0.120 0.094 0.118 0.099 0.091 0.152 0.066 0.168
## [37] 0.084 0.039 0.050 0.093 0.151 0.097 0.121 0.091 0.047 0.243 0.124 0.102
## [49] 0.118 0.101
```

## 24 LPC-PSSM
<font size="4">This feature uses Linear predictive coding algorithm for each column of PSSM. So for
producing this feature vector "lpc" function from "phontools" R-package is used which produces a 14-dimensional vector for each column, since PSSM has 20 column eventually it will be obtained a 20*14=280 dimensional feature vector for each PSSM.</font> 
#### Usage of this feature in PSSMCOOL package:
```
 w<-LPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
```
##  [1]  1.0000  0.9084  0.8390  0.5575  0.3737  0.4299  0.4070  0.4765  0.4369
## [10]  0.2830  0.3552  0.2212  0.1509  0.0496  1.0000  0.8380  0.7532  0.7463
## [19]  0.6135  0.3585  0.3674  0.2649  0.1610  0.2800  0.2347  0.2110  0.0763
## [28]  0.1057  1.0000  0.6738  0.6541  0.5619  0.5620  0.4152  0.4450  0.4506
## [37]  0.1943  0.1519  0.1683  0.0216 -0.1052  0.0215  1.0000  0.8553  0.8707
## [46]  0.6323  0.6506  0.7144  0.5384  0.5318
```

## 25 PSSM400
<font size="4">To generate this feature vector, for each of the standard amino acids, we find the positions containing that amino acid in the protein and separate the corresponding rows in the PSSM, to get a submatrix. Now, for the generated matrix, we calculate the average of its columns, and therefore, for each amino acid, a vector of length 20 is obtained. Finally, by putting these 20 vectors together, a feature vector of length 400 for each protein can be obtained. For example figure 16 shows the PSSM rows corresponding to amino acid S.</font>

![](vignettes/figures/pssm400.jpg)

#### Usage of this feature in PSSMCOOL package:
```
q<-pssm400(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
head(q, n = 50)
```
```
##  [1] 3.5000 2.2778 2.1667 2.1667 2.1667 2.6111 2.2778 2.8333 2.5556 2.0556
## [11] 2.1111 2.3889 2.1111 2.2778 1.9444 2.7222 2.3333 1.5000 2.1111 2.1667
## [21] 1.0556 1.7778 1.1111 0.8333 1.0000 1.1111 0.8889 0.8333 1.2222 0.7778
## [31] 0.8889 1.1667 1.0000 0.7778 0.7222 1.1667 1.0556 0.5556 0.8889 0.7778
## [41] 3.3889 2.2222 3.6667 3.6667 1.7778 2.8889 2.8333 3.7778 2.4444 2.0556
```

## 26 PSSM-BLOCK
<font size="4"> In this feature at first PSSM is divided to Blocks based on Number N which user imports.
Then for each Block the mean of columns is computed to get 20-dimensional vector, eventually by
appending these vectors to each other final feature vector is obtained.</font> 

#### Usage of this feature in PSSMCOOL package:
```
as<-PSSMBLOCK(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),5)
head(as, n = 50)
```
```
##  [1] 0.3773 0.2965 0.3371 0.3049 0.1816 0.4356 0.3395 0.2866 0.3109 0.2576
## [11] 0.2622 0.3150 0.2776 0.3187 0.2383 0.3846 0.2755 0.1791 0.2180 0.2808
## [21] 0.4032 0.2792 0.3606 0.2637 0.1736 0.5488 0.2988 0.3013 0.3996 0.2927
## [31] 0.3335 0.2389 0.3550 0.4466 0.2563 0.4349 0.2933 0.2153 0.2450 0.2929
## [41] 0.3530 0.3048 0.3146 0.3467 0.1921 0.3139 0.3719 0.2744 0.2156 0.2261
```

## 27 PSSM-SD
<font size="4">To generate this feature vector, at first for the column j, the sum of the total numbers in this column is calculated and denoted by <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT_j%7D">. Then, starting from the first row of this column, the numbers are added one by one together to reach a number less than or equal to 25 percent of <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT_j%7D">. Now the number of components used to calculate this sum is denoted by <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BI_j%5E1%7D"> and stored. Now, starting from the first row of this column again, the numbers are added one by one together to reach a number less than or equal to half of <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT_j%7D"> (50%), then we show the number of components to calculate this sum with <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BI_j%5E2%7D"> and store it. In the same way for column j, starting from the last row of this column, we start adding each elements together to reach a number less than or equal to 25% of the number <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT_j%7D">, and denote the number of these components by ![I_j^3](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BI_j%5E3%7D). In the next step, starting from the last row with summing of each element in this column to reach a number less than or equal to 50% of <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BT_j%7D">, the number ![I_j^4](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboldsymbol%7BI_j%5E4%7D) is also obtained. Therefore, 4 numbers are obtained for each column, and since the PSSM has 20 columns, for each protein a feature vector of length 80 is obtained. Figure 17 shows these steps schematically.</font>
![](vignettes/figures/screens/pssm-sd1.JPG)
#### Usage of this feature in PSSMCOOL package:
```
ww<-PSSM_SD(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ww, n = 50)
```
```
## [[1]]
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
## [1,]   54   23   15   13   31   10   14   22   18    44    52    16    62    77
## [2,]   70   57   70   32   46   13   54   66   89    71    79    32    84    94
## [3,]   88   83   90   88  104  121   87   86  110   107   108    63   109   116
## [4,]   71   57   70   32   46  110   54   66   89    72    79    32    84    97
##      [,15] [,16] [,17] [,18] [,19] [,20]
## [1,]    24    32    32    36    38    37
## [2,]    62    80    68    72    77    62
## [3,]   112   108    94    98   105    99
## [4,]    77    80    68    75    77    68
## 
## [[2]]
##  [1]  54  23  15  13  31  10  14  22  18  44  52  16  62  77  24  32  32  36  38
## [20]  37  70  57  70  32  46  13  54  66  89  71  79  32  84  94  62  80  68  72
## [39]  77  62  88  83  90  88 104 121  87  86 110 107 108  63 109 116 112 108  94
## [58]  98 105  99  71  57  70  32  46 110  54  66  89  72  79  32  84  97  77  80
## [77]  68  75  77  68
```

## 28 PSSM-SEG
<font size="4"> This feature, similar to the previous feature, divides each column into four parts and calculates the values for each column. Then, using the following equations, it calculates the values of Segmented Auto Covariance Features. The final feature vector length will be of length 100.</font>

![PSSM-Seg_{n,j}=\frac{1}{(I_j^n-m)}\sum_{i=1}^{I_j^n-m}(P_{i,j}-P_{ave,j})(P_{(i+m),j}-P_{ave,j})\\
(n=1,2,3,4 \ and \ j=1,...,20 \ and \ m\in \{1,2,...,11\})](https://latex.codecogs.com/svg.latex?%5Clarge%20PSSM-Seg_%7Bn%2Cj%7D%3D%5Cfrac%7B1%7D%7B%28I_j%5En-m%29%7D%5Csum_%7Bi%3D1%7D%5E%7BI_j%5En-m%7D%28P_%7Bi%2Cj%7D-P_%7Bave%2Cj%7D%29%28P_%7B%28i&plus;m%29%2Cj%7D-P_%7Bave%2Cj%7D%29%5C%5C%20%28n%3D1%2C2%2C3%2C4%20%5C%20and%20%5C%20j%3D1%2C...%2C20%20%5C%20and%20%5C%20m%5Cin%20%5C%7B1%2C2%2C...%2C11%5C%7D%29)


<font size="4">In the above equation, ![P_{ave,j}](https://latex.codecogs.com/svg.latex?%5Cinline%20P_%7Bave%2Cj%7D) represents the mean of column j in the PSSM and the number m is somehow a distance factor for each segment. Using the above equation, the feature of length 80 is obtained. Now the feature vector PSSM_AC is calculated using the previous factor m with length 20 and is added to the previous vector to get the final feature vector of length 100.</font>

![PSSM-AC_{m,j}=\frac{1}{(L-m)}\sum_{i=1}^{L-m}(P_{i,j}-P_{ave,j})(P_{(i+m),j}-P_{ave,j})](https://latex.codecogs.com/svg.latex?%5Clarge%20PSSM-AC_%7Bm%2Cj%7D%3D%5Cfrac%7B1%7D%7B%28L-m%29%7D%5Csum_%7Bi%3D1%7D%5E%7BL-m%7D%28P_%7Bi%2Cj%7D-P_%7Bave%2Cj%7D%29%28P_%7B%28i&plus;m%29%2Cj%7D-P_%7Bave%2Cj%7D%29)

<font size="4"> L Represents the total length of the protein. </font>
#### Usage of this feature in PSSMCOOL package:
```
q<-pssm_seg(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),3)
head(q, n = 50)

```
```
##  [1]  0.0342  0.0306  0.0257  0.0307  0.0316  0.0238  0.0170  0.0238  0.0348
## [10]  0.0308  0.0219  0.0308  0.0706  0.0358  0.0275  0.0358  0.0343  0.0322
## [19] -0.0013  0.0322 -0.0070  0.0398  0.0602  0.0649  0.0691  0.0482  0.0426
## [28]  0.0482  0.0091  0.0332  0.0243  0.0332  0.0526  0.0400  0.0364  0.0400
## [37]  0.0931  0.0559  0.0311  0.0555  0.1182  0.0769  0.0584  0.0769  0.0283
## [46]  0.0371  0.0175  0.0371  0.0287  0.0242
```

## 29 SOMA-PSSM
<font size="4"> This feature also considers each of the columns of the PSSM as a time series. If L represents the length of the protein, then the column j of the matrix can be thought of <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20y%28i%29%2C%20%5C%20i%3D1%2C2%2C...%2CL" style="width:120px;height:15px;"> as a time series. The SOMA algorithm is implemented in two steps using the following equations on the PSSM:
First, the moving average <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Csmall%20%5Coverline%7By_n%7D%28i%29" style="width:30px;height:15px;"> for the time series <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20y%28i%29" style="width:30px;height:15px;"> is calculated according to the following equation:</font>

<img src="https://latex.codecogs.com/svg.latex?%5Clarge%20%5Coverline%7By_n%7D%28i%29%3D%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bk%3D0%7D%5E%7Bn-1%7Dy%28i-k%29">

<font size="4">Where n is the size of the moving average window, and if it tends to zero, the moving average will tend to original series, in other words: if ![{n \to 0}](https://latex.codecogs.com/svg.latex?%5Cinline%20%7Bn%20%5Cto%200%7D) then ![{\overline{y_n}(i)\to y(i)}](https://latex.codecogs.com/svg.latex?%5Cinline%20%7B%5Coverline%7By_n%7D%28i%29%5Cto%20y%28i%29%7D). Next, for a moving average window size n which ![2\leq n<L](https://latex.codecogs.com/svg.latex?%5Cinline%202%5Cleq%20n%3CL), the second-order difference of the time series <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20y%28i%29" style="width:30px;height:15px;"> with respect to the moving average <img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Csmall%20%5Coverline%7By_n%7D%28i%29" style="width:30px;height:15px;"> is defined according to the following equation:</font>

![\sigma_{MA}^2=\frac{1}{L-n}\sum_{i=n}^L[y(i)-\overline{y_n}(i)]^2](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Csigma_%7BMA%7D%5E2%3D%5Cfrac%7B1%7D%7BL-n%7D%5Csum_%7Bi%3Dn%7D%5EL%5By%28i%29-%5Coverline%7By_n%7D%28i%29%5D%5E2)

<font size="4"> The number n must be smaller than the length of the smallest protein in the database under study. In the paper used by this algorithm, the length of the smallest protein is 10 and therefore the number n will vary from 2 to 9, so according to above Equation, by putting the numbers ![\sigma_{MA}^2](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Csigma_%7BMA%7D%5E2) next to each other, 8 numbers are obtained for each column, and therefore the final feature vector will be of length 160.</font>

#### Usage of this feature in PSSMCOOL package:
```
 w<-SOMA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
```
##  [1] 0.2362 0.1898 0.2082 0.2175 0.1262 0.3151 0.2126 0.1880 0.2021 0.1703
## [11] 0.1830 0.1946 0.1593 0.2495 0.1656 0.2278 0.1300 0.1307 0.1247 0.1810
## [21] 0.7849 0.5708 0.7035 0.6571 0.3412 1.1477 0.6830 0.5934 0.6513 0.5323
## [31] 0.6020 0.5920 0.4988 0.8121 0.4674 0.7992 0.4326 0.3434 0.3525 0.5541
## [41] 1.6475 1.1293 1.4382 1.2852 0.6546 2.4681 1.3966 1.2121 1.3366 1.1003
```
## 30 SVD-PSSM
<font size="4"> Singular value decomposition is a general-purpose matrix factorization approach
that has many useful applications in signal processing and statistics. In this feature ![SVD](https://latex.codecogs.com/svg.latex?%5Cinline%20SVD) is
applied to a matrix representation of a protein aimed to reduce its dimensionality.
Given an input matrix Mat with dimensions ![N*M \enspace SVD](https://latex.codecogs.com/svg.latex?%5Cinline%20N*M%20%5Censpace%20SVD) is used to calculate its factorization
of the form: ![Mat=U\Sigma V](https://latex.codecogs.com/svg.latex?%5Cinline%20Mat%3DU%5CSigma%20V) where ![\Sigma](https://latex.codecogs.com/svg.latex?%5Cinline%20%5CSigma) is a diagonal matrix whose diagonal
entries are known as the singular values of Mat. The resulting descriptor is the ordered set of singular values: ![SVD\in\mathcal{R}^L](https://latex.codecogs.com/svg.latex?%5Cinline%20SVD%5Cin%5Cmathcal%7BR%7D%5EL) where ![L=min(M,N)](https://latex.codecogs.com/svg.latex?%5Cinline%20L%3Dmin%28M%2CN%29). since the PSSM has 20 columns, the final feature vector would be of length 20 .</font>
#### Usage of this feature in PSSMCOOL package:
```
 w<-SVD_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 20)
```
```
##  [1] 16.312  8.469  5.364  4.757  4.254  3.710  3.687  3.207  2.943  2.687
## [11]  2.565  2.262  2.244  2.073  1.832  1.755  1.621  1.507  1.441  1.061
```

# useCase implementation in R
#### # *Installing PSSMCOOL and loading it*
```
# install.packages("PSSMCOOL")
# library(PSSMCOOL)
```
#### # *setting up working environment and downloading necessary files from GitHub*
```
current_directory <- "F:\\article400\\javad2\\saber"
setwd(current_directory)
```
#### # *Downloading the required PSSM files*
```
pssm_url <- 'https://github.com/BioCool-Lab/PSSMCOOL/raw/main/classification-code-data/all_needed_pssms90.zip' 
download.file(pssm_url, './all_needed_pssm90.zip', method = 'auto', quiet = FALSE) 
unzip('all_needed_pssm90.zip', exdir = 'all_needed_pssm90') 
PSSM_directory <- 'all_needed_pssm90/all_needed_pssms90/'
```
#### # *Downloading positive data and loading it to R*
```
url <- "https://raw.githubusercontent.com/BioCool-Lab/PSSMCOOL/main/classification-code-data/positive.csv"
download.file(url, './PositiveData.csv')
positive_data <- read.csv("./PositiveData.csv", header = TRUE)
```
#### # *Downloading negative data and loading it to R*
```
url <- "https://raw.githubusercontent.com/BioCool-Lab/PSSMCOOL/main/classification-code-data/negative.csv" 
download.file(url, './NegativeData.csv')
negative_data <- read.csv("./NegativeData.csv", header = TRUE)
```
####  ###############################---*Positive feature extraction*---####################################
####  # *Feature extraction*
```
positiveFeatures<- c() 
for(i in 1:dim(positive_data)[1]) { 
  ff<-FPSSM2(paste0(PSSM_directory, positive_data[i,1],'.fasta.pssm'), 
             paste0(PSSM_directory, positive_data[i,2],'.fasta.pssm'), 20) 
  positiveFeatures<-rbind(positiveFeatures, ff) 
}
```
#### # *Adding row names and class*
```
positiveFirstColumn <- c()
for(i in 1:dim(positive_data)[1]) { 
  dd <- paste(positive_data[i,1], '-' ,positive_data[i,2]) 
  positiveFirstColumn <- rbind(positiveFirstColumn, dd) 
}
```
```
pos_class <- rep("Interaction", dim(positiveFeatures)[1]) 
positiveFeatures2 <- cbind(positiveFirstColumn, positiveFeatures, pos_class) 
```
####   ###############################---*Negative feature extraction*---#####################################
####   # *Feature extraction*
```
negativeFeatures <- c() 
for(i in 1:dim(negative_data)[1]) { 
  ff2<-FPSSM2(paste0(PSSM_directory, negative_data[i,1],'.fasta.pssm'), 
              paste0(PSSM_directory, negative_data[i,2],'.fasta.pssm'), 20) 
  negativeFeatures<-rbind(negativeFeatures, ff2) 
}
```
#### # *Adding row names and class*
```
negativeFirstColumn <- c() 
for(i in 1:dim(negative_data)[1]) { 
  dd2 <- paste(negative_data[i,1], '-' ,negative_data[i,2]) 
  negativeFirstColumn <- rbind(negativeFirstColumn, dd2) 
} 
```
```
neg_class <- rep("Non.Interaction", dim(negativeFeatures)[1]) 
negativeFeatures2 <- cbind(negativeFirstColumn, negativeFeatures, neg_class) 
```
#### # *Merging two feature vectors*
```
mainDataSet <- rbind(positiveFeatures2, negativeFeatures2) 
```

####  ###############################---*Preparing data set for model training*---###########################
####  # *In the following we are going to carry out classification on the data we have prepared so far (mainDataSet)*
####  # *First we need to install and load caret package and its dependencies*
```
install.packages('caret', dependencies = TRUE) 
library(caret) 
bmp.R2.submission.data.df <- as.data.frame(mainDataSet) 
colnames(bmp.R2.submission.data.df)[1] <- "interactions" 
dim(bmp.R2.submission.data.df)#1730  102 
```
#### # *Assigning the Uniprot IDs for each protein pairs to the row name*
```
rownames(bmp.R2.submission.data.df) <- bmp.R2.submission.data.df$interactions
```
#### # *Removing the Uniprot IDs*
```
bmp.R2.submission.data.df <-bmp.R2.submission.data.df[,-1] 
View(bmp.R2.submission.data.df)
colnames(bmp.R2.submission.data.df) <- c(paste0('Frt', 1: dim(positiveFeatures)[2]), 'Class') 
dim(bmp.R2.submission.data.df)#1730  101 
table(bmp.R2.submission.data.df$Class)
```
######   # *Interaction--Non-Interaction*
######   #     865-------------865 
```
bmp.R2.submission.data.df$Class <- 
  as.factor(bmp.R2.submission.data.df$Class) 
write.csv(bmp.R2.submission.data.df, 'DataSet.csv') 
```
#### ###############################---*Training model with three classifier*---#############################

#### # *setting.the.trainControl*
```
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
```
#### # ####################---*setting cross validation parameters*---########################################
```
trainControl.for.PSSM <- setting.the.trainControl.3() 
```

#### # ###################---*10-fold cross-validation using "Bagged CART (treebag)" classifier*---#######################
```
cross.validation.bulit.model.treebag <- 
  train(Class ~ ., data = bmp.R2.submission.data.df, 
        method = "treebag", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE) 
print(cross.validation.bulit.model.treebag$results) 
```
######  # parameter---Accuracy-----Kappa-----AccuracySD----KappaSD
######  # 1---none---0.9965351---0.9930707---0.005582867---0.01116413


#### # #################---*10-fold cross-validation using "Single C5.0 Tree (C5.0Tree)" classifier*---#####################
```
cross.validation.bulit.model.C5.0Tree <- 
  train(Class ~ ., data = bmp.R2.submission.data.df, 
        method = "C5.0Tree", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE) 
print(cross.validation.bulit.model.C5.0Tree$results) 
```
######  # parameter---Accuracy----Kappa----AccuracySD----KappaSD
######  # 1----none---0.9976911---0.9953822---0.004028016---0.008056142

#### # ##############---*10-fold cross-validation using "Partial Least Squares (pls)" classifier*---#########################
```
cross.validation.bulit.model.pls <-
  train(Class ~ ., data = bmp.R2.submission.data.df, 
        method = "pls", 
        trControl = trainControl.for.PSSM, 
        verbose = FALSE) 
print(cross.validation.bulit.model.pls$results) 
```
######  #---ncomp---Accuracy----Kappa----AccuracySD----KappaSD
######  #1----1----0.5005948---0.00231924--0.01330889--0.02666192
######  #2----2----0.5070671---0.01413974--0.03944775--0.07848283
######  #3----3----0.5324142---0.06473979--0.02790055--0.05655902

