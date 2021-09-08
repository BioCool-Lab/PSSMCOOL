# PSSMCOOL

This package includs all features that has been extracted from Position-Specific Scoring Matrix (PSSM) in POSSUM website: http://possum.erc.monash.edu/ and some of other features
from several articles. each function in this package corresponds to one feature that has extracted from PSSM Matrix.

###### **Install from GitHub :  devtools::install_github("Alireza9651501005/PSSMCOOL")**

vignettes/figures/screens/pssm_ac.JPG


# Preface
<br></br>
<font size="4"> Feature extraction or feature encoding is a fundamental step in the construction of high-quality machine learning-based models. Specifically, this is a key step for determining the effectiveness of the trained models in bioinformatics applications. In the last two decades, a variety of feature encoding schemes have been proposed in order to exploit useful patterns from protein sequences. Such schemes are often based on sequence information or physicochemical properties of amino acids. Although direct features derived from sequences themselves (such as amino acid compositions, dipeptide compositions, and counting of k-mers) are regarded as essential for training models, an increasing number of studies have shown that evolutionary information in the form of PSSM profiles is much more informative than sequence information alone. Accordingly, PSSM-based feature descriptors have been commonly used as indispensable primary features to construct models, filling a major gap in the current bioinformatics research. For example, PSSM-based feature descriptors have successfully improved the prediction performance of structural and functional properties of proteins across a wide spectrum of bioinformatics applications. These predictions can be applied for protein fold recognition and the prediction of protein structural classes, protein-protein interactions, protein subcellular localization, RNA-binding sites, and protein functions. But at the same time, there is no comprehensive, simple tool in the R programming language for extracting all of these features from the PSSM and displaying it in the output. PSSMCOOL package here is developed in R for these purposes. First, in Figure 1 which is a table that lists all of the features implemented in this package with their feature-lengths is brought. Then each one of these features will be explained in full detail.</font>

<br></br>

<font size="4"> PSSMCOOL Package is currently available on CRAN website:</font>

[https://cran.r-project.org/web/packages/PSSMCOOL/index.html](https://cran.r-project.org/web/packages/PSSMCOOL/index.html)

<br></br>

<font size="4"> for issues about this package :</font>

[https://github.com/Alireza9651501005/PSSMCOOL/issues](https://github.com/Alireza9651501005/PSSMCOOL/issues)

<br></br>
![](vignettes/figures/feature_table.jpg)
<pre>                  Figure 1: List of implemented features in PSSMCOOL package </pre>

<br>*feature vector length depends on the choice of parameter</br>  **these features produce Matrix of features which its dimension depends on choice of parameter

<br></br>

```{r}
library(PSSMCOOL)
```

# PSSM-AC
<font size="4">This feature, which stands for auto-covariance transformation, for column j, calculates the average of this column as shown in Figure 2, Then subtracts the resulting number from the elements on the rows i and i + g of this column, and finally multiplies them and calculates the sum by changing the variable i from 1 to L-g. Because the variable j changes between 1 and 20 and the variable g changes between 1 and 10, eventually a feature vector of length 200 will be obtained [@zou2013accurate].</font>

<br></br>
![](vignettes/figures/pssm_ac.jpg)
<br></br>
<pre>                   Figure 2: process of extracting PSSM-AC feature vector from PSSM Matrix </pre>

![](vignettes/figures/screens/pssm_ac.JPG)

# DPC-PSSM
<font size="4">This feature stands for dipeptide composition, which multiplies the values that are located in two consecutive rows and two different columns. Having calculated these values for different rows and columns, they are summed. Next, for both columns, the sum is divided by L-1. Since the result depends on two different columns, eventually a feature vector of length 400, according to Figure 3 and following equation ,will be obtained.</font>
<br></br>

![](vignettes/figures/dpc-pssm.jpg)
<pre>                   Figure 3: process of extracting DPC-PSSM feature vector from PSSM Matrix </pre>
<br></br>

![](vignettes/figures/screens/dpc_pssm.JPG)


# Trigram-PSSM
<font size="4">This feature vector is of length 8000, which is extracted from the PSSM. If we multiply elements available in three consecutive rows and three different columns of the PSSM by each other, and apply this to all rows (all three consecutive rows) and then sum these numbers, eventually one of the elements of feature vector with length 8000 corresponding to the three selected columns will be obtained. Because we have 20 different columns, The final feature vector will be of length 8000 = 20 * 20 * 20. Figure 4: shows these steps. For example, in this figure for three marked rows and columns, the numbers obtained from the intersection of these rows and columns marked with a blue dotted circle around them, are multiplied to each other. </font>

![](vignettes/figures/trigram.jpg)

<pre>                   Figure 4: process of extracting trigram-PSSM feature vector from PSSM </pre>
<br></br>

![](vignettes/figures/screens/trigram.JPG)
<br></br>

# Pse-PSSM
<font size="4"> The length of this feature vector is 320. The first 20 numbers of this feature vector are the mean of 20 columns in PSSM, and the next values for each column are the mean squares of the difference between the elements of row i and i + lag in this column. Because the lag value varies between 1 and 15, the final feature vector will have a length of 320. Figure 5: and following equation shows the process of this function and the corresponding mathematical equation, respectively.</font>

![](vignettes/figures/pse-pssm.jpg)
<pre>                  Figure 5: process of extracting Pse-PSSM feature vector from PSSM </pre>
<br></br>

![](vignettes/figures/screens/pse_pssm.JPG)

# k-Separated-bigram-PSSM
<font size="4"> This feature is almost identical to the DPC feature, and in fact, the DPC feature is part of this feature (for k = 1) and for two different columns, it considers rows that have distance k.</font>

![](vignettes/figures/k-separated.jpg)
<pre>                  Figure 5: process of extracting Pse-PSSM feature vector from PSSM </pre>
<br></br>

![](vignettes/figures/screens/k-separated.JPG)


# EDP-EEDP-MEDP
<font size="4"> In this group of features, in order to use uniform dimensions to show proteins of different lengths, in the first step, the average evolutionary score between adjacent residues is calculated using the following equations:</font>

![](vignettes/figures/screens/edp.JPG)
![](vignettes/figures/EEDP.jpg)
<pre>                  Figure 7: process of extracting EDP-EEDP-MEDP feature vectors from PSSM </pre>
<br></br>


```{r eedp,echo=FALSE,fig.cap="Figure 7: process of extracting EDP-EEDP-MEDP feature vectors from PSSM",out.width = '70%'}
knitr::include_graphics("figures/EEDP.jpg")
```
<br></br>
<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 as<-EDP_MEDP(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GS61.txt.pssm"))
head(as, n = 50)
```
<br></br>

# AB-PSSM
<font size="4"> This feature consists of two types of feature vectors. At first, each protein sequence is divided into 20 equal parts, each of which is called a block, and in each block, the row vectors of the PSSM
related to that block are added together. The resulting final vector is divided by the length of that block, which is equal to 5% of protein length. Finally, by placing these 20 vectors side by side, the first feature vector of length 400 is obtained. The second feature for each amino acid in each column is the average of the positive numbers in that column and for each block, and these 20 values, corresponding to 20 blocks, are placed next to each other, and therefore for each of the 20 types of amino acids, a vector of length 20 is obtained, and by placing these together,the second feature vector of length 400, is obtained. Figure 8 represents this process [@cheol2010position].</font>

```{r ab-pssm,echo=FALSE,fig.cap="Figure 8: process of extracting AB-PSSM feature vectors from PSSM",out.width = '70%'}
knitr::include_graphics("figures/AB-PSSM.jpg")
```
<br></br>
<font size="4"> usage of this feature in PSSMCOOL package:</font>

```{r}
  zz<- AB_PSSM(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"))
head(zz[1], n = 50)
```
<br></br>

# AATP-TPC
<font size="4"> In this feature, at first, a TPM matrix is constructed from the PSSM, which has represented by a vector corresponding to the following equation:</font>

\begin{equation}
Y_{TPM}=(y_{1,1},y_{1,2},...,y_{1,20},...,y_{i,1},...,y_{i,20},...,y_{20,1},...,y_{20,20})^T
\end{equation}

<br></br>
<font size="4">Where the components are as follows:</font>
<br></br>

\begin{equation}
y_{i,j}=(\sum_{k=1}^{L-1}P_{k,i}\times P_{k+1,j})/(\sum_{j=1}^{20}\sum_{k=1}^{L-1}P_{k+1,j}\times P_{k,i}) \\  1\leq{i,j}\leq{20}
\end{equation}
<br></br>
<font size="4"> In the above equation, the numerator is the same as the equation related to DPC-PSSM feature without considering its coefficient. By placing these components together, a TPC feature vector of length 400 is obtained, and if we add the AAC feature vector of length 20 which is the average of columns of the PSSM to the beginning of this vector, AATP feature vector of length 420 is obtained [@zhang2012using].</font>
<br></br>
<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 as<-AATP_TPCC(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GQS7.txt.pssm"))
head(as, n = 50)
```
<br></br>

# CS-PSe-PSSM
<font size="4"> This feature consists of a combination of several types of features, and in general, the obtained feature vector would be of length 700. Here all parts of this feature are described separately.</font>
<br></br>

## CSAAC

<font size="4"> First, the consensus sequence is obtained from the PSSM according to the following equation, then from this consensus sequence next feature vectors are obtained </font>
<br></br>
\begin{equation}
\alpha(i)=argmax{P_{i,j}:\ 1\leq j\leq {20}} \quad ,1\leq i\leq L
\end{equation}

<font size="4">where $\alpha(i)$ is the index for the largest element in row i  of PSSM, and the ith component in the consensus sequence equal to the $\alpha(i)$'th amino acid in the standard amino acid alphabet which the column names of PSSM are labeled by them. Now, using the following equation, the feature vector of length 20 is obtained:</font>
<br></br>
\begin{equation}
CSAAC=\frac{n(j)}{L}\quad ,1\leq j\leq{20}
\end{equation}
<font size="4"> Here $n(j)$ shows the number of j-th amino acid occurrences in the consensus sequence </font>
<br></br>

## CSCM

<font size="4">This feature vector is obtained using the following equation from the consensus sequence:</font>

\begin{equation}
CSCM=\frac{\sum_{j=1}^{n_i}n_{i,j}}{L(L-1)}\quad ,1\leq i\leq{20}\quad ,1\leq j\leq L
\end{equation}
<font size="4"> Here $n(i)$ shows the number of i-th amino acid occurrences in the consensus sequence and $n_{i,j}$ indicates the j-th position of the i-th amino acid in the consensus sequence.</font>

## Segmented PsePSSM features

<font size="4">Here the PSSM is divided into n segments, which corresponds to dividing the initial protein sequence into n segments,  if n = 2:</font>

\begin{equation}
L_1=round(L/2) \quad , L_2=L-L_1
\end{equation}

<font size="4">Where L represents the length of the initial protein and indexed L's indicate the length of the first and second segments, respectively. Now, using the following equations, the feature vector components of length 200 are obtained</font>

\begin{equation}
\alpha_j^\lambda=\left\{\begin{array}{ll}\frac{1}{L_1}\sum_{i=1}^{L_1}P_{i,j} & j=1,2,...,20,\lambda=0 \\
\frac{1}{{L_1}-\lambda}\sum_{i=1}^{{L_1}-\lambda}(P_{i,j}-P_{i+\lambda,j})^2 & j=1,2,...,20,\lambda=1,2,3,4\end{array}\right.
\end{equation}
\begin{equation}
\beta_j^\lambda=\left\{\begin{array}{ll}\frac{1}{L-L_1}\sum_{i=L_1+1}^{L}P_{i,j} & j=1,2,...,20,\lambda=0 \\
\frac{1}{L-L_1-\lambda}\sum_{i=L_1+1}^{L-\lambda}(P_{i,j}-P_{i+\lambda,j})^2 & j=1,2,...,20,\lambda=1,2,3,4\end{array}\right.
\end{equation}
<font size="4">And if n = 3 then we have:</font>

\begin{equation}
L_1=round(L/3) \quad , L_2=2L_1 \quad ,L_3=L-2L_1
\end{equation}

<font size="4">Therefore, using the following equations, the  components of a feature vector with length 180 are obtained as follows:</font>

\begin{equation}
\theta_j^\lambda=\left\{\begin{array}{ll}\frac{1}{L_1}\sum_{i=1}^{L_1}P_{i,j} & j=1,2,...,20,\lambda=0 \\
\frac{1}{L_1-\lambda}\sum_{i=1}^{L_1-\lambda}(P_{i,j}-P_{i+\lambda,j})^2 & j=1,2,...,20,\lambda=1,2\end{array}\right.
\end{equation}


\begin{equation}
\mu_j^\lambda=\left\{\begin{array}{ll}\frac{1}{L_1}\sum_{i=L_1+1}^{2L_1}P_{i,j} & j=1,2,...,20,\lambda=0 \\
\frac{1}{L_1-\lambda}\sum_{i=L_1+1}^{2L_1-\lambda}(P_{i,j}-P_{i+\lambda,j})^2 & j=1,2,...,20,\lambda=1,2\end{array}\right.
\end{equation}



\begin{equation}
v_j^\lambda=\left\{\begin{array}{ll}\frac{1}{L-2L_1}\sum_{i=2L_1+1}^{L}P_{i,j} & j=1,2,...,20,\lambda=0 \\
\frac{1}{L-2L_1-\lambda}\sum_{i=2L_1+1}^{L-\lambda}(P_{i,j}-P_{i+\lambda,j})^2 & j=1,2,...,20,\lambda=1,2\end{array}\right.
\end{equation}

<font size="4">In total, using the previous feature vector, a feature vector of length 380 is obtained for this group.</font>

<br></br>

## Segmented ACTPSSM features
<font size="4">In this group, using the previous equations and the following equations, the feature vector of length 280 could be obtained. when n=2:</font>

\begin{equation}
AC1_j^{lg}=\frac{1}{L_1-lg}\sum_{i=1}^{L_1-lg}(P_{i,j}-\alpha_j^0)(P_{i+lg,j}-\alpha_j^0)\\ 
AC2_j^{lg}=\frac{1}{L-L_1-lg}\sum_{i=L_1+1}^{L-lg}(P_{i,j}-\beta_j^0)(P_{i+lg,j}-\beta_j^0)\\
j=1,2,...,20,lg=1,2,3,4
\end{equation}

<font size="4">when n=3:</font>

\begin{equation}
AC1_j^{lg}=\frac{1}{L_1-lg}\sum_{i=1}^{L_1-lg}(P_{i,j}-\theta_j^0)(P_{i+lg,j}-\theta_j^0)\\ 
AC2_j^{lg}=\frac{1}{L_1-lg}\sum_{i=L_1+1}^{2L_1-lg}(P_{i,j}-\mu_j^0)(P_{i+lg,j}-\mu_j^0)\\
AC3_j^{lg}=\frac{1}{L-2L_1-lg}\sum_{i=2L_1+1}^{L-lg}(P_{i,j}-v_j^0)(P_{i+lg,j}-v_j^0)\\
j=1,2,...,20,lg=1,2
\end{equation}

<br></br>
<font size="4">If we connect all these feature vectors together, we will get a feature vector with a length of 700, which is reduced by PCA method and is used as input for the support vector machine classifier [@liang2015prediction].</font>
<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 A<-CS_PSe_PSSM(system.file("extdata", "C7GSI6.txt.pssm", package="PSSMCOOL"),"total")
head(A, n = 50)
```
<br></br>

# D-FPSSM/S-FPSSM

<font size="4">If we sum the numbers of each column in the PSSM, we get a feature vector of length 20 as follows:</font>
<br></br>
\begin{equation}
D=(d_1,d_2,...,d_{20})
\end{equation}

<font size="4">If we remove the negative elements somehow of the PSSM and call the resulting new matrix FPSSM and then calculate this feature vector for FPSSM, the components of this vector will depend on the length of the original protein, so to eliminate this dependency, we normalize the components of this vector using the following equation:</font>

\begin{equation}
d_i=\frac{d_i-min}{max\times L}
\end{equation}

<font size="4">Where, min and max represent the smallest and largest values of the previous vector components, respectively, and L represents the length of the original protein. The second feature vector with length 400 is obtained as follows:</font>

\begin{equation}
S=(s_1^{(1)},s_2^{(1)},...,s_{20}^{(1)},s_1^{(2)},s_2^{(2)},...,s_{20}^{(2)},...,s_1^{(20)},s_2^{(20)},...,s_{20}^{(20)})
\end{equation}
<br></br>
<font size="4">If we name the columns of the FPSSM from $a_1$ to $a_{20}$ in the order from left to right, 
then $S_j^{(i)}$ is equal to the sum of those members in the j-th column in the FPSSM whose corresponding row amino acid is equal to $a_i$. Figure 9 schematically shows these steps [@zahiri2013ppievo]:</font>
<br></br>

```{r fpssm,echo=FALSE,fig.cap="Figure 9: process of making FPSSM and extracting corresponding feature vectors",out.width = '70%'}
knitr::include_graphics("figures/s-fpssm.jpg")
```

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 q<-FPSSM(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),20)
head(q, n = 50)
```
<br></br>

# SCSH2

<font size="4"> To generate this feature vector, the consensus sequence corresponding to the protein sequence is extracted using the PSSM. Then, by placing these two sequences next to each other, a matrix with dimensions of 2 * L will be created. In the next step, each component in the upper row of this matrix is connected to two components in the lower row of this matrix, and thus a graph similar to a bipartite graph could be created.
Now in this graph, each path of length 2 specifies a 3-mer and each path of length 1 denotes a 2-mer corresponding to these two sequences. Now if we consider a table consisting of two rows and 8000 columns so that the first row contains all possible 3-mers of 20 amino acids, then for every 3-mer obtained from this graph, we put number 1 below the corresponding cell With that 3-mer in the aforementioned table and 0 in other cells. so This gives us a vector of length 8000. For the 2-mers obtained from this graph, a vector of length 400 is obtained in a similar way. figures 10, 11 show these processes [@zahiri2014locfuse]. </font>
<br></br>
<br></br>


```{r scsh2,echo=FALSE,fig.cap="Figure 10: process of extracting scsh2 feature vector",out.width = '70%'}
knitr::include_graphics("figures/SCSH2.jpg")
```

<br></br>
<br></br>

```{r scshtable,echo=FALSE,fig.cap="Figure 11: tables of all 2-mers and all 3-mers",out.width = '70%'}
knitr::include_graphics("D:/Edit_thesis/thesis_pics/scshtable.jpg")
```

<br></br>
<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 zz<- scsh2(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"),2)
head(zz, n = 200)
```
<br></br>

# RPSSM

<font size="4">If we represent the PSSM as follows:</font>

\begin{equation}
D=(P_A,P_R,P_N,P_D,P_C,P_Q,P_E,P_G,P_H,P_I,P_L,P_K,P_M,P_F,P_P,P_S,P_T,P_W,P_Y,P_V)
\end{equation}

<font size="4">The indices will show the standard 20 amino acids. If we assume that our primary protein has length L, each of the above columns is as follows:</font>

\begin{equation}
P_A=(P_{1,A},P_{2,A},...,P_{L,A})^T
\end{equation}

<br></br>
<font size="4">Now, using the following equations, we merge the columns of the PSSM and obtain a matrix with dimensions $L\times 10$:</font>

<br></br>


\begin{equation}
P_1=\frac{P_F+P_Y+P_W}{3},P_2=\frac{P_M+P_L}{2}P_3=\frac{P_I+P_V}{2}\\
P_4=\frac{P_A+P_T+P_S}{2},P_5=\frac{P_N+P_H}{2},P_6=\frac{P_Q+P_E+P_D}{3}\\
P_7=\frac{P_R+P_K}{2},P_8=P_C,P_9=P_G,P_{10}=P_P
\end{equation}

<br></br>

\begin{equation}
RD=\begin{pmatrix}
-&1&2&3&4&5&6&7&8&9&10\\
a_1&p_{1,1}&p_{1,2}&p_{1,3}&p_{1,4}&p_{1,5}&p_{1,6}&p_{1,7}&p_{1,8}&p_{1,9}&p_{1,10}\\
a_2&p_{2,1}&p_{2,2}&p_{2,3}&p_{2,4}&p_{2,5}&p_{2,6}&p_{2,7}&p_{2,8}&p_{2,9}&p_{2,10}\\
\vdots&\vdots&\vdots&\vdots&\vdots&\vdots&\vdots&\vdots&\vdots&\vdots\\
a_L&p_{L,1}&p_{L,2}&p_{L,3}&p_{L,4}&p_{L,5}&p_{L,6}&p_{L,7}&p_{L,8}&p_{L,9}&p_{L,10}
\end{pmatrix}
\end{equation}

<br></br>

<font size="4">Now using this new matrix we get a feature vector of length 10 as follows: </font>

\begin{equation}
D_s=\frac{1}{L}\sum_{i=1}^L (p_{i,s}-\overline p_s)^2
\end{equation}

<font size="4"> where </font>

\begin{equation}
\overline p_s=\frac{1}{L}\sum_{i=1}^L p_{i,s} \quad ,s=1,2,...,10,\ i=1,2,...,L \ ,p_{i,s} \in RD
\end{equation}

<font size="4"> now using following equations; will create a feature vector of length 100 and by combining the feature vector of length 10 mentioned previously, the final feature vector of length 110 will be created [@ding2014protein].</font>

\begin{equation}
\begin{aligned}
x_{i,i+1}&=(p_{i,s}-\frac{p_{i,s}+p_{i+1,t}}{2})^2+(p_{i+1,t}-\frac{p_{i,s}+p_{i+1,t}}{2})^2\\
&=\frac{(p_{i,s}-p_{i+1,t})^2}{2}\quad i=1,2,...,L-1 \quad ,s,t=1,2,...,10 
\end{aligned}
\end{equation}

<br></br>

\begin{equation}
\begin{aligned}
D_{s,t}&=\frac{1}{L-1}\sum_{i=1}^{L-1}x_{i,i+1} \\
&=\frac{1}{L-1}\sum_{i=1}^{L-1}[(p_{i,s}-\frac{p_{i,s}+p_{i+1,t}}{2})^2+(p_{i+1,t}-\frac{p_{i,s}+p_{i+1,t}}{2})^2] \\
&=\frac{1}{L-1}\sum_{i=1}^{L-1}\frac{(p_{i,s}-p_{i+1,t})^2}{2} \quad ,s,t=1,2,...,10
\end{aligned}
\end{equation}

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 w<-rpssm(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
<br></br>

# CC-PSSM
<font size="4">This feature, which is similar to the PSSM-AC feature, stands for cross-covariance transformation. for column $j_1$, Calculates the average of this column as shown in Figure 12, and then subtract the result from the number on the i-th row in this column. Similarly, the feature calculates the average for the column $j_2$ and then subtracts the resulting number from the value on row i + g of this column and finally multiplies them. By changing the variable i from 1 to L-g, it calculates the sum of these, because the variable $j_1$ changes between 1 and 20 and the variable $j_2$ changes in the same interval (1,20) except for the number selected for the variable $j_1$, eventually feature vector of length 380 will be obtained [@dong2009new]. </font>

<br></br>
<br></br>


```{r ccpssm,echo=FALSE,fig.cap="Figure 12: process of extracting CC-PSSM feature vector",out.width = '70%'}
knitr::include_graphics("figures/cc-pssm.jpg")
```
<br></br>

\begin{equation}
CC-PSSM_{j1,j2}=\sum_{i=1}^{L-g}(P_{i,j1}-\overline {P_{j1}})(P_{i+g,j2}-\overline {P_{j2}})\\
1\leq j1,j2\leq 20
\end{equation}

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 aa<-pssm_cc(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),18)
head(aa, n = 50)
```
<br></br>

# Discrete cosine transform
<font size="4">Discrete cosine transforms can be described as follows:</font>

\begin{equation}
DCT(u,v)=\rho(u)\rho(v)\sum_{x=0}^{M-1}\sum_{y=0}^{N-1}f(x,y)\cos\frac{(2x+1)u\pi}{2M}\cos\frac{(2y+1)v\pi}{2N}\\
0\leq{u}\leq{M-1} \quad ,0\leq{v}\leq{N-1}
\end{equation}

<font size="4"> where:</font>

\begin{equation}
\rho(u)=\left\{\begin{array}{ll}\sqrt{\frac{1}{M}} &  u=0\\ 
\sqrt{\frac{2}{M}} & 1\leq{u}\leq{M-1}\end{array}\right.\\
\rho(v)=\left\{\begin{array}{ll}\sqrt{\frac{1}{N}} & v=0\\ 
\sqrt{\frac{2}{N}} & 1\leq{v}\leq{N-1}\end{array}\right.
\end{equation}

<font size="4"> In above Equation, the matrix $f(x,y)\in P^{N\times M}$ is the input signal and here represents the PSSM with dimensions $N\times 20$. According to Equation====, it is clear that the length of the resulting feature vector depends on the length of the original protein, so in most articles that have used this feature vector, the final feature vector DCT, which encodes a protein sequence by choosing the first 400 coefficients is obtained [@wang2017advancing].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 as<-Discrete_Cosine_Transform(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)
```
<br></br>

# Discrete Wavelet Transform
<font size="4">Wavelet transform (WT) is defined as the signal image $f(t)$ on the wavelet function according to the following equation:</font>

\begin{equation}
T(a,b)=\sqrt{\frac{1}{a}}\int_0^t f(t)\psi{(\frac{t-b}{a})}dt
\end{equation}

<font size="4">Where a is a scale variable and b is a transition variable. $\psi{(\frac{t-b}{a})}$ is the analyze wavelet function.  $T(a,b)$ is the conversion factor found for both specific locations on the signal as well as specific wavelet periods. Discrete wavelet transform can decompose amino acid sequences into coefficients in different states and then remove the noise component from the profile.Assuming that the discrete signal $f(t)$ is equal to $x[n]$, where the length of the discrete signal is equal to N, we have the following equations:</font>

\begin{equation}
y_{j,low}[n]=\sum_{k=1}^N x[k]g[2n-k]\\
y_{j,high}[n]=\sum_{k=1}^N x[k]h[2n-k]
\end{equation}

<font size="4">In these equations, g is the low-pass filter and h is the high-pass filter. $y_{low}[n]$ is the approximate coefficient (low-frequency components) of the signal and $y_{high}[n]$ is the exact coefficient (high-frequency components) of the signal.This decomposition is repeated to further increase the frequency resolution, and the approximate coefficients are decomposed with high and low pass filters and then are sampled lower. By increasing the level of j decomposition, we can see more accurate characteristics of the signal. We use level 4 DWT and calculate the maximum, minimum, mean, and standard deviation of different scales (4 levels of both coefficients High and low frequency). Because high-frequency components have high noise, only low-frequency components are more important.A schematic diagram of a level 4 DWT is shown in Figure 13:</font>
<br></br>
<br></br>

```{r dwt,echo=FALSE,fig.cap="Figure 13: Schematic diagram of a DWT with 4 levels",out.width = '70%'}
knitr::include_graphics("figures/dwt.jpg")
```

<br></br>
<font size="4">The PSSM has 20 columns. Therefore, the PSSM consists of 20 types of discrete signals (L-length). Therefore, we used the level 4 DWT as mentioned above to analyze these discrete signals from PSSM (each column) and to extract the PSSM-DWT feature vector from the PSSM, which is a feature vector of length 80 [@wang2019crystalm].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 as<-dwt_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)
```
<br></br>

# Disulfide_PSSM
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
table into training and testing data and predict the desired disulfide bonds between cysteines.Figure 14 shows a schematic of this process [@mapes2019residue]:</font>
<br></br>
<br></br>

```{r disulfid,echo=FALSE,fig.cap="Figure 14: The process of extracting disulfide-PSSM feature from the PSSM",out.width = '70%'}
knitr::include_graphics("figures/disulfid.jpg")
```

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
  aq<-disulfid(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(aq[,1:50])
```
<br></br>

# DP-PSSM
<font size="4">The extraction of this feature is obtained using the following equations from the PSSM:</font>


\begin{equation}
\begin{aligned}
P_{DP-PSSM}^{\alpha}&=[T',G']=[p_1,p_2,...,p_{40+40\times{\alpha}}]\\
T'&=[\bar{T}_1^P,\bar{T}_1^N,\bar{T}_2^P,\bar{T}_2^N,...,\bar{T}_{20}^P,\bar{T}_{20}^N]
\end{aligned}
\end{equation}

<br></br>

\begin{equation}
\left\{\begin{array}{ll}\bar{T}_j^P=\frac{1}{NP_j}\sum T_{i,j} &  ,if \ T_{i,j}\geq 0\\ 
\bar{T}_j^N=\frac{1}{NN_j}\sum T_{i,j} &  ,if \  T_{i,j}< 0
\end{array}\right.
\end{equation}

<br></br>

\begin{equation}
\begin{aligned}
G'&=[G_1,G_2,...,G_{20}]\\
G_j&=[\bar{\Delta}_{1,j}^P,\bar{\Delta}_{1,j}^N,\bar{\Delta}_{2,j}^P,\bar{\Delta}_{2,j}^N,...,\bar{\Delta}_{\alpha,j}^P,\bar{\Delta}_{\alpha,j}^N]
\end{aligned}
\end{equation}

<br></br>

\begin{equation}
\left\{\begin{array}{ll}\bar{\Delta}_{k,j}^P=\frac{1}{NDP_j}\sum [T_{i,j}-T_{i+k,j}]^2 &  ,if \ T_{i,j}-T_{i+k,j}\geq 0\\ 
\bar{\Delta}_{k,j}^N=\frac{-1}{NDN_j}\sum [T_{i,j}-T_{i+k,j}]^2 &  ,if \  T_{i,j}-T_{i+k,j}< 0
\end{array}\right.\\
0<k\leq{\alpha}
\end{equation}

<br></br>

<font size="4">In the above equations, $T_{i,j}$ represents the value on the i-th row and j-th column of the normalized PSSM, which is denoted by $M_T$. This matrix is constructed from the PSSM using the following equations:</font>

<br></br>

\begin{equation}
mean_i=\frac{1}{20}\sum_{i=1}^{20}E_{i,k}\\
STD_i=\sqrt{\frac{\sum_{u=1}^{20}[E_{i,u}-mean_i]^2}{20}}\\
T_{i,j}=\frac{E_{i,j}-mean_i}{STD_i}
\end{equation}

<font size="4">In above equations $\overline{T}_j^P$ represents the mean of the positive values of $\{T_{i,j}|i=1,2,...,L\}$ and $\overline T_j^N$ represents the mean of the negative values of the above set, which in fact the above set represents the j-th column of the matrix $M_T$ the expression $NP_j$ indicates the number of positive values of set $\{T_{i,j}|i=1,2,...,L\}$ and $NN_j$ is related to the number of negative values of the mentioned set. It is clear that this feature vector arises from the connection of two vectors $G'$,$T'$. According to the equations, it is clear that the length of the first feature vector is 40 and the length of the second feature vector is $\alpha\times 40$, which by selecting 2 in the used article, a feature vector of length 120 is created from the PSSM [@juan2009predicting].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 ss<-DP_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ss, n = 50)
```
<br></br>

# DFMCA-PSSM
<font size="4">In this feature, each of the columns of the PSSM is considered as a non-static time series. Assuming that $\{x_i\}$,$\{y_i\}$ for i=1,2,...,L represent two different columns of the PSSM, then two cumulative time series X,Y of these two columns are obtained according to the following equations:</font>

<br></br>

\begin{equation}
\left\{\begin{array}{ll}X_k=\sum_{i=1}^K x_i & k=1,2,...,L \\ 
Y_k=\sum_{i=1}^k y_i & k=1,2,...,L 
\end{array}\right.
\end{equation}

<font size="4">Now, using these two series, two backward moving average are obtained according to the following equations:</font>

\begin{equation}
\left\{\begin{array}{ll}\tilde{X}_{k,s}=\frac{1}{s}\sum_{i=-(s-1)}^0 X_{(k-i)} \\ 
\tilde{Y}_{k,s}=\frac{1}{s}\sum_{i=-(s-1)}^0 Y_{(k-i)} 
\end{array}\right.\\
1<s\leq L
\end{equation}

<font size="4">Finally, each element of the DFMCA feature vector is obtained using the above equations and the following formula:</font>

\begin{equation}
f_{DFMCA}^2(s)=\frac{1}{L-s+1}\sum_{k=1}^{L-s+1}(X_k-\tilde{X}_{k,s})(Y_k-\tilde{Y}_{k,s})
\end{equation}

<br></br>

<font size="4">According to the above equation, it is clear that each element of this feature vector is obtained by using two different columns of the PSSM, and since we have 20 different columns and the order of the columns does not matter, the length of the obtained feature vector will be equal to $\binom{20}{2}=\frac{20\times 19}{2}=190$ [@liang2018accurate]</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
  as<-DFMCA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7)
head(as, n = 50)
```
<br></br>

# grey_pssm_pseAAC

<font size="4">This function produces a feature vector of length 100. The first 20 components of this vector are the same as the normalized frequency of 20 standard amino acids in the protein. The second 20 components of this vector are the average of the 20 columns of the PSSM corresponding to the protein, and the grey system model method is used to define the next 60 components [@kabir2018improving]. If we show this feature vector with length 100 as follows:</font>

\begin{equation}
V=(\psi_1,\psi_2,...,\psi_{100})
\end{equation}

<font size="4"> then the first 20 components are as follows:</font>

\begin{equation}
\psi_i=f_i \quad (i=1,2,...,20)
\end{equation}

<font size="4">Where $f_i$ is the normalized frequency of type i amino acids of the 20 standard amino acids in the protein chain. If we denote the entries of the PSSM by $p_{i,j}$, then the next 20 components of this feature vector are obtained according to the following equation:</font>

\begin{equation}
\psi_{j+20}=\alpha_j \quad (j=1,2,...,20),\\
\alpha_j=\frac{1}{L}\sum_{i=1}^L p_{i,j}
\end{equation}

<font size="4"> the next 60 components are obtained by following equations:</font>

\begin{equation}
\psi_{j+40}=\delta_j \quad (j=1,2,...,60)
\end{equation}

<font size="4"> in this equation $\delta_j$'s are obtained as follows:</font>

\begin{equation}
\left\{\begin{array}{ll} \delta_{3j-2}=f_ja_1^j\\
\delta_{3j-1}=f_ja_2^j \quad (j=1,2,...,20) \\
\delta_{3j}=f_jb^j
\end{array}\right.
\end{equation}

<font size="4"> where $f_j$'s are as in above and $a_1^j$,\ $a_2^j$,\ $b^j$ are obtained as follows:</font>

\begin{equation}
\begin{bmatrix}a_1^j\\a_2^j\\b^j\end{bmatrix}=(B_j^TB_j)^{-1}B_j^TU_j \quad (j=1,2,...,20)
\end{equation}

<font size="4">where: </font>

\begin{equation}
B_j=\begin{bmatrix}
-p_{2,j}&-(p_{1,j}+0.5p_{2,j})&1\\
-p_{3,j}&-(\sum_{i=1}^2 p_{i,j}+0.5p_{3,j})&1\\
\vdots&\vdots&\vdots\\
-p_{L,j}&-(\sum_{i=1}^{L-1} p_{i,j}+0.5p_{L,j})&1
\end{bmatrix}
\end{equation}

<font size="4"> and:</font>

\begin{equation}
U_j=\begin{bmatrix}
p_{2,j}-p_{1,j}\\
p_{3,j}-p_{2,j}\\
\vdots\\
p_{L,j}-p_{L-1,j}
\end{bmatrix}
\end{equation}

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 as<-grey_pssm_pseAAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(as, n = 50)
```
<br></br>

# Smoothed_pssm 

<font size="4">This feature has been used to predict RNA binding sites in proteins, and therefore a specific feature vector is generated for each residue. To generate this feature vector, a matrix called smoothed_pssm is first created from the PSSM using a parameter ws called the smooth window size, which usually has a default value of 7 and can Change between 3, 5, 7, 9 and 11. The i-th row in the smoothed_PSSM is obtained by summing ws row vectors around the i-th row, which of course corresponds to the i-th row in the protein. If we show this problem with a mathematical equation, then we will have:</font>

\begin{equation}
V_{smoothed_i}=V_{i-\frac{(ws-1)}{2}}+...+V_i+...+V_{i+\frac{(ws-1)}{2}}
\end{equation}

<font size="4">In this regard, $V_i$'s represent the row vectors of the PSSM. To obtain the first and last rows of the smoothed_PSSM corresponding to the N-terminal and C-terminal of the protein, zero vectors are added to the beginning and end of the PSSM. The Figure 15 represents these processes schematically.</font>

<br></br>

```{r smooth,echo=FALSE,fig.cap="Figure 15: process of smoothed-PSSM generation, (A) represents the PSSM and (B) represents the smoothed-PSSM",out.width = '70%'}
knitr::include_graphics("figures/smoothed.jpg")
```

<br></br>

<font size="4">Now, using another parameter w called the slider window size, which its default value is 11 and can change in the interval (3,41)(ste=2), for the resid $\alpha_i$, the feature vector obtained from the smoothed-PSSM will be in the form as follows:</font>

\begin{equation}
(V_{smoothed_{i-\frac{(w-1)}{2}}},...,V_{smoothed_i},...,V_{smoothed_{i+\frac{(w-1)}{2}}})
\end{equation}

<font size="4">Here, as in the previous case, if the residue in question is the first or last residue of the protein, number of $(\frac{w-1}{2})$ zero vectors are added to the beginning or end of the smoothed matrix to obtain the feature vector. Therefore, the parameter w will determine the length of the feature vector per residue. If its value is 11, the length of the obtained feature vector will be equal to 220 = 11 * 20. eventually, the feature vector values are normalized between -1 and 1 [@cheng2008predicting].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 w<-smoothed_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7,11,c(2,3,8,9))
head(w[,1:50], n = 50)
```
<br></br>

# Kiderafactor
<font size="4">For producing this feature vector similar to the smoothed-PSSM feature, firstly PSSM is smoothed by appending zero vectors to its head and tail and utilizing sliding window of size odd, then this smoothed PSSM is condensed by the Kidera factors to produce feature vector for each residue [@fang2015condensing].</font>

<br></br>

<font size="4"> usage of this feature in PSSMCOOL package:</font>

```{r}
  w<-kiderafactor(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),c(2,3,8,9))
head(w[,1:50], n = 50)
```
<br></br>

# MBMGACPSSM
<font size="4">In this feature three different autocorrelation descriptors based on PSSM are adopted, which include: normalized Moreau-Broto autocorrelation, Moran autocorrelation and Geary autocorrelation
descriptors.Autocorrelation descriptor is a powerful statistical tool and defined based on the distribution of amino acid properties along the sequence, which measures the correlation between two
residues separated by a distance of d in terms of their evolution scores [@liang2015prediction2].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 w<-MBMGACPSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
<br></br>

# LPC-PSSM
<font size="4">This feature uses Linear predictive coding algorithm for each column of PSSM. So for
producing this feature vector "lpc" function from "phontools" R-package is used which produces a 14-dimensional vector for each column, since PSSM has 20 column eventually it will be obtained a 20*14=280 dimensional feature vector for each PSSM [@li2014pssp].</font> 

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 w<-LPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
<br></br>

# PSSM400
<font size="4">To generate this feature vector, for each of the standard amino acids, we find the positions containing that amino acid in the protein and separate the corresponding rows in the PSSM, to get a submatrix. Now, for the generated matrix, we calculate the average of its columns, and therefore, for each amino acid, a vector of length 20 is obtained. Finally, by putting these 20 vectors together, a feature vector of length 400 for each protein can be obtained. For example figure 16 shows the PSSM rows corresponding to amino acid S [@nanni2014empirical].</font>


```{r pssm400,echo=FALSE,fig.cap="Figure 16: process of extracting PSSM400 feature vector, which for amino acid S, represents the corresponding rows in PSSM",out.width = '70%'}
knitr::include_graphics("figures/pssm400.jpg")
```

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 q<-pssm400(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
head(q, n = 50)
```
<br></br>

# PSSM-BLOCK
<font size="4"> In this feature at first PSSM is divided to Blocks based on Number N which user imports.
Then for each Block the mean of columns is computed to get 20-dimensional vector, eventually by
appending these vectors to each other final feature vector is obtained [@an2017computational].</font> 

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 as<-PSSMBLOCK(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),5)
head(as, n = 50)
```
<br></br>

# PSSM-SD
<font size="4">To generate this feature vector, at first for the column j, the sum of the total numbers in this column is calculated and denoted by $T_j$. Then, starting from the first row of this column, the numbers are added one by one together to reach a number less than or equal to 25 percent of $T_j$. Now the number of components used to calculate this sum is denoted by $I_j^1$ and stored. Now, starting from the first row of this column again, the numbers are added one by one together to reach a number less than or equal to half of $T_j$ (50%), then we show the number of components to calculate this sum with $I_j^2$ and store it. In the same way for column j, starting from the last row of this column, we start adding each elements together to reach a number less than or equal to 25% of the number $T_j$, and denote the number of these components by $I_j^3$. In the next step, starting from the last row with summing of each element in this column to reach a number less than or equal to 50% of $T_j$, the number $I_j^4$ is also obtained. Therefore, 4 numbers are obtained for each column, and since the PSSM has 20 columns, for each protein a feature vector of length 80 is obtained [@dehzangi2014segmentation]. Figure 17 shows these steps schematically.</font>

<br></br>

```{r pssmsd,echo=FALSE,fig.cap="Figure 17: process of extracting PSSM-SD feature vector values for column j",out.width = '70%'}
knitr::include_graphics("figures/pssmsd.jpg")
```

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 ww<-PSSM_SD(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(ww, n = 50)
```
<br></br>

# PSSM-SEG
<font size="4"> This feature, similar to the previous feature, divides each column into four parts and calculates the values for each column. Then, using the following equations, it calculates the values of Segmented Auto Covariance Features. The final feature vector length will be of length 100.</font>

\begin{equation}
PSSM-Seg_{n,j}=\frac{1}{(I_j^n-m)}\sum_{i=1}^{I_j^n-m}(P_{i,j}-P_{ave,j})(P_{(i+m),j}-P_{ave,j})\\
(n=1,2,3,4 \ and \ j=1,...,20 \ and \ m\in \{1,2,...,11\})
\end{equation}

<font size="4">In the above equation, $P_{ave,j}$ represents the mean of column j in the PSSM and the number m is somehow a distance factor for each segment. Using the above equation, the feature of length 80 is obtained. Now the feature vector PSSM_AC is calculated using the previous factor m with length 20 and is added to the previous vector to get the final feature vector of length 100 [@dehzangi2014segmentation].</font>

\begin{equation}
PSSM-AC_{m,j}=\frac{1}{(L-m)}\sum_{i=1}^{L-m}(P_{i,j}-P_{ave,j})(P_{(i+m),j}-P_{ave,j})
\end{equation}

<font size="4"> L Represents the total length of the protein. </font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 q<-pssm_seg(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),3)
head(q, n = 50)
```
<br></br>

# SOMA-PSSM
<font size="4"> This feature also considers each of the columns of the PSSM as a time series. If L represents the length of the protein, then the column j of the matrix can be thought of $y(i), \ i=1,2,...,L$ as a time series. The SOMA algorithm is implemented in two steps using the following equations on the PSSM:
- First, the moving average $\overline{y_n}(i)$ for the time series $y(i)$ is calculated according to the following equation:</font>

\begin{equation}
\overline{y_n}(i)=\frac{1}{n}\sum_{k=0}^{n-1}y(i-k)
\end{equation}

<font size="4">Where n is the size of the moving average window, and if it tends to zero, the moving average will tend to original series, in other words: if ${n \to 0}$ then ${\overline{y_n}(i)\to y(i)}$. Next, for a moving average window size n which $2\leq n<L$, the second-order difference of the time series $y(i)$ with respect to the moving average $\overline{y_n}(i)$ is defined according to the following equation:</font>

\begin{equation}
\sigma_{MA}^2=\frac{1}{L-n}\sum_{i=n}^L[y(i)-\overline{y_n}(i)]^2
\end{equation}

<font size="4"> The number n must be smaller than the length of the smallest protein in the database under study. In the paper used by this algorithm, the length of the smallest protein is 10 and therefore the number n will vary from 2 to 9, so according to above Equation, by putting the numbers $\sigma_{MA}^2$ next to each other, 8 numbers are obtained for each column, and therefore the final feature vector will be of length 160 [@liang2017predict].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 w<-SOMA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 50)
```
<br></br>

# SVD-PSSM
<font size="4"> Singular value decomposition is a general-purpose matrix factorization approach
that has many useful applications in signal processing and statistics. In this feature SVD is
applied to a matrix representation of a protein aimed to reduce its dimensionality.
Given an input matrix Mat with dimensions N*M SVD is used to calculate its factorization
of the form: $Mat=U\Sigma V$ where $\Sigma$ is a diagonal matrix whose diagonal
entries are known as the singular values of Mat. The resulting descriptor is the ordered set of singular values: $SVD\in\mathcal{R}^L$ where L=min(M,N). since the PSSM has 20 columns, the final feature vector would be of length 20 [@nanni2014empirical].</font>

<br></br>

<font size="4"> Usage of this feature in PSSMCOOL package:</font>

```{r}
 w<-SVD_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
head(w, n = 20)
```

<br></br>

# References

<br></br>
