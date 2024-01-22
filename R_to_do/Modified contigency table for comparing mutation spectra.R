##################################################################
# Functions
# Run all the functions so that they are in your global environment
my.holm.sidak <- function(P){
  m <- length(P)
  if(m > 1){
    Psort <- matrix(c(P, 1:m), 2, m, byrow = TRUE)
    Psort <- Psort[, order(Psort[1, ])]
    for (i in 1:m) {
      adjust <- m + 1 - i
      Psort[1, i] <- pmin(1, (1 - ((1 - Psort[1, i])^adjust)))
    }
    Psort <- Psort[, order(Psort[2, ])]
    P.adjust <- Psort[1, ]
    P.adjust
  } else{
    #No adjustment as there is only one comparison
    P.adjust <- P
    P.adjust
  }
}

#G2 Statistic - Likelihood Ratio Statistic
G2 <- function(x, monte.carlo = FALSE, n.sim = 10000, seed = 1234){
  N <- sum(x)
  r <- apply(x, 1, sum)
  c <- apply(x, 2, sum)
  e <- r %*% t(c)/N
  
  G2 <- 0
  for(k in 1:ncol(x)){
    flag <- x[,k] > 0
    G2 <- G2 + t(x[flag,k]) %*% log(x[flag,k]/e[flag,k])	
  }
  G2 <- 2*G2
  
  R <- nrow(x)-1
  df <- R * (ncol(x)-1)
  
  if(monte.carlo == FALSE){
    if(N/R > 20){
      p.value <- 1-pchisq(G2, df)
    } else{
      p.value <- 1-pf(G2/R, R, N-df)
    }
  } else {
    #Monte Carlo
    #Generate random rxc tables
    set.seed(seed)
    r <- apply(x, 1, sum)
    c <- apply(x, 2, sum)
    rtbl <- r2dtable(1, r, c)
    
    ref.dist <- rep(0, n.sim)
    for(k in 1:length(ref.dist)){
      x <- r2dtable(1, r, c)[[1]]
      N <- sum(x)
      r <- apply(x, 1, sum)
      c <- apply(x, 2, sum)
      e <- r %*% t(c)/N
      
      G2.t <- 0
      for(j in 1:ncol(x)){
        flag <- x[,j] > 0
        G2.t <- G2.t + t(x[flag,j]) %*% log(x[flag,j]/e[flag,j])	
      }
      
      ref.dist[k] <- 2*G2.t
    } 
    flag <- ref.dist >= G2[1,1]
    p.value <- length(ref.dist[flag])/10000
  }
  
  
  data.frame(G2 = G2, p.value = p.value)
}


G2ab <- function(x, R, T, N, monte.carlo = FALSE, n.sim = 10000, seed = 1234){
  
  flag <- x > 0
  S1 <- t(x[flag]) %*% log(x[flag])
  S2 <- sum(apply(x, 1, sum) * log(apply(x, 1, sum)))
  S3 <- sum(apply(x, 2, sum) * log(apply(x, 2, sum)))
  S4 <- sum(x) * log(sum(x))
  
  G2ab <- 2*(S1 - S2 - S3 + S4)
  
  df <- (R-1)*(T-1)
  
  if(monte.carlo == FALSE){
    
    if(N/R > 20){
      p.value <- 1-pchisq(G2ab, df)
    } else{
      p.value <- 1-pf(G2ab/df, df, N-df)
    }
  } else {
    #Monte Carlo
    #Generate random rxc tables
    set.seed(seed)
    r <- apply(x, 1, sum)
    c <- apply(x, 2, sum)
    rtbl <- r2dtable(1, r, c)
    
    ref.dist <- rep(0, n.sim)
    for(k in 1:length(ref.dist)){
      x <- r2dtable(1, r, c)[[1]]
      flag <- x > 0
      S1 <- t(x[flag]) %*% log(x[flag])
      S2 <- sum(apply(x, 1, sum) * log(apply(x, 1, sum)))
      S3 <- sum(apply(x, 2, sum) * log(apply(x, 2, sum)))
      S4 <- sum(x) * log(sum(x))
      
      ref.dist[k] <- 2*(S1 - S2 - S3 + S4)
      
    } 
    flag <- ref.dist >= G2ab[1,1]
    p.value <- length(ref.dist[flag])/10000
  }
  
  data.frame(G2ab = G2ab, p.value = p.value)
}

##################################################################
# Perform the comparisons
# Groups to compare are rows, mutation subtypes are columns

# {row,col}

#Compare the spectrum between germ cells (ST) and bone marrow (BM)
# load in your data
dat <- read.delim("~/Masters/PRC data/Germ Cells/Sequencing Output/Contract 1 Control and High dose/R_suptypes_compare_BM2.txt")
names(dat)

#Dose 0
ST0 <- as.numeric(dat[1,2:10]) # row 1, columns 2-10
BM0 <- as.numeric(dat[5,2:10]) # row 5, columns 2-10
STvsBM0 <- G2(cbind(ST0,BM0)) # compare
STvsBM0 # see results (but we need to run the holm-sidak correction on the pvalues first)
# G2      p.value
# 93.88126 1.110223e-16

#Dose 1
ST1 <- as.numeric(dat[2,2:10])
BM1<-as.numeric(dat[6,2:10])
STvsBM1 <- G2(cbind(ST1,BM1))
STvsBM1

#Dose 2
ST2 <- as.numeric(dat[3,2:10])
BM2<-as.numeric(dat[7,2:10])
STvsBM2 <- G2(cbind(ST2,BM2))
STvsBM2

#Dose 3
ST3 <- as.numeric(dat[4,2:10])
BM3<-as.numeric(dat[8,2:10])
STvsBM3 <- G2(cbind(ST3,BM3))
STvsBM3


# Perform the holm-sidak correction for multiple comparison on the p-values
my.holm.sidak(c(as.numeric(STvsBM0[2]), as.numeric(STvsBM1[2]), as.numeric(STvsBM2[2]),as.numeric(STvsBM3[2])))  
# 1.110223e-16 0.000000e+00 0.000000e+00 0.000000e+00

###################
# Compare ST to BM all doses
dat <- read.delim("~/Masters/PRC data/Germ Cells/Sequencing Output/Contract 1 Control and High dose/R_PRC_ST vs BM_subtypes.txt")

#Dose 0
ST0<- as.numeric(dat[1,2:10])
ST1<- as.numeric(dat[2,2:7])
ST2<- as.numeric(dat[3,2:7])
ST3<- as.numeric(dat[4,2:7])

BM0 <- as.numeric(dat[5,2:10])
BM1 <- as.numeric(dat[6,2:7])
BM2 <- as.numeric(dat[7,2:7])
BM3 <- as.numeric(dat[8,2:7])

STvsBM.0 <- G2(cbind(ST0,BM0))
STvsBM.1 <- G2(cbind(ST1,BM1))
STvsBM.2 <- G2(cbind(ST2,BM2))
STvsBM.3 <- G2(cbind(ST3,BM3))

STvsBM.0 
#G2      p.value
#1 93.88126 1.110223e-16

STvsBM.1
#G2 p.value
#1 124.6071       0

STvsBM.2
#G2 p.value
#1 147.3143       0

STvsBM.3
#G2 p.value
#1 112.6308       0


my.holm.sidak(c(as.numeric(STvsBM.0[2]), as.numeric(STvsBM.1[2]), as.numeric(STvsBM.2[2]),as.numeric(STvsBM.3[2])))




