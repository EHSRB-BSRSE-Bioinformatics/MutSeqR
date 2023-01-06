

####################################################################################
# Generalized Linear Mixed Models
####################################################################################
#Library
library(lme4)
library(HLMdiag)
library(ggplot2)
library(Cairo)
library(robustlmm)
library(DHARMa)
library(car)
library(multcomp)
library(doBy)
library(writexl)

###################################################################################
####################################################################################
#Data
####################################################################################
#in case you need to add in contig info
genic_regions <- read.delim("~/GitHub/Duplex-Sequencing/per mouse/input files/DS_targets.txt")
sampledata <- read.delim("~/GitHub/Duplex-Sequencing/per mouse/input files/sampledata.txt", sep="\t", header=T, row.names=1)
sampledata$Dose <- factor(sampledata$Dose)
sampledata$sample <- row.names(sampledata)

dat <- read_excel("Mut_count_by_target_clonal.xlsx")

#add dose column
dat <- dplyr::left_join(dat, sampledata)
names(dat)[2]<-"description"
#add contig column and location relative to genes
dat <- dplyr::left_join(dat, genic_regions, by=c("description"), mode="inner")

as.factor<-dat$Region

#Indicate dose group of each sample
dat$dose <- "D0"
dat$dose[dat$sample %in% c("DNA02446", "DNA02451", "DNA02452", "DNA02453", "DNA02454", "DNA02455")] <- "D1"
dat$dose[dat$sample %in% c("DNA02456", "DNA02457", "DNA02458", "DNA02459", "DNA02460", "DNA02461")] <- "D2"
dat$dose[dat$sample %in% c("DNA02462", "DNA02463", "DNA02464", "DNA02465", "DNA02466", "DNA02467")] <- "D3"

#Indicate transcription status of targets
dat$Type <- "intergenic"
dat$Type[dat$Region %in% c("chr1.2", "chr3", "chr6", "chr7", "chr9", "chr12", "chr13", "chr15", "chr19")] <- "genic"
dat %>% mutate(mut_freq = (mut_depth)/total_depth)
head(dat)

#    sample mut_depth total_depth mut_freq Region dose
#1 dna00973         0    22443807 0.00e+00   chr1   D0
#2 dna00974         1    30189507 3.31e-08   chr1   D0
#3 dna00975         6    25372720 2.36e-07   chr1   D0
#4 dna00976         2    26740225 7.48e-08   chr1   D0
#5 dna00977         1    24928760 4.01e-08   chr1   D0
#6 dna00978         1    25082020 3.99e-08   chr1   D0
length(table(dat$sample))

m <- glmer(cbind(mut_depth,total_depth) ~ dose*Region + (1|sample), family = "binomial", data = dat, control = glmerControl(check.conv.grad = .makeCC("warning", tol = 3e-3, relTol = NULL)))

summary(m)

Anova(m)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: cbind(mut_depth, total_depth)
#              Chisq Df Pr(>Chisq)    
#dose        1310.92  3  < 2.2e-16 ***
#Region      1376.86 19  < 2.2e-16 ***
#dose:Region  124.13 57  6.838e-07 ***
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

dat$Resid <- residuals(m)
dat[abs(residuals(m)) == max(abs(residuals(m))),]


X11()
par(las = 2)

boxplot(dat$mut_freq ~ as.factor(paste(dat$dose, dat$Region)), col = "yellow", xlab = "")

hist(residuals(m), main = "Residuals")
qqnorm(residuals(m))
abline(h = -3, col = "red", lwd = 2)
abline(h = 3, col = "red", lwd = 2)

#####Point Estimates By Dose and By Target####################################
nCoef <- length(names(coef(m)$sample))
a <- c(1, rep(0, nCoef - 1))
lambda <- NULL

d <- paste("dose", unique(dat$dose), sep = "")
r <- paste("Region", unique(dat$Region), sep = "")
indx <- 0

for(k in 1:length(d)){
  for(j in 1:length(r)){
    b <- a
    #Dose
    flag <- names(coef(m)$sample) == d[k]
    if(length(b[flag]) > 0){
      b[flag] <- 1
    }
    #Region
    flag <- names(coef(m)$sample) == r[j]
    if(length(b[flag]) > 0){
      b[flag] <- 1
    }
    #Interaction
    flag <- names(coef(m)$sample) == paste(d[k], r[j], sep = ":")
    if(length(b[flag]) > 0){
      b[flag] <- 1
    }
    
    lambda <- rbind(lambda, b)
    indx <- indx + 1
    row.names(lambda)[indx] <- paste(d[k], r[j], sep = ":")
  }
}
colnames(lambda) <- names(coef(m)$sample) 


A <-esticon(m,lambda)

A <- as.data.frame(A)
A$estimate <- exp(A$estimate)
delta <- A$estimate^2
A$lwr <- exp(A$lwr)
A$upr <- exp(A$upr)
A$std.error <- sqrt(delta*A$std.error^2)
A <- A[,-c(3,4,5,6)]
colnames(A) <- c("Estimate", "Std.Err", "Lower", "Upper")

write_xlsx(A, "MF_estimates_by_dose_target_2.xlsx")

#######################################################
##########Comparisons between dose by target###########
nCoef <- length(names(coef(m)$sample))
a <- rep(0, nCoef)
lambda <- NULL

d <- paste("dose", unique(dat$dose), sep = "")
d <- d[-1]

r <- paste("Region", unique(dat$Region), sep = "")
indx <- 0

for(j in 1:length(r)){
  for(k in 1:length(d)){
    b <- a
    #Dose
    flag <- names(coef(m)$sample) == d[k]
    if(length(b[flag]) > 0){
      b[flag] <- 1
    }
    
    #Interaction
    flag <- names(coef(m)$sample) == paste(d[k], r[j], sep = ":")
    if(length(b[flag]) > 0){
      b[flag] <- 1
    }
    
    lambda <- rbind(lambda, b)
    indx <- indx + 1
    row.names(lambda)[indx] <- paste(d[k], "vsD0", r[j], sep = ":")
  }
}
colnames(lambda) <- names(coef(m)$sample) 

B <-esticon(m,lambda)

B <- as.data.frame(B)
B$estimate <- exp(B$estimate)
delta <- B$estimate^2
B$lwr <- exp(B$lwr)
B$upr <- exp(B$upr)
B$std.error <- sqrt(delta*B$std.error^2)
B <- B[,-5]
colnames(B) <- c("Estimate", "Std.Err", "Obs.T", "p.value", "df", "Lower", "Upper")

B

B$adjP <- 0

#Holm-Sidak
for(j in 1:length(r)){
  flag <- unlist(lapply(rownames(B), function(x) strsplit(x, ":")[[1]][3])) == r[j]
  B$adjP <- my.holm.sidak(B$p.value)
}

B <- B[,c(1:3,5,4,8,6:7)] 
write_xlsx(B, "Comparison_by_dose_target_2.xlsx")

######################################################
#######################################################
# intergenic vs genic 
#Point Estiamtes
ni <- 11
ng <-  9

#genic
g <- c("Regionchr1.2", "Regionchr3", "Regionchr6", "Regionchr7", "Regionchr9", "Regionchr12",
       "Regionchr13", "Regionchr15", "Regionchr19")

a <- c(1, rep(0, nCoef - 1))

lambda <- NULL

d <- paste("dose", unique(dat$dose), sep = "")
r <- paste("Region", unique(dat$Region), sep = "")

for(k in 1:length(d)){
  b <- a
  #Dose
  flag <- names(coef(m)$sample) == d[k]
  if(length(b[flag]) > 0){
    b[flag] <- 1
  }
  #Region
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == r[j] & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ng
    }
  }
  
  #Interaction
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == paste(d[k], r[j], sep = ":") & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ng
    }
  }
  
  lambda <- rbind(lambda, b)
  
}

#intergenic
a <- c(1, rep(0, nCoef - 1))
g <- c("Regionchr1", "Regionchr2", "Regionchr4", "Regionchr5", "Regionchr8", "Regionchr10",
       "Regionchr11", "Regionchr14", "Regionchr16", "Regionchr17", "Regionchr18")

for(k in 1:length(d)){
  b <- a
  #Dose
  flag <- names(coef(m)$sample) == d[k]
  if(length(b[flag]) > 0){
    b[flag] <- 1
  }
  #Region
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == r[j] & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ni
    }
  }
  
  #Interaction
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == paste(d[k], r[j], sep = ":") & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ni
    }
  }
  
  lambda <- rbind(lambda, b)
  
}

A <-esticon(m,lambda)

A <- as.data.frame(A)
A$estimate <- exp(A$estimate)
delta <- A$estimate^2
A$lwr <- exp(A$lwr)
A$upr <- exp(A$upr)
A$std.error <- sqrt(delta*A$std.error^2)
A <- A[,-c(3,4,5,6)]
colnames(A) <- c("Estimate", "Std.Err", "Lower", "Upper")
A

#Comparison - Intergenic vs Genic

d0 <- lambda[5,] - lambda[1,]
d1 <- lambda[6,] - lambda[2,]
d2 <- lambda[7,] - lambda[3,]
d3 <- lambda[8,] - lambda[4,]
lambda <- rbind(d0, d1, d2, d3)

B <-esticon(m,lambda)

B <- as.data.frame(B)
B$estimate <- exp(B$estimate)
delta <- B$estimate^2
B$lwr <- exp(B$lwr)
B$upr <- exp(B$upr)
B$std.error <- sqrt(delta*B$std.error^2)
B <- B[,-5]
colnames(B) <- c("Estimate", "Std.Err", "Obs.T", "p.value", "df", "Lower", "Upper")

B

my.holm.sidak(B$p.value)

######################################
#In Heterochromatin - No
#Point Estiamtes
ni <- 7
ng <- 13

#genic
g <- c("Regionchr1", "Regionchr3", "Regionchr4", "Regionchr5", "Regionchr6", "Regionchr7", "Regionchr8",
       "Regionchr9", "Regionchr10", "Regionchr13", "Regionchr15", "Regionchr17", "Regionchr18")

a <- c(1, rep(0, nCoef - 1))

lambda <- NULL

d <- paste("dose", unique(dat$dose), sep = "")
r <- paste("Region", unique(dat$Region), sep = "")

for(k in 1:length(d)){
  b <- a
  #Dose
  flag <- names(coef(m)$sample) == d[k]
  if(length(b[flag]) > 0){
    b[flag] <- 1
  }
  #Region
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == r[j] & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ng
    }
  }
  
  #Interaction
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == paste(d[k], r[j], sep = ":") & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ng
    }
  }
  
  lambda <- rbind(lambda, b)
  
}

#Yes
a <- c(1, rep(0, nCoef - 1))
g <- c("Regionchr1.2", "Regionchr2", "Regionchr11", "Regionchr12", "Regionchr14", "Regionchr16", "Regionchr19")

for(k in 1:length(d)){
  b <- a
  #Dose
  flag <- names(coef(m)$sample) == d[k]
  if(length(b[flag]) > 0){
    b[flag] <- 1
  }
  #Region
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == r[j] & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ni
    }
  }
  
  #Interaction
  for(j in 1:length(r)){
    flag <- names(coef(m)$sample) == paste(d[k], r[j], sep = ":") & r[j] %in% g
    if(length(b[flag]) > 0){
      b[flag] <- 1/ni
    }
  }
  
  lambda <- rbind(lambda, b)
  
}

A <-esticon(m,lambda)

A <- as.data.frame(A)
A$estimate <- exp(A$estimate)
delta <- A$estimate^2
A$lwr <- exp(A$lwr)
A$upr <- exp(A$upr)
A$std.error <- sqrt(delta*A$std.error^2)
A <- A[,-c(3,4,5,6)]
colnames(A) <- c("Estimate", "Std.Err", "Lower", "Upper")
A
#Comparison - Heterochromatin vs Euchromatin

d0 <- lambda[5,] - lambda[1,]
d1 <- lambda[6,] - lambda[2,]
d2 <- lambda[7,] - lambda[3,]
d3 <- lambda[8,] - lambda[4,]
lambda <- rbind(d0, d1, d2, d3)

B <-esticon(m,lambda)

B <- as.data.frame(B)
B$estimate <- exp(B$estimate)
delta <- B$estimate^2
B$lwr <- exp(B$lwr)
B$upr <- exp(B$upr)
B$std.error <- sqrt(delta*B$std.error^2)
B <- B[,-5]
colnames(B) <- c("Estimate", "Std.Err", "Obs.T", "p.value", "df", "Lower", "Upper")

B


my.holm.sidak(B$p.value)
write_xlsx(B, "comp_Het_eu.xlsx")
#########