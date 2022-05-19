library(car)
library(multcomp)
library(doBy)
library(readxl)
##################################################################
#Reading in the data
##################################################################

PRC_Mutation_Frequency_R <- read_excel("Masters/PRC data/PRC Data Sheets/R sheets/PRC Mutation Frequency R.xlsx")
PRC_Mutation_Frequency_R <- read_excel("~/Masters/PRC data/PRC Data Sheets/R sheets/PRC Mutation Frequency R_Clones.xlsx")
View(PRC_Mutation_Frequency_R)
dat<-PRC_Mutation_Frequency_R
head(dat)
dat$dose <- as.factor(dat$dose)


##################################################################
# Generalized Linear Model
# quasibinomial in order to account for over dispersion
##################################################################
#for lacZ data#
m <- glm(cbind(mutants,plaques) ~ dose, family = "quasibinomial", data = dat)
dat$MF <- dat$mutants/dat$plaques

#For DS data#
m <- glm(cbind(mut_depth,total_depth) ~ dose, family = "quasibinomial", data = dat)
dat$MF <- dat$mut_depth/dat$total_depth


dat$Resid <- m$residuals

#removing obs with residuals greater than 4 in absolute value
print(dat[abs(dat$Resid) == max(abs(dat$Resid)),])
#  dose Animal.ID mutants   plaques           MF    Resid
#3    0  dna00975      92 555207445 1.657038e-07 0.316941

#CLONAL DATA
# sample   description mut_depth total_depth    mut_freq    lower_ci    upper_ci dose           MF Resid
#<chr>    <chr>           <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <fct>       <dbl> <dbl>
#  1 DNA02444 PRC 121           129   728699669 0.000000177 0.000000130 0.000000188 0     0.000000177 0.255

hist(dat$Resid, main = "Residuals", col = "yellow")
library(ggplot2)
qqnorm(dat$Resid)
abline(h = -3, col = "red", lwd = 2)
abline(h = 3, col = "red", lwd = 2)
qqline(dat$Resid)
##################################################################
# Contrast matrix for the point estimates
# You'll need to modify the matrix to you study, number of dose groups
# and the dose ordering
##################################################################
coef(m)
lambda <- rbind(
  c(1,   0,   0,   0),
  c(1,   1,   0,   0),
  c(1,   0,   1,   0),  
  c(1,   0,   0,   1))  

#estican Contrasts for glm. Computes weighted sums of the estimated regression parameters
A <-esticon(m,lambda)

###Wald Statistics
A <- as.data.frame(A)

#take the exponential function of the rows?
A$estimate <- exp(A$estimate)
delta <- A$estimate^2
A$lwr <- exp(A$lwr)
A$upr <- exp(A$upr)
A$std.error <- sqrt(delta*A$std.error^2)
#remove extra columns (statistic, p value, beta0, df)
A <- A[,-c(3,4,5,6)]
colnames(A) <- c("Estimate", "Std.Err", "Lower", "Upper")
rownames(A) <- c("Control", "D1", "D2", "D3")
A

#            Estimate      Std.Err        Lower        Upper
#Control 1.258248e-07 1.040774e-08 1.058844e-07 1.495205e-07
#D1      3.321608e-07 1.632949e-08 2.997865e-07 3.680314e-07
#D2      6.754085e-07 2.255404e-08 6.299628e-07 7.241327e-07
#D3      1.041332e-06 2.771075e-08 9.851031e-07 1.100770e-06


##################################################################
# Contrast matrix for the pairwise comparisons
# You'll need to modify the matrix to you study, number of dose groups
# and the dose ordering
##################################################################
lambda <- rbind(
  c(0,   1,   0,   0),
  c(0,   0,   1,   0),  
  c(0,   0,   0,   1))

B <-esticon(m,lambda)
B <- as.data.frame(B)

B$estimate <- exp(B$estimate)
delta <- B$estimate^2
B$lwr <- exp(B$lwr)
B$upr <- exp(B$upr)
B$std.error <- sqrt(delta*B$std.error^2)

B <- B[,-5]
colnames(B) <- c("Estimate", "Std.Err", "Obs.T", "p.value", "df", "Lower", "Upper")
rownames(B) <- c("D1.to.Cont", "D2.to.Cont", "D3.to.Cont")
B
write_xlsx(B, "MF compairsons lacZ.xlsx")
#Estimate   Std.Err    Obs.T      p.value df    Lower    Upper
#D1.to.Cont 2.639868 0.2540152 10.08836 2.728173e-09 20 2.159792 3.226655
#D2.to.Cont 5.367849 0.4788250 18.83837 3.375078e-14 20 4.456458 6.465628
#D3.to.Cont 8.276045 0.7191164 24.32194 2.220446e-16 20 6.904084 9.920639

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
my.holm.sidak(B$p.value)
