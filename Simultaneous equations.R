

rm(list = ls())
library(BayesianTools)

#Temperory plots
ysTP <- read.table("TP.csv", header = TRUE, sep = ",")
#Permanent plots
ysPP <- read.table("PP.csv", header = TRUE, sep = ",")
#Stem analysis data
ysSA <- read.table("SA.csv", header = TRUE, sep = ",")

#Restrict ranges for effective searching
para.prior <- read.csv(file = 'Prior.csv')

# Heavy tailed noraml distribution
Sivia_log <- function(diff, sd) {
  diff[which(abs(diff) <= 1e-6)] <- 1e-6
  R2 <- (diff / sd) ^ 2
  prob <- 1 / (sd * (pi * 2) ^ 0.5) * (1 - exp(-R2 / 2)) / R2
  log(prob)
}

#Likelihood
likelihood <- function(pValues) {
  # pValues<-para.prior$Guess
  # pValues<-para.prior$Max
  # SDI model
  pred_lnN_tp <- pValues[1] - pValues[2] * log(ysTP$D)
  pred_lnN_pp <- pValues[1] - pValues[2] * log(ysPP$D)
  ysTP$SDI <- ysTP$N * (ysTP$D / 20) ^ pValues[2]
  ysPP$SDI <- ysPP$N * (ysPP$D / 20) ^ pValues[2]
  # Height model
  pred_Hgc_tp <- pValues[3] * exp(-pValues[4] / ysTP$A)
  pred_Hgc_pp <- pValues[3] * exp(-pValues[4] / ysPP$A)
  pred_Hgc_sa <- pValues[3] * exp(-pValues[4] / ysSA$A)
  
  ysTP$SCI <- ysTP$H * exp(pValues[4] / ysTP$A - pValues[4] / 30)
  ysPP$SCI <- ysPP$H * exp(pValues[4] / ysPP$A - pValues[4] / 30)
  # Basal area model
  pred_G_tp <-
    pValues[5] * ysTP$SCI ^ pValues[6] * exp(-pValues[7] * (ysTP$SDI / 1000) ^
                                               (-pValues[8]) / ysTP$A)
  pred_G_pp <-
    pValues[5] * ysPP$SCI ^ pValues[6] * exp(-pValues[7] * (ysPP$SDI / 1000) ^
                                               (-pValues[8]) / ysPP$A)
  # Diameter model
  pred_D_tp <-
    pValues[9] * ysTP$SCI ^ pValues[10] * exp(-pValues[11] * (ysTP$SDI / 1000) ^
                                                (pValues[12]) / ysTP$A)
  pred_D_pp <-
    pValues[9] * ysPP$SCI ^ pValues[10] * exp(-pValues[11] * (ysPP$SDI / 1000) ^
                                                (pValues[12]) / ysPP$A)
  # likelihood
  loglikelihood <-
    sum(Sivia_log(
      diff = pred_lnN_tp - log(ysTP$N),
      sd = pValues[13] + pValues[14] * log(ysTP$N)
    )) +
    sum(Sivia_log(
      diff = pred_lnN_pp - log(ysPP$N),
      sd = pValues[15] + pValues[16] * log(ysPP$N)
    )) +
    sum(Sivia_log(
      diff = pred_Hgc_tp - ysTP$H,
      sd = pValues[17] + pValues[18] * log(ysTP$H)
    )) +
    sum(Sivia_log(
      diff = pred_Hgc_pp - ysPP$H,
      sd = pValues[19] + pValues[20] * log(ysPP$H)
    )) +
    sum(Sivia_log(
      diff = pred_Hgc_sa - ysSA$H,
      sd = pValues[21] + pValues[22] * log(ysSA$H)
    )) +
    sum(Sivia_log(
      diff = pred_G_tp - ysTP$G,
      sd = pValues[23] + pValues[24] * log(ysTP$G)
    )) +
    sum(Sivia_log(
      diff = pred_G_pp - ysPP$G,
      sd = pValues[25] + pValues[26] * log(ysPP$G)
    )) +
    sum(Sivia_log(
      diff = pred_D_tp - ysTP$D,
      sd = pValues[27] + pValues[28] * log(ysTP$D)
    )) +
    sum(Sivia_log(
      diff = pred_D_pp - ysPP$D,
      sd = pValues[29] + pValues[30] * log(ysPP$D)
    ))
  if (is.na(loglikelihood)) {
    loglikelihood <- -1e10
  }
  return(loglikelihood)
}

# Define the prior. Here uninform distriuction was defined to narrow the searching ranges
prior <-
  createUniformPrior(lower = para.prior$Min, upper = para.prior$Max)
# Create Bayesian Setup
BCSetup  <-
  createBayesianSetup(
    likelihood,
    prior,
    best = para.prior$Guess,
    names = as.character(para.prior$Name),
    parallel = F
  )

settings <- list(iterations = 3e6,
                 nrChains = 1,
                 thin = 10)
cali <- runMCMC(BCSetup, sampler = "DREAMzs", settings = settings)
save(cali, file = "cali.RData")
# Check convergence
# plot(cali)
gelmanDiagnostics(cali)
MAP(cali)
# save summary of parameters
n <- dim(getSample(cali, thin = 10))[1]
Psample <-
  # three chains, for each take the latter half
  getSample(cali, start = floor(n * 10 / 3 / 2), thin = floor(n * 10 / 3000)) # three chains, for each take the latter half
parameters <- data.frame(
  name = colnames(Psample),
  Mean = apply(Psample, 2, mean),
  SD = apply(Psample, 2, sd),
  MAP = MAP(cali)[[1]]
)
write.csv(parameters, file = 'parameters.csv')

#### variance based uncertainty partitioning
library(sensitivity)
Psample <-
  getSample(cali, start = floor(n * 10/6), thin = 10) 
dim(Psample)
X1<-Psample[1:5000,1:12]
X2<-Psample[5001:10000,1:12]

pred_H<-function(pValues){
  y <- pValues[,3] * exp(-pValues[,4] / 40)
  y
}

pred_G<-function(pValues){
   y <- pValues[,5] * 11 ^ pValues[,6] * exp(-pValues[,7] * (600 / 1000) ^
                                             (-pValues[,8]) / 40) 
   y
}

pred_D<-function(pValues){
  y <- pValues[,9] * 11 ^ pValues[,10] * exp(-pValues[,11] * (600 / 1000) ^
                                              (-pValues[,12]) / 40) 
  y
}

pred_N<-function(pValues){
  y1 <- pValues[,5] * 11 ^ pValues[,6] * exp(-pValues[,7] * (600 / 1000) ^
                                              (-pValues[,8]) / 40) 
  y2 <- pValues[,9] * 11 ^ pValues[,10] * exp(-pValues[,11] * (600 / 1000) ^
                                               (pValues[,12]) / 40) 
  y<-y1/(y2*pi/40000)
  y
}


sobol.H <- sobolSalt(model = pred_H, X1, X2, scheme="A", nboot = 1000000)
sobol.H

sobol.G <- sobolSalt(model = pred_G, X1, X2, scheme="A", nboot = 1000000)
sobol.G

sobol.D <- sobolSalt(model = pred_D, X1, X2, scheme="A", nboot = 100000)
sobol.D

sobol.N <- sobolSalt(model = pred_N, X1, X2, scheme="A", nboot = 100000)
sobol.N

sobol.T<-data.frame(H=sobol.H$T[,1],
           G=sobol.G$T[,1],
           D=sobol.D$T[,1],
           N=sobol.N$T[,1])
write.csv(sobol.T,file = 'Sobol total effect.csv')

sobol.F<-data.frame(H=sobol.H$S[,1],
                    G=sobol.G$S[,1],
                    D=sobol.D$S[,1],
                    N=sobol.N$S[,1])
write.csv(sobol.F,file = 'Sobol first order effect.csv')
