###################################################
#### This file replicates the simulation results ##
############## in Feng (2024)  ####################
########### Author: Yingjie Feng ##################
########## Last updated: 04/18/2024 ################

rm(list = ls())
library(Rfast)
library(foreach)
library(doParallel)
library(doRNG)
library(irlba)
library(RcppNumerical)
source("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\Feng-2024_SimFuns.R")

# param
rep <- 500
par <- read.csv("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\Feng-2024_SimModels.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))

for (j in 1:6) {
n   <- par$n[j]
p   <- par$p[j]
model <- funlist[[par$hdmodel[j]]]
nlam   <- par$nlam[j]
K      <- n^(4/5)
Kseq   <- floor(K * c(0.5, 1, 1.5))
const  <- par$const[j]


###############################################
# calculate true value
p1 <- integrate(px, 0, 1)$value
integrand <- function(alpha) {
  return(mu0(alpha) * px(alpha))
}
theta0 <- (integrate(integrand, 0, 1)$value)/p1
#####################################################

#ptm <- proc.time()
## simulation
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

output <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('Rfast','irlba','RcppNumerical'),
                   .combine=rbind) %dorng% {
                     output <- sim(i, n, p, model, Kseq, nlam, theta0, const)
                     output   # (5*rep) by (length(Kseq)+1) matrix
                   }

stopCluster(cl)

#proc.time() - ptm
###################
write.table(output, paste("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\rawoutput_par", j, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)

}
