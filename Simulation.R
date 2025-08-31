rm(list = ls())
library(Rfast)
library(foreach)
library(doParallel)
library(doRNG)
library(irlba)
library(RcppNumerical)
source("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Simulation\\Simfunctions.R")

# param
rep <- 500
par <- read.csv("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Simulation\\SimModels.csv", header = T, colClasses=rep("numeric", 4))

for (j in 1:6) {
  n   <- par$n[j]
  T0   <- par$T0[j]
  model <- funlist[[par$hdmodel[j]]]
  
  
  
  ###############################################
  # calculate theta0
  p1 <- integrate(px, 0, 1)$value
  integrand <- function(alpha) {
    return(mu0(alpha) * px(alpha))
  }
  theta0 <- (integrate(integrand, 0, 1)$value) / p1
  ###############################################
  
  # 开启并行
  cl <- makeCluster(detectCores() - 2)  
  registerDoParallel(cl)
  
  output <- foreach (i = 1:rep, 
                     .options.RNG=1234, 
                     .packages=c('Rfast')) %dorng% {
                       sim(n=n, T0=T0, hdmodel=model, theta0=theta0)
                     }
  
  stopCluster(cl)
  
  # transfer to data.frame
  output_df <- do.call(rbind, lapply(output, as.data.frame))
  
  # save the results
  write.table(output_df, 
              paste0("output_deaner_par", j, ".txt"), 
              sep = ",", row.names = F, col.names = T)
}