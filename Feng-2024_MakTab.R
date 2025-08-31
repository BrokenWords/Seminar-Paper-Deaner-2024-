#########################################################
##### This file produces Table 3 in the main paper#######
########################################################
library(Hmisc)
library(Rfast)
source("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\Feng-2024_SimFuns.R")

############################################
# Table used in the main paper
par <- read.csv("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\Feng-2024_SimModels.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))


# 初始化结果矩阵
results <- matrix(NA, nrow=6, ncol=3)
colnames(results) <- c("Median Abs. Dev.", "95% CI Coverage", "Median CI Length")
rownames(results) <- paste0("Setting", 1:6)

for (j in 1:6) {
  
  n   <- par$n[j]
  p   <- par$p[j]
  model <- funlist[[par$hdmodel[j]]]
  nlam <- par$nlam[j]
  K    <- n^(4/5)
  Kseq <- floor(K * c(0.5, 1, 1.5))
  ###############################################
  # calculate true value
  p1 <- integrate(px, 0, 1)$value
  integrand <- function(alpha) {
    return(mu0(alpha)*px(alpha))
  }
  theta0 <- (integrate(integrand, 0, 1)$value)/p1
  #####################################################
  
  output <- as.matrix(read.table(paste("C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\rawoutput_par", j ,"txt", sep="."), sep = ","))
  
  rep <- 500
  rownum <- rep(1:5, rep)
  est   <- output[rownum==2,] 
  se    <- output[rownum==3,]
  rej   <- output[rownum==4,]
  ci    <- output[rownum==5,]
  k.cv  <- output[rownum==1,4]
  
  colid <- 4   # V4 = K_cv
  
  MAD <- median(abs(est[, colid] - theta0),na.rm = T)
  CR  <- 1 - mean(rej[, colid],na.rm = T)
  MAL <- median(ci[, colid],na.rm = T)
  
  table <- round(rbind(MAD, CR, MAL), 3)
  rownames(table) <- c("Median Abs. Dev.", "95% CI Coverage", "Median CI Length")
  results[j, ] <- c(MAD, CR, MAL)
}

write.csv(results, file = "C:\\Users\\12833\\OneDrive\\LMU\\第三学期\\SEMINAR\\Other Literature\\Feng's Paper\\results.csv", row.names = TRUE)
  