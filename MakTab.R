#########################################################
##### This file produces summary table for Deaner  ######
#########################################################
library(Hmisc)


# initialize the result matrix
results <- matrix(NA, nrow=6, ncol=3)
colnames(results) <- c("Median Abs. Dev.", "95% CI Coverage", "Median CI Length")
rownames(results) <- paste0("Setting", 1:6)

for (j in 1:6) {
  
  # load the result file
  output <- read.csv(paste0("output_deaner_par", j, ".txt"))
  
  # key columns
  MAD <- median(abs(output$bias), na.rm = TRUE)
  CR  <- 1 - mean(output$rej, na.rm = TRUE)
  MAL <- median(output$ci_length, na.rm = TRUE)
  
  results[j, ] <- c(MAD, CR, MAL)
}

# save the result
write.csv(results, 
          file = "results_deaner.csv", 
          row.names = TRUE)
