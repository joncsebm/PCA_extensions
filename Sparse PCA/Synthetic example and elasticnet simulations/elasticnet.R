# Initial parameters
sigma_x <- 2
sigma_y <- 1
rho <- 0.25
Sigma <- matrix(data = c(sigma_x ^ 2, sigma_x * sigma_y * rho, sigma_x * sigma_y * rho, sigma_y ^ 2), nrow = 2, ncol = 2)

  eigen <- eigen(x = Sigma, symmetric = TRUE)
  pca_ord <- list(variances = eigen$values, loadings = eigen$vectors)
  angle_PC1 <- (180 / pi) * atan(pca_ord$loadings[2, 1] / pca_ord$loadings[1, 1])
  angle_PC2 <- (180 / pi) * atan(pca_ord$loadings[2, 2] / pca_ord$loadings[1, 2])
  
# SPCA varying lasso penalty (population covariance matrix)
  library(elasticnet)
  results_lasso <- matrix(ncol = 3, nrow = 2 * length(1:300))
  for(i in 1:300){
    pca_spca_lasso <- spca(x = Sigma, K = 2, para = rep(x = 1e-4 * i, times = 2), type = "Gram", lambda = 1e-50)
    results_lasso[2 * i - 1, 1] <- 1e-4 * i
    results_lasso[2 * i, 1] <- results_lasso[2 * i - 1, 1]
    results_lasso[2 * i - 1, 2] <- pca_spca_lasso$loadings[1, 1]
    results_lasso[2 * i - 1, 3] <- pca_spca_lasso$loadings[1, 2]
    results_lasso[2 * i, 2] <- pca_spca_lasso$loadings[2, 1]
    results_lasso[2 * i, 3] <- pca_spca_lasso$loadings[2, 2]
  }
  results_lasso <- as.data.frame(results_lasso)
  names(results_lasso) <- c("penalisation", "PC1", "PC2")
  angle_PC1_lasso <- vector(mode = "numeric", length = length(1:300))
  angle_PC2_lasso <- vector(mode = "numeric", length = length(1:300))
  for(i in 1:300){
    angle_PC1_lasso[i] <- atan(results_lasso$PC1[2 * i] / results_lasso$PC1[2 * i - 1]) * (180 / pi)
    angle_PC2_lasso[i] <- atan(results_lasso$PC2[2 * i] / results_lasso$PC2[2 * i - 1]) * (180 / pi)
  }
  angle_PC1_lasso <- round(angle_PC1_lasso, digits = 5)
  angle_PC2_lasso <- round(angle_PC2_lasso, digits = 5)
  
# SPCA varying ridge penalty
  library(elasticnet)
  results_ridge <- matrix(ncol = 3, nrow = 2 * length(1:300))
  for(i in 1:300){
    pca_spca_ridge <- spca(x = Sigma, K = 2, para = rep(x = 1e-50, times = 2), lambda = 1e-4 * i, type = "Gram", eps.conv = 1e-10)
    results_ridge[2 * i - 1, 1] <- 1e-4 * i
    results_ridge[2 * i, 1] <- results_ridge[2 * i - 1, 1]
    results_ridge[2 * i - 1, 2] <- pca_spca_ridge$loadings[1, 1]
    results_ridge[2 * i - 1, 3] <- pca_spca_ridge$loadings[1, 2]
    results_ridge[2 * i, 2] <- pca_spca_ridge$loadings[2, 1]
    results_ridge[2 * i, 3] <- pca_spca_ridge$loadings[2, 2]
  }
  results_ridge <- as.data.frame(results_ridge)
  names(results_ridge) <- c("penalisation", "PC1", "PC2")
  angle_PC1_ridge <- vector(mode = "numeric", length = length(1:300))
  angle_PC2_ridge <- vector(mode = "numeric", length = length(1:300))
  for(i in 1:300){
    angle_PC1_ridge[i] <- atan(results_ridge$PC1[2 * i] / results_ridge$PC1[2 * i - 1]) * (180 / pi)
    angle_PC2_ridge[i] <- atan(results_ridge$PC2[2 * i] / results_ridge$PC2[2 * i - 1]) * (180 / pi)
  }
  angle_PC1_ridge <- round(angle_PC1_ridge, digits = 10)
  angle_PC2_ridge <- round(angle_PC2_ridge, digits = 10)
  
  data_frame <- data.frame(cbind(1e-4 * (1:300), angle_PC1_lasso, angle_PC1_ridge, angle_PC2_lasso, angle_PC2_ridge))
  names(data_frame) <- c("Penalisation", "PC1_lasso", "PC1_ridge", "PC2_lasso", "PC2_ridge")

  plot(x = data_frame$Penalisation, y = abs(angle_PC1 - data_frame$PC1_lasso), type = "l", xlab = expression(lambda),
       ylab = "Absolute angular error (º)", ylim = c(0, 10), main = "Angular error with lasso penalty", col = "red")
  
  lines(x = data_frame$Penalisation, y = abs(angle_PC2 - ifelse(data_frame$PC2_lasso != 90, data_frame$PC2_lasso, -data_frame$PC2_lasso)), col = "blue")
  legend("bottomright", legend = c("PC1", "PC2"), col = c("red", "blue"), lty = 1)
  
  plot(x = data_frame$Penalisation, y = abs(angle_PC1 - data_frame$PC1_ridge), type = "l", xlab = expression(lambda),
       ylab = "Absolute angular error (º)", ylim = c(-0.1, 0.1), main = "Angular error with ridge penalty", col = "red")
  lines(x = data_frame$Penalisation, y = abs(angle_PC2 - data_frame$PC2_ridge), col = "blue")
  legend("bottomright", legend = c("PC1", "PC2"), col = c("red", "blue"), lty = 1)
