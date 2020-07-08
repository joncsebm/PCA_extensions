Cov_matrix <- read.table(file = "matrix.txt", col.names = paste0("X", 1:10), row.names = paste0("X", 1:10))
eigen <- eigen(Cov_matrix, symmetric = TRUE)

spca <- list()

library(elasticnet)
lambdas <- c(0, 0.1, 1, 10, 100, 1000)
for(i in lambdas){
  spca[[paste0("lambda", i)]] <- spca(x = Cov_matrix, K = nrow(Cov_matrix), type = "Gram",
                                       sparse = "penalty", para = rep(x = i, times = nrow(Cov_matrix)))
  spca[[paste0("lambda", i)]]$cumpev <- spca[[paste0("lambda", i)]]$pev
  for(j in 2:10){
    spca[[paste0("lambda", i)]]$cumpev[j] <- spca[[paste0("lambda", i)]]$cumpev[j] + spca[[paste0("lambda", i)]]$cumpev[j - 1]
  }
}

library(latex2exp)

plot(100 * spca$lambda0$pev, col = 1, type = "l", ylab = "Explained variance (%)", xaxt = "n")
for(i in 2:10){
	lines(100 * spca[[paste0("lambda", lambdas[i])]]$pev, col = i)
}
axis(side = 1, at = 1:10, labels = paste0('PC', 1:10))
legend("topright", title = TeX("Lasso $\\lambda$"), legend = TeX(paste0("$\\lambda = ", lambdas, "$")), col = 1:6, lty = 1)

plot(100 * spca$lambda0$cumpev, col = 1, type = "l", ylab = "Cumulative explained variance (%)", xaxt = "n",
     ylim = 100 * c(0, 1))

for(i in 2:10){
    lines(100 * spca[[paste0("lambda", lambdas[i])]]$cumpev, col = i)
}
  axis(side = 1, at = 1:10, labels = paste0('PC', 1:10))
  legend("bottomright", title = TeX("Lasso $\\lambda$"), legend = lambdas, col = 1:6, lty = 1)

# library('xtable')
# xtable(x = cbind(eigen$vectors[, 1:2], spca$lambda0.1$loadings[, 1:2], spca$lambda1$loadings[, 1:2], spca$lambda10$loadings[, 1:2],
#        spca$lambda100$loadings[, 1:2], spca$lambda1000$loadings[, 1:2]))
# 
# xtable(x = cbind(eigen$vectors[, 3:4], spca$lambda0.1$loadings[, 3:4], spca$lambda1$loadings[, 3:4], spca$lambda10$loadings[, 3:4],
#                  spca$lambda100$loadings[, 3:4], spca$lambda1000$loadings[, 3:4]))


# # simulation
# n <- 1e7
# V1 <- rnorm(n = n, mean = 0, sd = sqrt(290))
# V2 <- rnorm(n = n, mean = 0, sd = sqrt(300))
# V3 <- -0.3 * V1 + 0.925 * V2 + rnorm(n = n)
# matriz_V <- cbind(V1, V1, V1, V1, V2, V2, V2, V2, V3, V3)
# matriz_epsilon <- mvtnorm::rmvnorm(n = n, sigma = diag(x = 1, nrow = 10, ncol = 10))
# X <- matriz_V + matriz_epsilon
# cov <- t(X) %*% X / n
# rownames(cov) <- paste('X', 1:10, sep = '')
# colnames(cov) <- rownames(cov)