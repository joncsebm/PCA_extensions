# number of different standard deviations (i.e. plots)
p <- 4

# standard deviations distribution
std_dev <- seq(from = 0, to = 4, length.out = p)

# number of points
n <- 200

# correlation coefficient between 'x' & 'y'
rho <- 0.2

# variable simulation

library(mvtnorm)
set.seed(1)
x <- rmvnorm(n = n, sigma = diag(rep(5, times = p)))
y <- rho * x + mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = diag(std_dev ^ 2))

# Points and PCs

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1) + 0.1)
for(i in 1:p){
  pca <- princomp(cbind(x[, i], y[, i]))
  plot(x[, i], y[, i], asp = 1, pch = 16, xlab = 'x', ylab = 'y', cex = 0.5)
  p0 <- c(mean(x[, i]), mean(y[, i]))
  v1 <- pca$sdev[1] ^ 2 * pca$loadings[, 1]
  v2 <- pca$sdev[2] ^ 2 * pca$loadings[, 2]
  arrows(x0 = p0[1], y0 = p0[2], x1 = p0[1] + 0.75 * v1[1],
         y1 = p0[2] + 0.75 * v1[2], col = 2, lwd = 3, length = 0.1)
  arrows(x0 = p0[1], y0 = p0[2], x1 = p0[1] + 0.75 * v2[1],
         y1 = p0[2] + 0.75 * v2[2], col = 3, lwd = 3, length = 0.1)
  legend("topleft", legend = c("PC1", "PC2"), lwd = 2, col = 2:3)
}


# x-axis projections

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1) + 0.1)
for(i in 1:p){
  plot(x[, i], rep(0, times = length(x[, i])), asp = 1, pch = 16, xlab = 'x', ylab = 'y', cex = 0.5)
}

# y-axis projections

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1) + 0.1)
for(i in 1:p){
  plot(rep(0, times = length(y[, i])), y[, i], asp = 1, pch = 16, xlab = 'x', ylab = 'y', ylim = range(y), cex = 0.5)
}

# PC's projections

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1) + 0.1)
for(i in 1:p){
  pca <- princomp(cbind(x[, i], y[, i]))
  PC1 <- pca$loadings
  p0 <- c(mean(x[, i]), mean(y[, i]))
  plot(x = p0[1] + (PC1[1] * (PC1 %*% rbind(x[, 1], y[, 1]))),
       y = p0[2] + (PC1[2] * (PC1 %*% rbind(x[, 1], y[, 1]))),
       asp = 1, pch = 16, xlab = 'x', ylab = 'y', xlim = range(x) / 1.5, ylim = range(y) / 1.5, cex = 0.5)
}

# Final plot with all projections

par(mfrow = c(3, 4), mar = c(2, 2, 2, 2))
for(i in 1:4){
  plot(x[, i], rep(0, 200), asp = 1, pch = 16, xlab = "x", ylab = "y", xlim = range(x) / 1.5, ylim = range(y) / 1.5, cex = 0.5)
}
  
for(i in 1:4){
  plot(rep(0, 200), y[, i], asp = 1, pch = 16, xlab = "x", ylab = "y", xlim = range(x) / 1.5, ylim = range(y) / 1.5, cex = 0.5)
}
  
for(i in 1:4){
  pca <- princomp(cbind(x[, i], y[, i]))
  PC1 <- pca$loadings
  p0 <- c(mean(x[, i]), mean(y[, i]))
  plot(x = p0[1] + (PC1[1] * (PC1 %*% rbind(x[, 1], y[, 1]))),
       y = p0[2] + (PC1[2] * (PC1 %*% rbind(x[, 1], y[, 1]))),
       asp = 1, pch = 16, xlab = 'x', ylab = 'y', xlim = range(x) / 1.5, ylim = range(y) / 1.5, cex = 0.5)
}