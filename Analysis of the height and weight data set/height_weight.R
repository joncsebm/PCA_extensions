# Load data
height_weight <- read.csv("height_weight.csv", sep = ",")
names(height_weight) <- c("Height", "Weight")
n <- length(height_weight$Height)

# Height histogram and normal fitting

xh <- min(height_weight$Height):max(height_weight$Height)
yh <- dnorm(xh, mean = mean(height_weight$Height), sd = sd(height_weight$Height) * (n - 1) / n)
hist(height_weight$Height, freq = FALSE, main = "Height density", xlab = "Height")
lines(xh, yh, type = 'l', col = 'red')

# Weight histogram and normal fitting

xw <- min(height_weight$Weight):max(height_weight$Weight)
yw <- dnorm(xw, mean = mean(height_weight$Weight), sd = sd(height_weight$Weight) * (n - 1) / n)
hist(height_weight$Weight, freq = FALSE, main = "Weight density", xlab = "Weight")
lines(xw, yw, type = 'l', col = 'red')


# Multivariate normal fitting
library(mvtnorm)
grid <- 500
X <- cbind(seq(from = min(height_weight$Height), to = max(height_weight$Height), length.out = grid),
           seq(from = min(height_weight$Weight), to = max(height_weight$Weight), length.out = grid))
mu <- c(mean(height_weight$Height), mean(height_weight$Weight))
# covariance matrix is computed with 1 / n instead of 1 / (n - 1)
Sigma <- (n - 1) / n * cov(cbind(height_weight$Height, height_weight$Weight))
Y <- dmvnorm(expand.grid(X[, 1], X[, 2]), mean = mu, sigma = Sigma)

#3-dimensional representations

contour(x = X[, 1], y = X[, 2], z = matrix(Y, nrow = grid, ncol = grid), nlevels = 10, col = viridis::viridis(n = 14),
        xlab = "Height", ylab = "Weight")
