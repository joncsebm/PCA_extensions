# univariate autorregresive AR(3)

a1 <- 0.3
a2 <- -0.2
a3 <- -0.1

n <- 1000
set.seed(1)
epsilon <- rnorm(n = n, sd = 0.2)
x <- vector(mode = "numeric", length = n)
x[1:3] <- epsilon[1:3]
for(t in 4:n){
  x[t] <- epsilon[t] - a1 * x[t - 1] - a2 * x[t - 2] - a3 * x[t - 3]
}

pacf(x)

plot.ts(x)


# univariate moving average MA(2)

b1 <- 0.7
b2 <- -0.2

n <- 1000
epsilon <- rnorm(n = n, sd = 0.2)
x <- vector(mode = "numeric", length = n)
x[1:3] <- epsilon[1:3]
for(t in 4:n){
  x[t] <- epsilon[t] + b1 * epsilon[t - 1] + b2 * epsilon[t - 2]
}

acf(x)

plot.ts(x)

# univariate autorregresive ARMA(1,1)

a1 <- 0.5
b1 <- 0.2

n <- 1000
epsilon <- rnorm(n = n, sd = 0.2)
x <- vector(mode = "numeric", length = n)
x[1] <- epsilon[1]
for(t in 2:n){
  x[t] <- epsilon[t] - a1 * x[t - 1] + b1 * epsilon[t - 1]
}

acf(x)

pacf(x)

plot.ts(x)

# multivariate moving average VMA(1)

b1 <- rbind(c(0.3, 0.1), c(-0.2, 0.7))
n <- 1000
library(mvtnorm)
epsilon <- rmvnorm(n = n, sigma = diag(0.2, 2))
x <- matrix(nrow = n, ncol = 2)
x[1, ] <- epsilon[1, ]
for(t in 2:n){
  x[t, ] <- epsilon[t, ] + b1 %*% epsilon[t - 1, ]
}

plot.ts(x)

# multivariate autoregressive VAR(1)

a1 <- rbind(c(0.2, 0.1), c(-0.1, 0.6))
n <- 1000
library(mvtnorm)
epsilon <- rmvnorm(n = n, sigma = diag(0.2, 2))
x <- matrix(nrow = n, ncol = 2)
x[1, ] <- epsilon[1, ]
for(t in 2:n){
  x[t, ] <- epsilon[t, ] - a1 %*% x[t - 1, ]
}

plot.ts(x)

# multivariate ARMA VARMA(1, 1)

a1 <- rbind(c(0.2, 0.1), c(-0.1, 0.6))
b1 <- rbind(c(0.3, 0.1), c(-0.2, 0.7))
n <- 1000
library(mvtnorm)
epsilon <- rmvnorm(n = n, sigma = diag(0.2, 2))
x <- matrix(nrow = n, ncol = 2)
x[1, ] <- epsilon[1, ]
for(t in 2:n){
  x[t, ] <- epsilon[t, ] - a1 %*% x[t - 1, ] + b1 %*% epsilon[t - 1, ]
}

plot.ts(x)