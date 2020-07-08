# for a given j, compute all the vector c_j by trapezoidal integration

c_j <- function(h, j, p){
  omega <- h$freq
  delta_x <- omega[2] - omega[1]
  steps <- length(omega)
  c_omega <- array(dim = c(p, p, steps))
  for(i in 1:steps){
    c_omega[, , i] <- eigen(h$fxx[, , i])$vectors
  }
  exp_array_1 <- array(rep(complex(real = cos(2 * pi * j * omega), imaginary = sin(2 * pi * j * omega)), each = p ^ 2), dim = c(p, p, steps))
  int_1 <- c_omega * exp_array_1
  exp_array_2 <- array(rep(complex(real = cos(2 * pi * j * omega), imaginary = -sin(2 * pi * j * omega)), each = p ^ 2), dim = c(p, p, steps))
  int_2 <- c_omega * exp_array_2
  integral_1 <- delta_x * (rowSums(int_1, dims = 2) - 0.5 * (int_1[, , 1] + int_1[, , steps]))
  integral_2 <- delta_x * (rowSums(int_2, dims = 2) - 0.5 * (int_2[, , 1] + int_2[, , steps]))
  return(integral_1 + integral_2)
}

# number of variables in the model
p <- 2
# AR order
k <- 1
# MA order
l <- 1

# model parameters

Sigma <- 0.01 * diag(p)
a1 <- rbind(c(0.3, -0.1), c(0.1, 0.9))
b1 <- rbind(c(0.8, 0.2), c(-0.2, 0.5))

a <- array(dim = c(p, p, k))
if(k > 0){
  for(i in 1:k){
    a[, , i] <- get(paste0('a', i))
  }
}
b <- array(dim = c(p, p, l))
if(l > 0){
  for(i in 1:l){
    b[, , i] <- get(paste0('b', i))
  }
}

# model simulation

size <- 200
library(mvtnorm)
set.seed(1)
epsilon <- rmvnorm(n = size, sigma = Sigma)
x <- epsilon
for(i in (max(k, l) + 1):size){
  if(k > 0){
    for(j in 1:k){
      x[i, ] <- x[i, ] - a[, , j] %*% x[i - j, ]
    }
  }
  if(l > 0){
    for(j in 1:l){
      x[i, ] <- x[i, ] + b[, , j] %*% epsilon[i - j, ]
    }
  }
}

library(astsa)
h_omega <- mvspec(x, plot = FALSE, spans = 40)

omega <- seq(from = 0, to = 0.5, length.out = length(h_omega$freq))
steps <- length(omega)
delta_x <- omega[2] - omega[1]
h <- array(dim = c(2, 2, steps))
for(i in 1:steps){
  alpha <- diag(2) + complex(real = cos(2 * pi * omega[i]), imaginary = -sin(2 * pi * omega[i])) * a1
  beta <- diag(2) + complex(real = cos(2 * pi * omega[i]), imaginary = -sin(2 * pi * omega[i])) * b1
  G <- solve(alpha) %*% beta
  h[, , i] <- G %*% Sigma %*% t(Conj(G)) / (2 * pi)
}

library(latex2exp)

par(mfcol = c(2, 2), mar = rep(4, 4) + c(0, 0.75, 0, 0))
plot(x = h_omega$freq, y = Re(h_omega$fxx[1, 1, ]), ylim = range(c(Re(h_omega$fxx[1, 1, ]), Im(h_omega$fxx[1, 1, ]))), type = "l", xlab = quote(omega),
     ylab = TeX("$\\hat{h}_{X}(\\omega)_{11}$"))
lines(x = h_omega$freq, y = Im(h_omega$fxx[1, 1, ]), col = 2)
plot(x = h_omega$freq, y = Re(h_omega$fxx[2, 1, ]), type = "l", xlab = quote(omega),
     ylab = TeX("$\\hat{h}_{X}(\\omega)_{21}$"), ylim = range(c(Re(h_omega$fxx[2, 1, ]), Im(h_omega$fxx[2, 1, ]))))
lines(x = h_omega$freq, y = Im(h_omega$fxx[2, 1, ]), col = 2)
plot(x = h_omega$freq, y = Re(h_omega$fxx[1, 2, ]), type = "l", xlab = quote(omega),
     ylab = TeX("$\\hat{h}_{X}(\\omega)_{12}$"), ylim = range(c(Re(h_omega$fxx[1, 2, ]), Im(h_omega$fxx[1, 2, ]))))
lines(x = h_omega$freq, y = Im(h_omega$fxx[1, 2, ]), col = 2)
plot(x = h_omega$freq, y = Re(h_omega$fxx[2, 2, ]), ylim = range(c(Re(h_omega$fxx[2, 2, ]), Im(h_omega$fxx[2, 2, ]))), type = "l", xlab = quote(omega),
     ylab = TeX("$\\hat{h}_{X}(\\omega)_{22}$"))
lines(x = h_omega$freq, y = Im(h_omega$fxx[2, 2, ]), col = 2)


mirror <- 0

# c_j and b_j coefficient for series reconstruction

c_j_computed <- list()
for(i in (1 - trunc((1 + mirror) * size)):trunc((2 + mirror) * size)){
  c_j_computed[[as.character(i)]] <- Conj(t(c_j(h = h_omega, j = i, p = p)))
}

b_j_computed <- list()
for(i in (1 - trunc((1 + mirror) * size)):trunc((2 + mirror) * size)){
  b_j_computed[[as.character(i)]] <- c_j(h = h_omega, j = i, p = p)
}

# y series computation

y <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    y[t, ] <- y[t, ] + c_j_computed[[as.character(t - i)]] %*% x[k, ]
  }
}

# reconstruction of original series

x_hat <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    x_hat[t, ] <- x_hat[t, ] + b_j_computed[[as.character(t - i)]] %*% y[k, ]
  }
}

# checking results

par(mfrow = c(2, 1), mar =  c(4, 4.75, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_1$"), ylim = c(min(x[, 1]), max(x[, 1]) + 0.1))
lines(Re(x_hat[, 1]), col = 2)
# legend("top", legend = c("Original", "Reconstruction"), lty = 1, col = 1:2)
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_2$"), ylim = c(min(x[, 2]), max(x[, 2]) + 0.1))
lines(Re(x_hat[, 2]), col = 2)
# legend("top", legend = c("Original", "Reconstruction"), lty = 1, col = 1:2)


layout(1)
plot(x = x[, 1], y = x[, 2], type = "l", xlab = TeX("$x_1$"), ylab = TeX("$x_2$"))
points(x = c(x[1, 1], x[size, 1]), y = c(x[1, 2], x[size, 2]), pch = c(16, 8), col = "blue")
abline(h = 0, v = 0)
lines(x = Re(x_hat[, 1]), y = Re(x_hat[, 2]), col = 2, lty = "dashed")
points(x = c(Re(x_hat[1, 1]), Re(x_hat[size, 1])), y = c(Re(x_hat[1, 2]), Re(x_hat[size, 2])), pch = c(16, 8), col = "green")
legend("topleft", legend = c("x", "x_hat"), col = 1:2, lty = c("solid", "dashed"))


# library(viridis)

# layout(cbind(1, 2))
# plot(1, type = "n", xlab = "x1", ylab = "x2", xlim = range(c(x[, 1])), ylim = range(c(x[, 2])), asp = 1)
# col <- viridis::inferno(n = size, direction = -1, alpha = 0.75)
# for(i in 1:(length(x[, 1]) - 1)){
#   segments(x0 = x[i, 1], x1 = x[i + 1, 1], y0 = x[i, 2], y1 = x[i + 1, 2], col = col[i], cex = 0.5)
# }
#
# plot(1, type = "n", xlab = "x_hat_1", ylab = "x_hat_2", xlim = range(Re(x_hat[, 1])), ylim = range(Re(x_hat[, 2])), asp = 1)
# col <- viridis::inferno(n = size, direction = -1, alpha = 0.75)
# for(i in 1:(length(x[, 1]) - 1)){
#   segments(x0 = x[i, 1], x1 = x[i + 1, 1], y0 = x[i, 2], y1 = x[i + 1, 2], col = col[i], cex = 0.5)
# }