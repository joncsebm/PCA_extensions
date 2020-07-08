# Function to compute the cross product of two vectors

cross <- function(v1, v2){
  return(c(v1[2] * v2[3] - v1[3] * v2[2], v1[3] * v2[1] - v1[1] * v2[3], v1[1] * v2[2] - v1[2] * v2[1]))
}

# Eigenvectors of spectral matrix for a given omega

c_omega <- function(a, b, omega, Sigma, coefficient, p, k, l){
  # MA transfer function
  alpha <- diag(p)
  beta <- diag(p)
  if(k > 0){
    for(i in 1:k){
      alpha <- alpha + a[, , i] * complex(real = cos(2 * pi * omega * i), imaginary = -sin(2 * pi * omega * i))
    }
  }
  if(l > 0){
    for(i in 1:l){
      beta <- beta + b[, , i] * complex(real = cos(2 * pi * omega * i), imaginary = -sin(2 * pi * omega * i))
    }
  }
  # ARMA transfer function
  G <- solve(alpha) %*% beta
  # MA spectral matrix
  h <- G %*% Sigma %*% t(Conj(G)) / (2 * pi)
  # eigen analysis
  vectors <- eigen(h)$vectors
  # returns eigenvectors in case of c(omega) and their conjugates in case of b(omega)
  for(i in 1:p){
    vectors[, i] <- vectors[, i] * sign(Re(vectors[1, i]))
  }
  if(coefficient == "c"){
    return(vectors)
  }else{
    return(Conj(vectors))
  }
}

# For a given j, compute all the vector c_j by trapezoidal integration

c_j <- function(cs_omega, j, omega, delta_x, steps){
  exp_array <- array(rep(complex(real = cos(2 * pi * j * omega), imaginary = sin(2 * pi * j * omega)), each = p ^ 2), dim = c(p, p, steps))
  int <- cs_omega * exp_array
  return(delta_x * (rowSums(int, dims = 2) - 0.5 * (int[, , 1] + int[, , steps])))
}

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# orthonormal base (1)

v1 <- c(1, 2, 3) / sqrt(14)
v2 <- c(3, 2, -7 / 3) * 3 / sqrt(166)
v3 <- cross(v1 = v1, v2 = v2)

# coefficient matrix computation

a1 <- matrix(0, 3, 3)
lambda <- c(0.9, 0.01, 0.001)

for(i in 1:3){
  v <- get(paste0("v", i))
  a1 <- a1 + lambda[i] * (v %*% t(v))
}

# orthonormal base (2)

v1 <- c(1, 1, 1) / sqrt(3)
v2 <- c(1, 1, -2) / sqrt(6)
v3 <- cross(v1 = v1, v2 = v2)

# variances matrix computation

Sigma <- matrix(0, 3, 3)
lambda <- c(100, 1, 0.01)

for(i in 1:3){
  v <- get(paste0("v", i))
  Sigma <- Sigma + lambda[i] * (v %*% t(v))
}

# number of variables in the model
p <- 3
# AR order
k <- 1
# MA order
l <- 0

# information from model extraction

a <- array(dim = c(p, p, k))
if(k > 0){
  for(i in 1:k){
    a[, , i] <- get(paste0("a", i))
  }
}
b <- array(dim = c(p, p, l))
if(l > 0){
  for(i in 1:l){
    b[, , i] <- get(paste0("b", i))
  }
}

# model simulation

size <- 500
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

# Function c(omega) evaluated for trapezoidal integration

delta_x <- 0.0005
omega <- seq(-0.5, 0.5, delta_x)
steps <- length(omega)
cs_omega_c <- array(dim = c(p, p, steps))
cs_omega_b <- array(dim = c(p, p, steps))
for(i in 1:steps){
  cs_omega_c[, , i] <- c_omega(a = a, b = b, omega = omega[i], Sigma = Sigma, coefficient = "c", p = p, k = k, l = l)
  cs_omega_b[, , i] <- c_omega(a = a, b = b, omega = omega[i], Sigma = Sigma, coefficient = "b", p = p, k = k, l = l)
}

mirror <- 0

# c_j and b_j coefficient for series reconstruction

c_j_computed <- list()
for(i in (1 - trunc((1 + mirror) * size)):trunc((2 + mirror) * size)){
  c_j_computed[[as.character(i)]] <- Conj(t(c_j(cs_omega = cs_omega_c, j = i, omega = omega, delta_x = delta_x, steps = steps)))
}

b_j_computed <- list()
for(i in (1 - trunc((1 + mirror) * size)):trunc((2 + mirror) * size)){
  b_j_computed[[as.character(i)]] <- c_j(cs_omega = cs_omega_b, j = i, omega = omega, delta_x = delta_x, steps = steps)
}

library(latex2exp)

par(mfrow = c(3, 1), mar =  c(4, 4, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$x_1$"))
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$x_2$"))
plot(x[, 3], type = "l", xlab = "Time", ylab = TeX("$x_3$"))


# number of principal components
n <- 1

# y series computation

y_1 <- matrix(0, nrow = size, ncol = n)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    y_1[t, ] <- y_1[t, ] + c_j_computed[[as.character(t - i)]][1:n, ] %*% x[k, ]
  }
}

# reconstruction of original series

x_hat_1 <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    if(n == 1){
      x_hat_1[t, ] <- x_hat_1[t, ] + b_j_computed[[as.character(t - i)]][, 1] * y_1[k, ]
    }else{
      x_hat_1[t, ] <- x_hat_1[t, ] + b_j_computed[[as.character(t - i)]][, 1:n] %*% y_1[k, ]
    }
  }
}

# series reconstruction 1 PC

par(mfrow = c(3, 1), mar =  c(4, 4.75, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_1$"), ylim = c(min(x[, 1]), max(x[, 1]) + 0.1))
lines(Re(x_hat[, 1]), col = 2)
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_2$"), ylim = c(min(x[, 2]), max(x[, 2]) + 0.1))
lines(Re(x_hat[, 2]), col = 2)
plot(x[, 3], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_3$"), ylim = c(min(x[, 3]), max(x[, 3]) + 0.1))
lines(Re(x_hat[, 3]), col = 2)


# number of principal components
n <- 2

# y series computation

y_2 <- matrix(0, nrow = size, ncol = n)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    y_2[t, ] <- y_2[t, ] + c_j_computed[[as.character(t - i)]][1:n, ] %*% x[k, ]
  }
}

# reconstruction of original series

x_hat_2 <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    if(n == 1){
      x_hat_2[t, ] <- x_hat_2[t, ] + b_j_computed[[as.character(t - i)]][, 1] * y_2[k, ]
    }else{
      x_hat_2[t, ] <- x_hat_2[t, ] + b_j_computed[[as.character(t - i)]][, 1:n] %*% y_2[k, ]
    }
  }
}

# series reconstruction 2 PC

par(mfrow = c(3, 1), mar =  c(4, 4.75, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_1$"), ylim = c(min(x[, 1]), max(x[, 1]) + 0.1))
lines(Re(x_hat[, 1]), col = 2)
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_2$"), ylim = c(min(x[, 2]), max(x[, 2]) + 0.1))
lines(Re(x_hat[, 2]), col = 2)
plot(x[, 3], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_3$"), ylim = c(min(x[, 3]), max(x[, 3]) + 0.1))
lines(Re(x_hat[, 3]), col = 2)

# number of principal components

n <- 3

# y series computation

y_3 <- matrix(0, nrow = size, ncol = n)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    y_3[t, ] <- y_3[t, ] + c_j_computed[[as.character(t - i)]][1:n, ] %*% x[k, ]
  }
}

# reconstruction of original series

x_hat_3 <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    if(n == 1){
      x_hat_3[t, ] <- x_hat_3[t, ] + b_j_computed[[as.character(t - i)]][, 1] * y_3[k, ]
    }else{
      x_hat_3[t, ] <- x_hat_3[t, ] + b_j_computed[[as.character(t - i)]][, 1:n] %*% y_3[k, ]
    }
  }
}

# series reconstruction 3 PC

par(mfrow = c(3, 1), mar =  c(4, 4.75, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_1$"), ylim = c(min(x[, 1]), max(x[, 1]) + 0.1))
lines(Re(x_hat[, 1]), col = 2)
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_2$"), ylim = c(min(x[, 2]), max(x[, 2]) + 0.1))
lines(Re(x_hat[, 2]), col = 2)
plot(x[, 3], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_3$"), ylim = c(min(x[, 3]), max(x[, 3]) + 0.1))
lines(Re(x_hat[, 3]), col = 2)

ylim1 <- c(-20, 30)
ylim2 <- c(-40, 40)
ylim3 <- c(-60, 60)

# Figure for original series
par(mfrow = c(3, 1), mar = c(2, 4, 2, 2) + 0.1)
for(i in 1:3){
  plot(x[, i], type = "l", xlab = "", ylab = TeX(paste0("x_", i)), ylim = get(paste0("ylim", i)))
}

for(j in 1:3){
  par(mfrow = c(3, 1), mar = c(2, 4.5, 2, 2) + 0.1)
  for(i in 1:3){
    plot(Re(get(paste0("x_hat_", j))[, i]), type = "l", xlab = "", ylab = TeX(paste0("$\\hat{x}_", i, "$")), ylim = get(paste0("ylim", i)))
  }
}