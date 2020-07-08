# eigenvectors of spectral matrix for a given omega

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

# Gaussian integration
  
c_j <- function(cs_omega, weights, j, omega, steps){
  exp_array <- array(rep(complex(real = cos(2 * pi * j * omega), imaginary = sin(2 * pi * j * omega)), each = p ^ 2), dim = c(p, p, steps))
  int <- array(dim = c(p, p, steps))
  for(i in 1:steps){
    int[, , i] <- weights[i] * cs_omega[, , i] * exp_array[, , i]
  }
  return(rowSums(int, dims = 2))
}

# number of variables in the model
p <- 2
# AR order
k <- 1
# MA order
l <- 0

# model parameters

Sigma <- 0.1 * diag(p)
a1 <- rbind(c(0.5, 0.2), c(0.1, 0.75))
b1 <- rbind(c(0.7, 0.2), c(0.1, 0.75))

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

size <- 100
library(mvtnorm)
set.seed(2)
epsilon <- rmvnorm(n = size, sigma = Sigma)
x <- epsilon
for(i in (max(k, l) + 1):size){
  if(k > 0){
    for(j in 1:k){
      x[i, ] <- x[i, ] + a[, , j] %*% x[i - j, ]
    }
  }
  if(l > 0){
    for(j in 1:l){
      x[i, ] <- x[i, ] + b[, , j] %*% epsilon[i - j, ]
    }
  }
}

# function c(omega) evaluated for Gaussian integration

steps <- 1280
omega <- sphunif::Gauss_Legen_nodes(a = -0.5, b = 0.5, N = steps)
weights <- sphunif::Gauss_Legen_weights(a = -0.5, b = 0.5, N = steps)
cs_omega_c <- array(dim = c(p, p, steps))
cs_omega_b <- array(dim = c(p, p, steps))
for(i in 1:steps){
  cs_omega_c[, , i] <- c_omega(a = a, b = b, omega = omega[i], Sigma = Sigma, coefficient = "c", p = p, k = k, l = l)
  cs_omega_b[, , i] <- c_omega(a = a, b = b, omega = omega[i], Sigma = Sigma, coefficient = "b", p = p, k = k, l = l)
}

# c_j and b_j coefficient for series reconstruction

c_j_computed <- list()
for(i in (1 - size):(2 * size)){
  c_j_computed[[as.character(i)]] <- Conj(t(c_j(cs_omega = cs_omega_c, weights = weights, j = i, omega = omega, steps = steps)))
}
b_j_computed <- list()
for(i in (1 - size):(2 * size)){
  b_j_computed[[as.character(i)]] <- c_j(cs_omega = cs_omega_b, weights = weights, j = i, omega = omega, steps = steps)
}

# y series computation

y <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in 1:size){
    y[t, ] <- y[t, ] + c_j_computed[[as.character(t - i)]] %*% x[i, ]
  }
}

# reconstruction of original series

x_hat <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in 1:size){
    x_hat[t, ] <- x_hat[t, ] + b_j_computed[[as.character(t - i)]] %*% y[i, ]
  }
}