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

# for a given j, compute all the vector c_j by trapezoidal integration

c_j <- function(cs_omega, j, omega, delta_x, steps){
  exp_array <- array(rep(complex(real = cos(2 * pi * j * omega), imaginary = sin(2 * pi * j * omega)), each = p ^ 2), dim = c(p, p, steps))
  int <- cs_omega * exp_array
  return(delta_x * (rowSums(int, dims = 2) - 0.5 * (int[, , 1] + int[, , steps])))
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

# spectral density matrix plot

delta_x <- 0.0005
omega <- seq(0, 0.5, delta_x)
steps <- length(omega)
fxx <- array(dim = c(2, 2, steps))
c_omega_plot <- array(dim = c(2, 2, steps))
for(i in 1:steps){
  alpha <- diag(2) + complex(real = cos(2 * pi * omega[i]), imaginary = -sin(2 * pi * omega[i])) * a1
  beta <- diag(2) + complex(real = cos(2 * pi * omega[i]), imaginary = -sin(2 * pi * omega[i])) * b1
  G <- solve(alpha) %*% beta
  fxx[, , i] <- G %*% Sigma %*% t(Conj(G)) / (2 * pi)
  c_omega_plot[, , i] <- eigen(fxx[, , i])$vectors
  for(j in 1:2){
    c_omega_plot[, j, i] <- c_omega_plot[, j, i] * sign(Re(c_omega_plot[1, j, i]))
  }
}

library(latex2exp)

# plots for spectral matrix and eigen analysis

par(mfcol = c(2, 2), mar = c(4, 4, 1, 1) + 0.1)
plot(x = omega, y = Re(fxx[1, 1, ]), type = "l", xlab = quote(omega), ylab = TeX("$h_{X}( \\omega )_{11}$"), ylim = range(c(0, Re(fxx[1, 1, ]))))
lines(x = omega, y = Im(fxx[1, 1, ]), col = 2)
plot(x = omega, y = Re(fxx[2, 1, ]), ylim = range(c(Re(fxx[2, 1, ]), Im(fxx[2, 1, ]))), type = "l", xlab = quote(omega),
     ylab = TeX("$h_{X}( \\omega )_{21}$"))
lines(x = omega, y = Im(fxx[2, 1, ]), col = 2)
plot(x = omega, y = Re(fxx[1, 2, ]), ylim = range(c(Re(fxx[1, 2, ]), Im(fxx[1, 2, ]))), type = "l", xlab = quote(omega),
     ylab = TeX("$h_{X}( \\omega )_{12}$"))
lines(x = omega, y = Im(fxx[1, 2, ]), col = 2)
plot(x = omega, y = Re(fxx[2, 2, ]), type = "l", xlab = quote(omega), ylab = TeX("$h_{X}( \\omega )_{22}$"), ylim = range(c(0, Re(fxx[2, 2, ]))))
lines(x = omega, y = Im(fxx[2, 2, ]), col = 2)

par(mfcol = c(2, 2), mar = c(4, 4, 1, 1) + 0.1)
plot(x = omega, y = Re(c_omega_plot[1, 1, ]), type = "l", xlab = quote(omega), ylab = TeX("$c(\\omega)_{11}$"), ylim = c(0, 1))
lines(x = omega, y = Im(c_omega_plot[1, 1, ]), col = 2)
plot(x = omega, y = Re(c_omega_plot[2, 1, ]), ylim = range(c(Re(c_omega_plot[2, 1, ]), Im(c_omega_plot[2, 1, ]))), type = "l", xlab = quote(omega),
     ylab = TeX("$c(\\omega)_{21}$"))
lines(x = omega, y = Im(c_omega_plot[2, 1, ]), col = 2)
plot(x = omega, y = Re(c_omega_plot[1, 2, ]), ylim = range(c(Re(c_omega_plot[1, 2, ]), Im(c_omega_plot[1, 2, ]))), type = "l", xlab = quote(omega),
     ylab = TeX("$c(\\omega)_{12}$"))
lines(x = omega, y = Im(c_omega_plot[1, 2, ]), col = 2)
plot(x = omega, y = Re(c_omega_plot[2, 2, ]), type = "l", xlab = quote(omega), ylab = TeX("$c(\\omega)_{22}$"), ylim = c(-1, 1))
lines(x = omega, y = Im(c_omega_plot[2, 2, ]), col = 2)

# simulation
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

par(mfrow = c(2, 1), mar =  c(4, 4, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$x_1$"))
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$x_2$"))

acf(x)

pacf(x)

# function c(omega) evaluated for trapezoidal integration

delta_x <- 0.0005
omega <- seq(-0.5, 0.5, delta_x)
steps <- length(omega)
cs_omega_c <- array(dim = c(p, p, steps))
cs_omega_b <- array(dim = c(p, p, steps))
for(i in 1:steps){
  cs_omega_c[, , i] <- c_omega(a = a, b = b, omega = omega[i], Sigma = Sigma, coefficient = "c", p = p, k = k, l = l)
  cs_omega_b[, , i] <- c_omega(a = a, b = b, omega = omega[i], Sigma = Sigma, coefficient = "b", p = p, k = k, l = l)
}

# integrand analysis fo Gaussian integration

for(j in c(0, 1, 5)){
  exp_array <- array(rep(complex(real = cos(2 * pi * j * omega), imaginary = sin(2 * pi * j * omega)), each = p ^ 2), dim = c(p, p, steps))
  int <- cs_omega_c * exp_array
  pdf(file = paste0("C:\\Users\\jontx\\Google Drive\\TFG Jon\\Memoria\\imagenes\\c", j, ".pdf"), width = 4, height = 4)
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(x = omega, y = Re(int[1, 1, ]), type = "l", xlab = quote(omega), ylab = TeX("$c(\\omega)e^{2\\pi ij\\omega}$"),
         ylim = range(c(Re(int), Im(int))))
    lines(x = omega, y = Im(int[1, 1, ]), lty = "dashed")
    lines(x = omega, y = Re(int[2, 1, ]), col = 2)
    lines(x = omega, y = Im(int[2, 1, ]), col = 2, lty = "dashed")
    lines(x = omega, y = Re(int[1, 2, ]), col = 3)
    lines(x = omega, y = Im(int[1, 2, ]), col = 3, lty = "dashed")
    lines(x = omega, y = Re(int[2, 2, ]), col = 4)
    lines(x = omega, y = Im(int[2, 2, ]), col = 4, lty = "dashed")
  dev.off()
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

j <- (1 - trunc((1 + mirror) * size)):trunc((2 + mirror) * size)
c_j_plot <- array(dim = c(2, 2, length(j)))
for(i in 1:length(j)){
  c_j_plot[, , i] <- c_j_computed[[as.character(j[i])]]
}


# c_j figure

par(mfrow = c(2, 2), mar =  c(4.1, 4, 1, 1))
plot(1, type = "n", xlim = c(-9, 9), ylim = range(Re(c_j_plot[1, 1, ]), Im(c_j_plot[1, 1, ])), xlab = "j", ylab = TeX("$c_{j}$"))
abline(v = 0, h = 0, col = "gray")
points(x = j[j < 10 & j > -10], y = Re(c_j_plot[1, 1, j < 10 & j > -10]), pch = 16)
points(x = j[j < 10 & j > -10], y = Im(c_j_plot[1, 1, j < 10 & j > -10]), col = 2, pch = 16)
plot(1, type = "n", xlim = c(-9, 9), ylim = range(Re(c_j_plot[2, 1, ]), Im(c_j_plot[2, 1, ])), xlab = "j", ylab = TeX("$c_{j}$"))
abline(v = 0, h = 0, col = "gray")
points(x = j[j < 10 & j > -10], y = Re(c_j_plot[2, 1, j < 10 & j > -10]), pch = 16)
points(x = j[j < 10 & j > -10], y = Im(c_j_plot[2, 1, j < 10 & j > -10]), col = 2, pch = 16)
plot(1, type = "n", xlim = c(-9, 9), ylim = range(Re(c_j_plot[1, 2, ]), Im(c_j_plot[1, 2, ])), xlab = "j", ylab = TeX("$c_{j}$"))
abline(v = 0, h = 0, col = "gray")
points(x = j[j < 10 & j > -10], y = Re(c_j_plot[1, 2, j < 10 & j > -10]), pch = 16)
points(x = j[j < 10 & j > -10], y = Im(c_j_plot[1, 2, j < 10 & j > -10]), col = 2, pch = 16)
plot(1, type = "n", xlim = c(-9, 9), ylim = range(Re(c_j_plot[2, 2, ]), Im(c_j_plot[2, 2, ])), xlab = "j", ylab = TeX("$c_{j}$"))
abline(v = 0, h = 0, col = "gray")
points(x = j[j < 10 & j > -10], y = Re(c_j_plot[2, 2, j < 10 & j > -10]), pch = 16)
points(x = j[j < 10 & j > -10], y = Im(c_j_plot[2, 2, j < 10 & j > -10]), col = 2, pch = 16)

# number of principal components

n <- p

# y series computation

y <- matrix(0, nrow = size, ncol = n)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    y[t, ] <- y[t, ] + c_j_computed[[as.character(t - i)]][1:n, ] %*% x[k, ]
  }
}

par(mfrow = c(2, 1), mar =  c(4, 4, 2, 1))
plot(Re(y[, 1]), type = "l", xlab = "Time", ylab = TeX("$y_1$"))
plot(Re(y[, 2]), type = "l", xlab = "Time", ylab = TeX("$y_2$"))

# reconstruction of original series
x_hat <- matrix(0, nrow = size, ncol = p)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    if(n == 1){
      x_hat[t, ] <- x_hat[t, ] + b_j_computed[[as.character(t - i)]][, 1] * y[k, ]
    }else{
      x_hat[t, ] <- x_hat[t, ] + b_j_computed[[as.character(t - i)]][, 1:n] %*% y[k, ]
    }
  }
}

# series reconstruction

par(mfrow = c(2, 1), mar =  c(4, 4.75, 2, 1))
plot(x[, 1], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_1$"), ylim = c(min(x[, 1]), max(x[, 1]) + 0.1))
lines(Re(x_hat[, 1]), col = 2)
# legend("top", legend = c("Original", "Reconstruction"), lty = 1, col = 1:2)
plot(x[, 2], type = "l", xlab = "Time", ylab = TeX("$\\hat{x}_2$"), ylim = c(min(x[, 2]), max(x[, 2]) + 0.1))
lines(Re(x_hat[, 2]), col = 2)
# legend("top", legend = c("Original", "Reconstruction"), lty = 1, col = 1:2)

# 2D representation

layout(1)
plot(x = x[, 1], y = x[, 2], type = "l", xlab = TeX("$x_1$"), ylab = TeX("$x_2$"))
points(x = c(x[1, 1], x[size, 1]), y = c(x[1, 2], x[size, 2]), pch = c(16, 8), col = "blue")
abline(h = 0, v = 0)
lines(x = Re(x_hat[, 1]), y = Re(x_hat[, 2]), col = 2, lty = "dashed")
points(x = c(Re(x_hat[1, 1]), Re(x_hat[size, 1])), y = c(Re(x_hat[1, 2]), Re(x_hat[size, 2])), pch = c(16, 8), col = "green")
legend("topleft", legend = c("x", "x_hat"), col = 1:2, lty = c("solid", "dashed"))

# # colored evolution
# 
# library(viridis)
# 
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