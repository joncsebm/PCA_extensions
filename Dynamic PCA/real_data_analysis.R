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

library(astsa)
library(latex2exp)
x <- fmri1[, -1]
size <- nrow(x)
p <- ncol(x)

par(mfrow = c(2, 4), mar = c(4, 5, 2, 1))
for(i in 1:8){
  plot(x[, i], xlab = "Time", ylab = colnames(x)[i], type = "l")
}


h_omega <- mvspec(x, plot = FALSE, spans = 9)
par(mfcol = c(3, 3), mar = c(4, 5, 2, 1))
for(i in 1:3){
  for(j in 1:3){
    plot(x = h_omega$freq, y = Re(h_omega$fxx[i, j, ]), ylim = range(c(Re(h_omega$fxx[i, j, ]), Im(h_omega$fxx[i, j, ]))),
         type = "l", xlab = quote(omega), ylab = TeX(paste0("$\\hat{h}_{", i, j, "}(\\omega)$")))
    lines(x = h_omega$freq, y = Im(h_omega$fxx[i, j, ]), col = 2)
  }
}

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
colnames(x_hat) <- colnames(x)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    x_hat[t, ] <- x_hat[t, ] + b_j_computed[[as.character(t - i)]] %*% y[k, ]
  }
}

par(mfrow = c(2, 4), mar = c(4, 5, 2, 1))
for(i in 1:8){
  plot(x[, i], xlab = "Time", ylab = colnames(x)[i], type = "l")
  lines(Re(x_hat[, i]), col = 2)
}