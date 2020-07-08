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
x <- fmri1[, -1]
size <- nrow(x)
p <- ncol(x)

h_omega <- mvspec(x, plot = FALSE, spans = 9)

mirror <- 0

# number of principal components

n <- 1

q <- 1000
predictions <- 72

# c_j and b_j coefficient for series reconstruction

c_j_computed <- list()
for(i in (1 - size - predictions):(2 * (size + predictions))){
  c_j_computed[[as.character(i)]] <- Conj(t(c_j(h = h_omega, j = i, p = p)))
}
b_j_computed <- list()
for(i in (1 - size - predictions):(2 * (size + predictions))){
  b_j_computed[[as.character(i)]] <- c_j(h = h_omega, j = i, p = p)
}

# y series computation

y <- matrix(0, nrow = size, ncol = n)
for(t in 1:size){
  for(i in trunc(-mirror * size):trunc((1 + mirror) * size)){
    k <- ifelse(i <= 0, -i + 1, ifelse(i > size, 2 * size - i, i))
    y[t, ] <- y[t, ] + c_j_computed[[as.character(t - i)]][1:n, ] %*% x[k, ]
  }
}

plot(Re(y), xlab = "Time", ylab = "y", type = "l")
lines(Im(y), col = 2)

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

library(forecast)
re_y_arima <- auto.arima(ts(Re(y)[1:115], frequency = 32), d = 1, D = 1, trace = TRUE)

library(latex2exp)

plot(diff(diff(Re(y), lag = 32)), type = "l", xlab = "Time", ylab = TeX("$\\nabla\\nabla_{32}Re(y_t)$"))

acf(diff(diff(Re(y), lag = 32)))

pacf(diff(diff(Re(y), lag = 32)))

x_hat_sim <- array(0, dim = c(size + predictions, p, q))

for(l in 1:q){
  y_simulated <- simulate(re_y_arima, future = TRUE, nsim = 13 + predictions, seed = l)
  y_simulated <- c(Re(y)[1:115], y_simulated)
  # reconstruction of original series
  for(t in 1:(size + predictions)){
    for(i in 1:(size + predictions)){
        x_hat_sim[t, , l] <- x_hat_sim[t, , l] + b_j_computed[[as.character(t - i)]][, 1] * y_simulated[i]
    }
  }
}

y_simulated <- matrix(nrow = 13 + predictions, ncol = q)
for(i in 1:q){
  y_simulated[, i] <- simulate(re_y_arima, future = TRUE, nsim = 13 + predictions, seed = i)
}
x_int <- 116:(size + predictions)
y_int_sup <- vector(mode = "numeric", length = 13 + predictions)
y_int_inf <- vector(mode = "numeric", length = 13 + predictions)
  for(l in 1:(13 + predictions)){
    y_int_sup[l] <- quantile(y_simulated[l, ], probs = 0.95)
    y_int_inf[l] <- quantile(y_simulated[l, ], probs = 0.05)
  }
plot(Re(y)[1:115], type = "l", xlim = c(0, size + predictions), xlab = "Time", ylab = "y", ylim = range(y_simulated))
for(i in 1:q){
  # points(x_int, y_simulated[, i], col = gray(0.8), cex = 0.25, pch = 16)
  lines(x_int, y_simulated[, i], col = gray(0.8), cex = 0.1)
}
lines(x_int, y_int_sup, col = "deepskyblue4", cex = 0.5)
lines(x_int, y_int_inf, col = "deepskyblue4", cex = 0.5)
lines(1:(size + predictions), c(as.vector(re_y_arima$fitted), as.vector(predict(re_y_arima, n.ahead = 13 + predictions)$pred)), col = 2)
lines(116:(size + predictions), Re(y)[116:(size + predictions)])
abline(v = c(116, size), lty = "dashed", col = c(3, 6))
legend("topleft", legend = c("Original data", "ARIMA fitting", "ARIMA forecast simulations", "90% confidence interval"), lty = 1,
       col = c(1, 2, gray(0.8), "deepskyblue4"))
legend("bottomleft", legend = c("Data used in the ARIMA fitting", "Only forecast from this point"), lty = "dashed", col = c(3, 6))

par(mfcol = c(4, 2), mar = c(2, 4, 2, 2) + 0.1)
for(j in 1:8){
  plot(1, 1, type = "n", ylab = colnames(x)[j], xlim = c(0, size + predictions), ylim = range(Re(x_hat_sim[116:(size + predictions), j, ])))
  for(i in 1:q){
    # points(116:(size + predictions), Re(x_hat_sim[116:(size + predictions), j, i]), col = gray(0.8), pch = 16, cex = 0.25)
    lines(116:(size + predictions), Re(x_hat_sim[116:(size + predictions), j, i]), col = gray(0.9))
  }
  x_int <- 116:(size + predictions)
  y_int_sup <- vector(mode = "numeric", length = length(116:(size + predictions)))
  y_int_inf <- vector(mode = "numeric", length = length(116:(size + predictions)))
  
  for(l in x_int){
    y_int_sup[l - 115] <- quantile(Re(x_hat_sim[l, j, ]), probs = 0.9)
    y_int_inf[l - 115] <- quantile(Re(x_hat_sim[l, j, ]), probs = 0.1)
  }
  lines(x_int, y_int_sup, col = "deepskyblue4", cex = 0.5)
  lines(x_int, y_int_inf, col = "deepskyblue4", cex = 0.5)
  lines(x[, j])
  abline(v = c(116, size), lty = "dashed")
}