# Function to estimate the tangent point between a circle and a ellipsis

ridge_tangent_b <- function(x0){
  d <- 1000
  i <- 1
  t <- seq(from = -pi / 2, to = 0, length.out = 1000)
  b <- 1
  while(abs(d) > 1e-20){
    x <- 2 * b * cos(t)
    y <- b * sin(t)
    d <- min(sqrt((b * sqrt(2) * cos(t) + b * 0.5 * sqrt(2) * sin(t) + x0[1]) ^ 2 +
                    (-b * sqrt(2) * cos(t) + b * 0.5 * sqrt(2) * sin(t) + x0[2]) ^ 2)) - 2.5
    if(abs(d) > 1e-20){
      if(d > 0){
        b <- b + 1 / i
      }else{
        b <- b - 1 / i
      }
      i <- i + 1
    }
  }
  return(b)
}


library(plotrix)

#ridge penalty

plot(x = -10:10, y = -10:10, type = 'n', asp = 1, xlab = expression(beta[1]), ylab = expression(beta[2]), main = "Ridge regression")
draw.circle(x = 0, y = 0, nv = 1e5, radius = 2.5, col = "slategray1", border = "white")
abline(h = 0, v = 0, col = gray(0.5))
beta_hat <- c(-5, 7.5)
points(x = beta_hat[1], y = beta_hat[2], pch = 16, cex = 0.75)
text(x = beta_hat[1] + 0.5, y = beta_hat[2] + 0.2, labels = expression(hat(beta)), font = 2)
a <- c(2.3, 4.3, 2 * ridge_tangent_b(beta_hat))
b <- a / 2
draw.ellipse(x = rep(beta_hat[1], times = length(a)), y = rep(beta_hat[2], times = length(a)), a = a, b = b, angle = -45, border = 'tomato2',
             nv = 1e2)


#lasso penalty

plot(x = -10:10, y = -10:10, type = 'n', asp = 1, xlab = expression(beta[1]), ylab = expression(beta[2]), main = "Lasso regression")
polygon(x = c(0, -2.5, 0, 2.5), y = c(2.5, 0, -2.5, 0), border = "white", col = "slategray1")
abline(h = 0, v = 0, col = gray(0.5))
beta_hat <- c(-5, 7.5)
points(x = beta_hat[1], y = beta_hat[2], pch = 16, cex = 0.75)
text(x = beta_hat[1] + 0.5, y = beta_hat[2] + 0.2, labels = expression(hat(beta)), font = 2)
a <- c(2.5, 4.5, 5*sqrt(2))
b <- a / 2
draw.ellipse(x = rep(beta_hat[1], times = length(a)), y = rep(beta_hat[2], times = length(a)), a = a, b = b, angle = -45, border = "tomato2",
             nv = 1e2)