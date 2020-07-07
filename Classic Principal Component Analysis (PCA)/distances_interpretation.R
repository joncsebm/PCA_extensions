# Initial parameters

n <- 30
sigma1 <- 1
sigma2 <- 1
rho <- 0.7
Sigma <- rbind(c(sigma1 ^ 2, rho * sigma1 * sigma2), c(rho * sigma1 * sigma2, sigma2 ^ 2))

# Generation of variables

set.seed(5)
x <- mvtnorm::rmvnorm(n = n, mean = c(0, 0), sigma = Sigma)
x[, 1] <- x[, 1] - mean(x[, 1])
x[, 2] <- x[, 2] - mean(x[, 2])

# Principal Components computation
pca <- princomp(x)

# Graphical representation

basic_plot <- function(x, pca){
  p0 <- c(mean(x[, 1]), mean(x[, 2]))
  v1 <- pca$sdev[1] ^ 2 * pca$loadings[, 1]
  v2 <- pca$sdev[2] ^ 2 * pca$loadings[, 2]
  plot(x, asp = 1, pch = 20, cex = 0.6, xlim = c(-3 * sigma1, 3 * sigma1), ylim = c(-3 * sigma2, 3 * sigma2),
       xlab = 'x1', ylab = 'x2')
  # x_1-axis
  arrows(x0 = p0[1], y0 = p0[2], x1 = p0[1] + 2 * sigma1, length = 0.1, lwd = 2)
  text(x = p0[1] + 2 * sigma1, y = p0[2] - 0.2 * sigma1, label = 'x1')
  # x_2-axis
  arrows(x0 = p0[1], y0 = p0[2], y1 = p0[2] + 2 * sigma2, length = 0.1, lwd = 2)
  text(x = p0[1] - 0.3 * sigma2, y = p0[2] + 2 * sigma2, label = 'x2')
}

principal_components_plot <- function(x, pca){
  p0 <- c(mean(x[, 1]), mean(x[, 2]))
  # PC's directions
  angles <- c(atan(pca$loadings[2, 1]/pca$loadings[1, 1]), atan(pca$loadings[2, 2]/pca$loadings[1, 2]))
  # PC1 direction line
  x_line <- seq(from = -5 * sigma1, to = 5 * sigma1, length.out = 1000)
  y_line <- tan(angles[1]) * x_line + p0[2] - p0[1] * tan(angles[1])
  lines(x = x_line, y = y_line, col = 'blue')
  legend('bottomright', legend = 'PC1', col = 'blue', lty = 1, title = 'Principal component')
}

principal_components_axis <- function(x, pca){
  p0 <- c(mean(x[, 1]), mean(x[, 2]))
  angles <- c(atan(pca$loadings[2, 1]/pca$loadings[1, 1]), atan(pca$loadings[2, 2]/pca$loadings[1, 2]))
  # PC1-axis
  arrows(x0 = p0[1], y0 = p0[2], x1 = p0[1] + 2.5 * sigma1 * norm(pca$loadings[, 1], type = '2') * cos(angles[1]),
         y1 = p0[2] + 2.5 * sigma1 * norm(pca$loadings[, 1], type = '2') * sin(angles[1]),
         length = 0.1, col = 2, lwd = 1.75)
  text(x = p0[1] + 2.5 * sigma1 * norm(pca$loadings[, 1], type = '2') * cos(angles[1]) - 0.4 * sigma1,
       y =  p0[2] + 2.5 * sigma1 * norm(pca$loadings[, 1], type = '2') * sin(angles[1]) - 0.2 * sigma1,
       label = 'PC1', col = 2, srt = angles[1] * 360 / (2 * pi))
  # PC2-axis
  arrows(x0 = p0[1], y0 = p0[2], x1 = p0[1] + 2.5 * sigma1 * norm(pca$loadings[, 2], type = '2') * cos(angles[2]),
         y1 = p0[2] + 2.5 * sigma1 * norm(pca$loadings[, 2], type = '2') * sin(angles[2]),
         length = 0.1, col = 3, lwd = 1.75)
  text(x = p0[1] + 2.5 * sigma1 * norm(pca$loadings[, 2], type = '2') * cos(angles[2]) - 0.3 * sigma1,
       y =  p0[2] + 2.5 * sigma1 * norm(pca$loadings[, 2], type = '2') * sin(angles[2]) + 0.05 * sigma1,
       label = 'PC2', col = 3, srt = angles[2] * 360 / (2 * pi))
}

# Regression lines

plot_regression_y_x <- function(x){
  coefs1 <- coef(lm(x[, 2] ~ x[, 1]))
  abline(coef = coefs1, col = 5)
}

plot_regression_x_y <- function(x){
  coefs2 <- coef(lm(x[, 1] ~ x[, 2]))
  abline(coef = c(-coefs2[1] / coefs2[2], 1 / coefs2[2]), col = 6)
}

regression_legend <- function(choice){
  if(choice == 1){
    legend('topleft', legend = 'x2 ~ x1', col = 5, lty = 1, title = 'Regression line')
  }else if(choice == 2){
    legend('topleft', legend = 'x1 ~ x2', col = 6, lty = 1, title = 'Regression line')
  }else{
  legend('topleft', legend = c('x1 ~ x2', 'x2 ~ x1'), col = c(5, 6), lty = 1, title = 'Regression lines')
  }
}

projections_principal_components <- function(x, pca){
  # Projections (Principal components)
  PC1 <- c(pca$loadings[1, 1], pca$loadings[2, 1])
  for(i in 1:n){
    segments(x0 = (PC1 * pca$scores[i, 1])[1],
             y0 = (PC1 * pca$scores[i, 1])[2],
             x1 = x[i, 1], y1 = x[i, 2], lwd = 0.5, lty = 4, col = 'gray')
  }
}

projections_y_x <- function(x){
  # Projections (Regression y~x)
  intercept <- lm(x[, 2] ~ x[, 1])$coefficients[1]
  slope <- lm(x[, 2] ~ x[, 1])$coefficients[2]
  for(i in 1:n){
    segments(x0 = x[i, 1], y0 = intercept + slope * x[i, 1],
             x1 = x[i, 1], y1 = x[i, 2], lwd = 0.5, lty = 4, col = 'gray')
  }
}

projections_x_y <- function(x){
  # Projections (Regression x~y)
  intercept <- lm(x[, 1] ~ x[, 2])$coefficients[1]
  slope <- lm(x[, 1] ~ x[, 2])$coefficients[2]
  for(i in 1:n){
    segments(x0 = x[i, 1], y0 = x[i, 2],
             x1 = intercept + slope * x[i, 2], y1 = x[i, 2], lwd = 0.5, lty = 4, col = 'gray')
  }
}

# PC's with projections

basic_plot(x, pca)
principal_components_plot(x, pca)
projections_principal_components(x, pca)

# Regression line 1 with projections

basic_plot(x, pca)
plot_regression_y_x(x)
regression_legend(1)
projections_y_x(x)


# Regression line 1 with projections

basic_plot(x, pca)
plot_regression_x_y(x)
regression_legend(2)
projections_x_y(x)


# Final plot

basic_plot(x, pca)
principal_components_plot(x, pca)
principal_components_axis(x, pca)
plot_regression_y_x(x)
plot_regression_x_y(x)
regression_legend(3)