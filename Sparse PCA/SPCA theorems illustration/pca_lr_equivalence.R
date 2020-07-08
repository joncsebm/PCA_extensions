# Initial parameters

l <- 250
mu <- c(0, 0)
sigma_x <- 1
sigma_y <- 1
rho_xy <- 0.8
Sigma <- matrix(data = c(sigma_x ^ 2, rep(sigma_x * sigma_y * rho_xy, times = 2), sigma_y ^ 2), nrow = 2, ncol = 2)

# Simulate data

library(mvtnorm)
set.seed(5)
data <- as.data.frame(rmvnorm(n = l, mean = mu, sigma = Sigma))
data$V1 <- data$V1 - mean(data$V1)
data$V2 <- data$V2 - mean(data$V2)

library(latex2exp)

visualizeFitLm3D <- function(data, lambda = 0, nGrid = 10, theta = 60, phi = 15, alpha = 0, decimals = 1){
  
  # choose between n = 1 (PC1) and n = 2 (PC2)
  data <- cbind(data, princomp(x = data)$scores[, 1], princomp(x = data)$scores[, 2])
  names(data) <- c("X", "Y", "Z1", "Z2")
  
  # PC directions
  loadings <- princomp(x = data[, c("X", "Y")])$loadings
  var_pc <- princomp(x = data[, c("X", "Y")])$sdev ^ 2
  
  # Estimate lm
  X <- as.matrix(data[, c("X", "Y")])
  Y1 <- data$Z1
  Y2 <- data$Z2
  library("glmnet")
  # beta_hat_lambda_1 <- solve(crossprod(X) + diag(lambda, nrow = 2, ncol = 2)) %*% t(X) %*% Y
  beta_hat_lambda_1 <- glmnet(x = X, y = Y1, lambda = lambda, alpha = alpha, thresh = 1e-20, standardize = FALSE)$beta
  beta_hat_lambda_2 <- glmnet(x = X, y = Y2, lambda = lambda, alpha = alpha, thresh = 1e-20, standardize = FALSE)$beta
  if(beta_hat_lambda_2[1] == 0 | beta_hat_lambda_2[2] == 0){
    plano_horizontal <- TRUE
  } else {
    plano_horizontal <- FALSE
  }
  # lrm <- lm(Z ~ 0 + X + Y, data = data)
  
  # Create and plot initial XYZ-grid
  sdX <- sd(data$X)
  sdY <- sd(data$Y)
  sdZ <- max(sd(data$Z1), sd(data$Z2))
  gX <- seq(min(data$X) - sdX, max(data$X) + sdX, length = nGrid)
  gY <- seq(min(data$Y) - sdY, max(data$Y) + sdY, length = nGrid)
  gZ <- seq(min(c(data$Z1, data$Z2)) - sdZ, max(c(data$Z1, data$Z2)) + sdZ, length = nGrid)
  
  # Plot data
  require(plot3D)
  panelFirst <- function(pmat) {
    # Plot points projection in the basal plane
    XY <- trans3D(data$X, data$Y, rep(gZ[1], length(data$X)), pmat = pmat)
    scatter2D(XY$x, XY$y, pch = 16, cex = 0.5, add = TRUE,
              colkey = FALSE, col = "darkslategrey")
  }
  
  library(latex2exp)
  
  # Plot
  
  # gridMat1 <- scatter3D(data$X, data$Y, data$Z1, pch = 16, theta = theta,
  #                       phi = phi, bty = "g", axes = FALSE, colkey = FALSE,
  #                       col = "indianred2", xlim = range(gX), ylim = range(gY),
  #                       zlim = range(gZ), panel.first = panelFirst,
  #                       nticks = nGrid, cex = 0.75, scale = TRUE, main = TeX(paste0("$\\lambda = ", sprintf(paste0("%.", decimals, "f"), lambda), "$")))
  gridMat1 <- scatter3D(data$X, data$Y, data$Z1, pch = 16, theta = theta,
                        phi = phi, bty = "g", axes = FALSE, colkey = FALSE,
                        col = "indianred2", xlim = range(gX), ylim = range(gY),
                        zlim = range(gZ), panel.first = panelFirst,
                        nticks = nGrid, cex = 0.75, scale = TRUE)
  gridMat2 <- scatter3D(data$X, data$Y, data$Z2, pch = 16, theta = theta,
                        phi = phi, bty = "g", axes = FALSE, colkey = FALSE,
                        col = "dodgerblue2", xlim = range(gX), ylim = range(gY),
                        zlim = range(gZ), cex = 0.75, add = TRUE)
  # text(x = trans3d(median(gX), gY[1], gZ[1], gridMat1), labels = "x", pos = 1)
  # text(x = trans3d(gX[1], median(gY), gZ[1], gridMat1), labels = "y", pos = 2)
  # text(x = trans3d(gX[1], gY[nGrid], median(gZ), gridMat1), labels = "z",pos = 2)
  
  # Regression plane
  grid_plane <- 20
  x_plane <- seq(from = min(data$X), to = max(data$X), length.out = grid_plane)
  y_plane <- seq(from = min(data$Y), to = max(data$Y), length.out = grid_plane)
  mesh <- mesh(x = x_plane, y = y_plane)
  z1_plane <- beta_hat_lambda_1["X", ] * mesh$x + beta_hat_lambda_1["Y", ] * mesh$y
  z2_plane <- beta_hat_lambda_2["X", ] * mesh$x + beta_hat_lambda_2["Y", ] * mesh$y
  surf3D(x = mesh$x, y = mesh$y, z = z1_plane, col = "indianred", alpha = 0.25, border = NA, add = TRUE)
  if(!plano_horizontal){
    surf3D(x = mesh$x, y = mesh$y, z = z2_plane, col = "dodgerblue", alpha = 0.1, border = NA, add = TRUE)
  } else {
    rect3D(x0 = min(data$X), y0 = min(data$Y), z0 = 0, x1 = max(data$X), y1 = max(data$Y), col = "dodgerblue",
           add = TRUE, alpha = 0.25)
  }
  
  norm_end_1 <- c(beta_hat_lambda_1["X", ], beta_hat_lambda_1["Y", ], -1) / norm(c(beta_hat_lambda_1["X", ], beta_hat_lambda_1["Y", ], -1), type = "2")
  arrows3D(x0 = 0, y0 = 0, z0 = 0, x1 = norm_end_1[1], y1 = norm_end_1[2], z1 = norm_end_1[3], add = TRUE, col = "indianred2")
  arrows3D(x0 = 0, y0 = 0, z0 = gZ[1], x1 = norm_end_1[1] / norm(as.matrix(norm_end_1[1:2]), type = "2"), y1 = norm_end_1[2] / norm(as.matrix(norm_end_1[1:2]), type = "2"), col = 1, add = TRUE)
  segments3D(x0 = norm_end_1[1], y0 = norm_end_1[2], z0 = gZ[1], z1 = norm_end_1[3], add = TRUE, col = "indianred1")
  
  norm_end_2 <- c(beta_hat_lambda_2["X", ], beta_hat_lambda_2["Y", ], -1) / norm(c(beta_hat_lambda_2["X", ], beta_hat_lambda_2["Y", ], -1), type = "2")
  arrows3D(x0 = 0, y0 = 0, z0 = 0, x1 = norm_end_2[1], y1 = norm_end_2[2], z1 = norm_end_2[3], add = TRUE, col = "dodgerblue2")
  arrows3D(x0 = 0, y0 = 0, z0 = gZ[1], x1 = norm_end_2[1] / norm(as.matrix(norm_end_2[1:2]), type = "2"), y1 = norm_end_2[2] / norm(as.matrix(norm_end_2[1:2]), type = "2"), col = 1, add = TRUE)
  segments3D(x0 = norm_end_2[1], y0 = norm_end_2[2], z0 = gZ[1], z1 = norm_end_2[3], add = TRUE, col = "dodgerblue1")
  
  # PC1
  arrows3D(x0 = 0, y0 = 0, z0 = gZ[1], x1 = 3 * loadings[1, 1], y1 = 3 * loadings[2, 1], add = TRUE, col = "red")
  # PC2
  arrows3D(x0 = 0, y0 = 0, z0 = gZ[1], x1 = 3 * loadings[1, 2], y1 = 3 * loadings[2, 2], add = TRUE, col = "blue")
}

par(mar = rep(0, 4), oma = c(1, 0, 1, 0))

visualizeFitLm3D(data = data, lambda = 0, theta = 60, phi = 20)

visualizeFitLm3D(data = data, lambda = 0.1, theta = 60, phi = 20)

visualizeFitLm3D(data = data, lambda = 0.2, theta = 60, phi = 20)

visualizeFitLm3D(data = data, lambda = 0.3, theta = 60, phi = 20)

# GIF simulations

library(magick)

par(mar = rep(0, 4), oma = c(1, 0, 1, 0))
graph_ridge <- image_graph()
for(i in c(seq(from = 0, to = 0.8, by = 0.05), seq(from = 1, to = 12, by = 0.5))){
  visualizeFitLm3D(data = data, alpha = 0, lambda = i, decimals = 2)
  segments(x0 = -0.46, y0 = -0.5, y1 = 0.28, col = "gray")
  segments(x0 = -0.47, x1 = -0.45, y0 = -0.5 + i / 12 * 0.78)
  
  text(x = rep(-0.4, 3), y = c(-0.5, 0.28, -0.5 + i / 12 * 0.78),
       labels = c("0.0", "12.0", sprintf("%.1f", i)))
}
dev.off()

image_write(image_animate(graph_ridge, fps = 10, optimize = TRUE), "ridge.gif")

par(mar = rep(0, 4), oma = c(1, 0, 1, 0))
graph_lasso <- image_graph()
for(i in c(seq(from = 0, to = 0.2, by = 0.01), 0.1 * 2:14)){
  visualizeFitLm3D(data = data, alpha = 1, lambda = i, decimals = 2)
  segments(x0 = -0.46, y0 = -0.5, y1 = 0.28, col = "gray")
  segments(x0 = -0.47, x1 = -0.45, y0 = -0.5 + i / 1.4 * 0.78)
  
  text(x = rep(-0.4, 3), y = c(-0.5, 0.28, -0.5 + i / 1.4 * 0.78),
       labels = c("0.00", "1.40", sprintf("%.2f", i)))
}
dev.off()

image_write(image_animate(graph_lasso, fps = 10, optimize = TRUE), "lasso.gif")