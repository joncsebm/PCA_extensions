# Load data
load("MNIST.RData")

# Function to show a given image in the data set

show_digit <- function(vec, col = gray(12:1 / 12), ...){
  image(matrix(vec, nrow = 28)[, 28:1], col = col, xaxt = "n", yaxt = "n", ann = FALSE, ...)
}

# Function to show a fictional score within a PC

show_pca <- function(pca, k = 1, lambda = 0, col = gray(12:1 / 12), ...){
  reconstruction <- pca$center + lambda * pca$loadings[, k]
  reconstruction <- pmin(pmax(reconstruction, 0), 1)
  image(matrix(reconstruction, nrow = 28)[, 28:1],
        col = col, zlim = c(0, 1), xaxt = "n", yaxt = "n", ann = FALSE, ...)
}

# Code for PCA and SPCA computation (not required, results are provided)

# Classic PCA analysis

pca_MNIST <- list()
for(i in c(2, 5, 8)){
    pca_MNIST[[paste0("n=", i)]] <- princomp(x = MNIST$x[MNIST$labels == i, ], cor = FALSE)
}

# Covariance estimation of the data set to speed up elasticnet computations

covariance_matrices <- list()
for(i in c(2, 5, 8)){
    covariance_matrices[[paste0('n=', i)]] <- cov(x = MNIST$x[MNIST$labels == i, ])
}

# SCPA from elasticnet (for digits 2, 5 and 8; for others modify the code)

library(elasticnet)
spca_MNIST <- list()
spca_MNIST[["n=2"]] <- spca(x = covariance_matrices[['n=2']], K = 1, para = 1e-1, type = 'Gram', sparse = 'penalty')
spca_MNIST[["n=2"]]$center <- pca_MNIST[['n=2']]$center
spca_MNIST[["n=5"]] <- spca(x = covariance_matrices[['n=5']], K = 1, para = 1e-1, type = 'Gram', sparse = 'penalty')
spca_MNIST[["n=5"]]$center <- pca_MNIST[['n=5']]$center
spca_MNIST[["n=8"]] <- spca(x = covariance_matrices[['n=8']], K = 1, para = 1e-1, type = 'Gram', sparse = 'penalty')
spca_MNIST[["n=8"]]$center <- pca_MNIST[['n=8']]$center

# Show an example for each digit

  par(mfrow = c(2, 5), mar = rep(0, 4))
  for(i in c(1:6, 8, 10, 14, 18)){
    show_digit(vec = MNIST$x[i, ])
  }

# Computation results are provided
  
load('prin_comp_analysis.RData')

  # Center of the data, that is, mean of each digit

  par(mfrow = c(2, 5), mar = rep(0, 4))
  for(i in 0:9){
    show_digit(vec = colMeans(MNIST$x[MNIST$labels == as.character(i), ]))
  }
  
  # show in layout different values of lambda for the k-th PC
  
  par(mfcol = c(3, 9), mar = rep(0, 4))
  for(i in seq(from = -10, to = 10, by = 2.5)){
    show_pca(pca_MNIST[["n=2"]], lambda = i, k = 1)
    show_pca(pca_MNIST[["n=5"]], lambda = i, k = 1)
    show_pca(pca_MNIST[["n=8"]], lambda = i, k = 1)
  }
  
  # show in layout different values of lambda for the SPC

  par(mfcol = c(3, 9), mar = rep(0, 4))
  for(i in seq(from = -10, to = 10, by = 2.5)){
    show_pca(spca_MNIST[['n=2']], lambda = i, k = 1)
    show_pca(spca_MNIST[['n=5']], lambda = i, k = 1)
    show_pca(spca_MNIST[['n=8']], lambda = i, k = 1)
  }

# # Animations
# 
# # Manipulate plot in function of lambda
# library(manipulate)
# 
# # first principal components directions
# 
# layout(mat = matrix(1:3, nrow = 1, ncol = 3, byrow = TRUE))
# manipulate({
#   for(i in c(2, 5, 8)){
#     show_pca(pca_MNIST[[paste('n=', i, sep = '')]], k = 3, lambda = lambda)
#   }
# }, lambda = slider(-10, 10, initial = 0))
# 
# library('magick')
# k_1 <- image_graph(width = 1600, height = 1200, res = 96)
# for(i in -10:10){
#   layout(mat = matrix(1:3, nrow = 1, ncol = 3, byrow = TRUE))
#   for(j in c(2, 5, 8)){
#     show_pca(spca_MNIST[[paste('n=', j, sep = '')]], k = 1, lambda = i)
#   }
# }
# dev.off()
# image_write(image_animate(k_1, fps = 4, optimize = TRUE), 'k_1.gif')
