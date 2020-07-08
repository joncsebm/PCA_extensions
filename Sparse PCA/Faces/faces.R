# Load data

load('faces.RData')
faces <- t(faces)

# # This piece of code reduces resolution averaging matrices of pixels of size 2x2 for faster results (there might be some code that needs to be adapted)

# faces_2 <- matrix(nrow = 2410, ncol = 2016)
# 
# for(i in 1:2410){
#   columna <- 0
#   for(j in 1:2016){
#     faces_2[i, j] <- mean(c(faces[i, (2 * j - 1 + columna):(2 * j + columna)], faces[i, (2 * j + 83 + columna):(2 * j + 84 + columna)]))
#     if(j %% 42 == 0){
#       columna <- columna + 84
#     }
#   }
# }

# Plotted faces from the data set are stored in "n"

n <- vector(mode = "numeric", length = 9)

#Plot nine random faces from the data set

  par(mfrow = c(3, 3), mar = rep(0, 4))
  set.seed(1)
  for(i in 1:9){
    n[i] <- sample(1:length(faces[, 1]), 1)
    image(matrix(rev(faces[n[i], ]), nrow = 84, ncol = 96), col = gray(0:255 / 255), xaxt = "n", yaxt = "n", ann = FALSE)
  }

# Use of the standard PCA function

faces.pca <- prcomp(faces, center = TRUE, scale. = TRUE)

# Use of the SPCA

faces.spca <- list()
library(sparsepca)
for(i in c(1e-4, 1e-3, 1e-2)){
  faces.spca[[as.character(i)]] <- sparsepca::spca(X = faces, k = 10, alpha = i, scale = TRUE)
}

# Show PCA results

  par(mfrow = c(3, 7), mar = rep(0, 4))
  
  for(i in 1:3){
    for(j in seq(from = 0.1, to = 0.9, length.out = 7)){
      image(matrix(rev(faces.pca$center + faces.pca$sdev[i] * quantile(faces.pca$x[, i], probs = j) * faces.pca$rotation[, i]),
                   nrow = 84, ncol = 96), col = gray(0:255 / 255), xaxt = 'n', yaxt = 'n')
    }
  }

# Choice of sparsity to show ("1e-04", "1e-03" or "1e-02")
sparsity <- "1e-04"

# Show SPCA results

  par(mfrow = c(5, 7), mar = rep(0, 4))
  for(i in 1:5){
    for(j in seq(from = 0.1, to = 0.9, length.out = 7)){
      image(matrix(rev(faces.spca[[sparsity]]$center + faces.spca[[sparsity]]$sdev[i] * quantile(faces.spca[[sparsity]]$scores[, i], probs = j) * faces.spca[[sparsity]]$loadings[, i]),
                   nrow = 84, ncol = 96), col = gray(0:255 / 255), xaxt = 'n', yaxt = 'n', xlab = paste0('PC', i))
    }
  }

# Scores densities analysis
plot_quantiles <- function(pca, density, i = 1){
  quantiles <- quantile(x = pca, probs = seq(from = 0.1, to = 0.9, length.out = 7))
  for(j in 1:7){
    segments(x0 = quantiles[j], y0 = 0, y1 = density$y[[which.min(abs(density$x - quantiles[j]))]], col = i, cex = 0.5, lty = 'dashed')
  }
}

# Scores densities

library(latex2exp)

  lim1 <- c(-150, 150)
  lim2 <- lim1
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 1) + 0.1)
  for(i in 1:4){
    density_spca <- list()
    for(j in c("1e-04", "1e-03", "1e-02")){
      density_spca[[j]] <- density(faces.spca[[j]]$scores[, i]) 
    }
    density_pca <- density(faces.pca$x[, i])
    if(i < 3){
      
      plot(density_pca$x, density_pca$y, ylim = range(c(density_pca$y, density_spca[["1e-04"]]$y, density_spca[["1e-03"]]$y, density_spca[["1e-02"]]$y)),
           main = TeX(paste0("PC_{", i, "} scores density")), type = "l", xlab = "", ylab = "", xlim = get(paste0("lim", i)))
      k <- 2
      for(j in c("1e-04", "1e-03", "1e-02")){
        lines(x = density_spca[[j]]$x, y = density_spca[[j]]$y, col = k)
        k <- k + 1
      }
      legend("topleft", title = TeX("$\\lambda$"), legend = c("0", TeX("10^{-4}"), TeX("10^{-3}"), TeX("10^{-2}")), lty = 1, col = 1:4)
      
    } else {
      plot(density_pca$x, density_pca$y, ylim = range(c(density_pca$y, density_spca[["1e-04"]]$y, density_spca[["1e-03"]]$y)),
           main = TeX(paste0("PC_{", i, "} scores density")), type = "l", xlab = "", ylab = "")
      k <- 2
      for(j in c("1e-04", "1e-03")){
        lines(x = density_spca[[j]]$x, y = density_spca[[j]]$y, col = k)
        k <- k + 1
      }
      legend("topleft", title = TeX("$\\lambda$"), legend = c("0", TeX("10^{-4}"), TeX("10^{-3}")), lty = 1, col = 1:3)
    }
  }