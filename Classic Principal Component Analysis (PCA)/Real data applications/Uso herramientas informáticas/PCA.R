#Load data

X <- read.csv("INE_herramientas_informaticas.csv", sep = ",")
rownames(X) <- X[, 1]
X <- X[, -1]

#Principal component analysis
pca <- princomp(X, cor = FALSE)

# standard deviations of each PC
barplot(pca$sdev, col = gray(0.95), ylab = 'Standard deviation')

# variances of each PC
barplot(pca$sdev ^ 2, col = gray(0.95), ylab = 'Variance')

# cumulative explained variance computation

cum_var <- vector(mode = 'numeric', length = 8)
cum_var[1] <- pca$sdev[1] ^ 2
for(i in 2:8){
  cum_var[i] <- cum_var[i - 1] + pca$sdev[i] ^ 2
}
total_var <- 18 / 19 * var(X[, 1])
for(i in 2:8){
  total_var <- total_var + 18 / 19 * var(X[, i])
}

# % of cumulative explained variance per PC

plot(x = -5, y = -5, xlim = c(1, 8), ylim = c(60, 100), xlab = 'Component', ylab = 'Cumulative % of explained variance')
points(x = 1:8, y = 100 * cum_var / total_var, pch = 16)

# # Heuristic rule to choose the number of PCs
# 
# rect(xleft = 0, xright = 9, ybottom = 70, ytop = 100, density = 20, col = gray(0.75), xlim = c(1, 8), ylim = c(40, 100))
# rect(xleft = 0, xright = 9, ybottom = 80, ytop = 100, density = 20, angle = -45, col = gray(0.75))

# biplot
biplot(pca)

# scores plot
plot(x = pca$scores[, 1], y = pca$scores[, 2], xlab = 'PC1', ylab = 'PC2', pch = 16, asp = 1, cex = 0.6)
