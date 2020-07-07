# Load data

X <- read.csv("INE_consumos_energeticos.csv", sep = ";")
rownames(X) <- X[, 1]
X <- X[, -1]

# Principal component analysis
pca <- princomp(X, cor = TRUE)

# Useful plots

plot(pca, type = "l")
biplot(pca, cex = 0.6)