# Load data

X <- read.csv("INE_emisiones_sector.csv", sep = ";")
rownames(X) <- X[, 1]
X <- X[, -1]

# Principal components
pca <- princomp(X, cor = TRUE)

# Useful plots

plot(pca, type = "l")
biplot(pca, cex = 0.75)