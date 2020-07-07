# initial parameters
vector_n <- 2 ^ (1:6)
n <- 1e5
mu_1 <- 0
mu_2 <- 0
Mu <- c(mu_1, mu_2)
sigma_1 <- 2
sigma_2 <- 1
rho <- 0.5
Sigma <- rbind(c(sigma_1 ^ 2, rho * sigma_1 * sigma_2), c(rho * sigma_1 * sigma_2, sigma_2 ^ 2))

# for this piece of code no simulations are necessary

##########################################################################################################################################################
  # plots with provided data from analysis
  load("data.RData")
  
  names(data_plot) <- c("Type", "Size", "Angle")
  data_plot_PCA <- data_plot[1:(n * length(vector_n)), 2:3]
  data_plot_lr <- data_plot[(n * length(vector_n) + 1):(2 * n * length(vector_n)), 2:3]
  
  means_PCA_per_size <- cbind(0.775 + 0:(length(vector_n) - 1),
                              aggregate(abs(data_plot_PCA$Angle), by = list(data_plot_PCA$Size), FUN = function(x) mean(as.numeric(as.character(x)))))
  
  names(means_PCA_per_size) <- c('X', 'Size', 'Mean_angle')

  means_lr_per_size <- cbind(1.225 + 0:(length(vector_n) - 1),
                             aggregate(abs(data_plot_lr$Angle), by = list(data_plot_lr$Size), FUN = function(x) mean(as.numeric(as.character(x)))))
  
  names(means_lr_per_size) <- c('X', 'Size', 'Mean_angle')
  
  library('ggplot2')
  ggplot(data = data_plot, aes(x = Size, y = Angle, fill = Type)) + geom_violin() + 
    ylab('Angle (rad)') + xlab('Sample size') + ggtitle('Principal component angle difference between \n population and sample') + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_vline(xintercept = 0.5 + 1:5, col = gray(0.3), linetype = 'dashed') +
    geom_point(data = means_PCA_per_size, aes(x = X, y = Mean_angle), inherit.aes = FALSE) +
    geom_point(data = means_lr_per_size, aes(x = X, y = Mean_angle), inherit.aes = FALSE)

##########################################################################################################################################################

# Code for obtaining data used for the plot 

  # theoretical computation of principal components
  
  lambdas <- c((sigma_1 ^ 2 + sigma_2 ^ 2 + sqrt(sigma_1 ^ 4 + sigma_2 ^ 4 + 2 * sigma_1 ^ 2 * sigma_2 ^ 2 * (2 * rho ^ 2 -1))) / 2,
               (sigma_1 ^ 2 + sigma_2 ^ 2 - sqrt(sigma_1 ^ 4 + sigma_2 ^ 4 + 2 * sigma_1 ^ 2 * sigma_2 ^ 2 * (2 * rho ^ 2 -1))) / 2)
  v1 <- (1 + (rho * sigma_1 * sigma_2 / (lambdas[1] - sigma_1 ^ 2)) ^ 2) ^ (-0.5) * c(rho * sigma_1 * sigma_2 / (lambdas[1] - sigma_1 ^ 2), 1)
  # v2 <- (1 + (rho * sigma_1 * sigma_2 / (lambdas[2] - sigma_1 ^ 2)) ^ 2) ^ (-0.5) * c(rho * sigma_1 * sigma_2 / (lambdas[2] - sigma_1 ^ 2), 1)
  angle_v1 <- atan(v1[2] / v1[1])
  
  set.seed(5)
  
  #simulations of PC's
  theta_nabs <- vector(mode = 'numeric', length = length(vector_n) * n)
  for (j in vector_n){
    for(i in 1:n){
      pca <- princomp(mvtnorm::rmvnorm(n = j, mean = Mu, sigma = Sigma), cor = FALSE)
      theta_nabs[(log2(j) - log2(vector_n[1])) * n + i] <- atan(pca$loadings[2, 1] / pca$loadings[1, 1]) - angle_v1
    }
  }
  
  #-------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #theoretical computation of regression angle
  
  m_r_p <- (rho * sigma_1 * sigma_2) / (sigma_1 ^ 2)
  angle_r_p <- atan(m_r_p)
  v_r_p <- c(1 / sqrt(1 + m_r_p ^ 2), m_r_p / sqrt(1 + m_r_p ^ 2))
  
  set.seed(5)
  
  #simulations of regression line
  
  gamma_nabs <- vector(mode = 'numeric', length = length(vector_n) * n)
  for (j in vector_n){
    for(i in 1:n){
      data <- mvtnorm::rmvnorm(n = j, mean = Mu, sigma = Sigma)
      gamma_nabs[(log2(j) - log2(vector_n[1])) * n + i] <- atan(lm(formula = data[, 2] ~ data[, 1])$coefficients['data[, 1]']) - angle_r_p
    }
  }
  
  data_plot <- data.frame(rep(c("PCA", "Linear regression"),each = n * length(vector_n)),
                          as.factor(rep(rep(vector_n, each = n), times = 2)),
                          c(theta_nabs, gamma_nabs))