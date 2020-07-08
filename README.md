# Extensions to Principal Component Analysis (PCA)
This repository provides several R scripts with simulations and applications to real data related to Principal Component Analysis (PCA) and two extensions: Sparse PCA (SPCA) and Dynamic PCA (DPCA). 
This repository is part of my Bachelor Thesis, also available.
# Repository structure
The repository is divided in folders. These folders contain different R scripts grouped by convenience. The folder names are descriptive and, if needed, the data sets used for each script are provided.
# Principal Component Analysis
Principal Component Analysis (PCA) is a well-known technique that takes data from high dimensions and reduce the number of variables trying to maintain most of the explained variance.

First folder contains scripts used to introduce basic statistical concepts and to apply classic PCA numerically.
# Sparse PCA
In the second place, we introduce Sparse PCA (SPCA) which aims to obtain sparse directions on which project the data favouring the analysis interpretability. We perform several numerical experiments and analyze two real data sets, showing the benefits of the extension.

Animation for ridge regression penalty
![Ridge regression penalty](https://github.com/joncsebm/PCA_extensions/blob/master/Sparse%20PCA/SPCA%20theorems%20illustration/ridge.gif)
Animation for the lasso regression penalty
![Lasso regression penalty](https://github.com/joncsebm/PCA_extensions/blob/master/Sparse%20PCA/SPCA%20theorems%20illustration/lasso.gif)
# Dynamic PCA
Lastly, Dynamic PCA (DPCA) is presented as a way of applying PCA to time series. We show how the projected information can be retrieved and we do several numerical examples with simulated and real data.

We also propose an algorithm that uses univariate forecasting techniques for forecasting multivariate time series, showing its benefits and limitations.
