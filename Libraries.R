#===========#
# Libraries #
#===========#

library(BayesLogit) # for rpg()
library(LearnBayes) # for computing inverse gamma: rigamma(), alternative to 1/rgamma()
library(MASS) # for mvrnorm()
library(pROC) # for ROC curves and AUC values
options(pROCProgress = list(name = "none")) # Disable pROC progress bar globally
