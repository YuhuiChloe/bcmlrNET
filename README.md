
# A generalized Bayesian approach to multiple changepoint analysis in dynamic networks

The bcmlrNET function gives the options to 
- detect a fixed number of changepoints when the number of changes is known OR
- detect the number of changes when the number of changes is unknown and used with the bcmlrNET_model_select function that fits bcmlr to first detect the number of changes and generates posterior samples for posterior inference. 

# Main arguments in the bcmlrNET or bcmlrNET_model_select functions: 
- A: a Tt*n*n 3D array with Tt observations of n*n network.
- S: a Tt*p matrix with p summary statistics of each of the Tt networks
- num_CP: The total number of changepoints the user anticipates (Notation in the paper: $k$)
- init: initialization method (The default initialization is to place the changepoints at evenly spaced positions of the data series.)
- prior_kappa: Choose "default" to use the prior (based on segment lengths) we recommend; Enter a function name if you want to use your function that takes in all possible changepoint locations and outputs a value of the joint  PMF of the prior distribution on the changepoints. 
- alpha_f: default to 0.1. (1-alpha_f) is the confidence level o a frequentist confidence internval (Notation in the paper: $\alpha$)
- thinning: default at 1. It leads to how much data is held out (Notation in the paper: $\zeta$)
- min_size: mininum segment length, default to 30. (Notation in the paper: $m$)
- sd_beta: value (variance) on the diagonal of the covariance matrix for regression coefficients, default to 3. (Related notation in the paper:$V_0$, a covaraince matrix)
- max_num_cp: a conservative upper bound of the number of changepoints

# Some outputs of the bcmlrNET_model_select function: 
- Kappa: posterior modes of the changepoint distribution (Notation in the paper: $\kappa$)
- Beta: posterior coefficients (from logistic regression or multinomial logistic regression) (Notation in the paper: $\beta$)
- num_cp_dist: posterior distribution of the number of changepoints
- num_cp: the posterior mode of num_cp_dist = the length of Kappa (Notation in the paper: $L$ (truth), $\hat{k}$ (estimate))
- max_cps: the posterior modes based on the max number of changepoints the user anticipates (Notation in the paper: $k$)
- probs_num_cp: posterior probabilities of all possible numbers of changepoints (up to the length of max_cps). The number with the highest posterior probabilities is the number of changes bcmlr deem there to be. 

