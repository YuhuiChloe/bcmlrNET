#===============================================#
# When the number of CPs is known: Use bcmlrNET #
#===============================================#
out = bcmlrNET(A_array, stats_mat, num_CP = 2, latent_dim = 5, init = "even", prior_kappa = "default", alpha_f = 0.1, 
         min_size = 30, sd_beta = 3, thinning = 1, 
         num_iter = 20, num_warmup =10, print_outputs = FALSE, 
         print_progress = TRUE, model_selection = TRUE)
str(out)
# Check posterior mode of changepoints
out$Kappa_mode

#===============================================#
# When the number of CPs is unknown: 
# Use bcmlrNET_model_select that calls bcmlrNET #
#===============================================#
out = bcmlrNET_model_select(A_array, stats_mat, latent_dim = 5, init = "even", thinning = 5, 
                            prior_kappa = "default", alpha_f = 0.1, threshold = 0.5, min_size = 30, 
                            max_num_cp = 4, num_iter = 10, num_warmup = 5, 
                            print_progress = TRUE, print_outputs = FALSE)

str(out)
# Check the distribution of fitted changepoints that correspond to true changepoints
hist(out$Kappa)
