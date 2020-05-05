rm(list=ls()) # clean R environment
cat("\f") # clear the console

source('./utils.R') # source the utils functions
source('./bs_greedy_search.R') # source the bs GS

set.seed(47) # for reproducibility

####### Data Generation
m_obs = 50 # number of observations
var_noise = 1 # noise of the data
data = make_test_Data(m_obs, var_noise) # generate data
bsIterations = c(40, 80, 100) # the amount of iterations we are going to use to compare

# AUC matrix where we are going to save our results for varying number of BS iterations
AUROC_matrix = matrix(nrow = 5, ncol = length(bsIterations)) 
for (j in 1:5) {
  results = vector()
  for (i in 1:length(bsIterations)) {
    ####### Call the BS Greedy Search Algorithm
    currbsIterations = bsIterations[i]  # num BS iterations
    start_time = Sys.time()
    bsResults = bs_greedy_search(data, currbsIterations) # get the bs output list
    end_time = Sys.time()
    elapsed_time = end_time - start_time
    cat('\nElapsed time for bootstrap Greedy Search: ', elapsed_time, '\n')
    
    ###### Evaluate the results
    # sum over all values outputed by the GS
    SCORES_undirected <- matrix(0,11,11)
    for (i in 1:currbsIterations) {
      SCORES_undirected <- (SCORES_undirected + bsResults[[i]] + t(bsResults[[i]]))
    }
    
    # divide by the number of iterations to get the score for a certain edge
    SCORES_undirected = SCORES_undirected/currbsIterations
    
    # Compute the AUROC
    true_Net = make_true_Net()
    true_Net = true_Net + t(true_Net)
    AUROC    = compute_AUROC(SCORES_undirected,true_Net)
    cat('The AUROC of the experiment was:', AUROC)
    results = c(results, AUROC) # append the value of the AUROC
  }
  # append to the results matrix
  AUROC_matrix[j,] = results
}

AUROC_matrix_means = colMeans(AUROC_matrix) # calculate the means of the experiment
AUROC_matrix_std = sqrt(apply(AUROC_matrix,2,var)) # calculate the std
library(ggplot2)
plot(AUROC_matrix_means, type = 'o') # starndard plot
AUROC_df = as.data.frame(AUROC_matrix_means) # to be able to plot with ggplot
mean_minus_std = AUROC_matrix_means - AUROC_matrix_std
mean_plus_std = AUROC_matrix_means + AUROC_matrix_std
p = ggplot(data=AUROC_df, aes(x=c(20, 40, 80), y=AUROC_df$AUROC_matrix_means)) +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin= mean_minus_std, ymax = mean_plus_std)) +
  labs(title = "Mean AUROC vs BS-Iterations") + 
  xlab('Number of BS Iterations') +
  ylab('Mean AUROC') +
  geom_point(size=3, shape=21, fill="white")  # 21 is filled circle
