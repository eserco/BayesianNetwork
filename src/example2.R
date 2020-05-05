rm(list=ls()) # clean R environment
cat("\f") # clear the console

source('./utils.R') # source the utils functions
source('./bs_greedy_search.R') # source the bs GS

set.seed(47) # for reproducibility

####### Data Generation
dataList = list()
idx = 0
for (i in seq(20,100,10)) {
  idx = idx + 1
  m_obs = i
  var_noise = 1 # noise of the data
  currData = make_test_Data(m_obs, var_noise) # generate data
  
  dataList[[idx]] = currData # append to the current dataList
}

# Matrix of results for now just 3 pts for mean + std
AUROC_matrix = matrix(nrow = 3, ncol = length(dataList)) 
for (j in 1:3){
  n_AUROC = vector()
  for (i in 1:length(dataList)) {
    data = dataList[[i]] # Get the current data
    ####### Call the BS Greedy Search Algorithm
    bsIterations = 40 # num BS iterations
    start_time = Sys.time()
    bsResults = bs_greedy_search(data, bsIterations) # get the bs output list
    end_time = Sys.time()
    elapsed_time = end_time - start_time
    cat('\nElapsed time for bootstrap Greedy Search: ', elapsed_time, '\n')
    
    ###### Evaluate the results
    # sum over all values outputed by the GS
    SCORES_undirected <- matrix(0,11,11)
    for (i in 1:bsIterations) {
      SCORES_undirected <- (SCORES_undirected + bsResults[[i]] + t(bsResults[[i]]))
    }
    
    # divide by the number of iterations to get the score for a certain edge
    SCORES_undirected = SCORES_undirected/bsIterations
    
    # Compute the AUROC
    true_Net = make_true_Net()
    true_Net = true_Net + t(true_Net)
    AUROC    = compute_AUROC(SCORES_undirected,true_Net)
    n_AUROC = c(n_AUROC, AUROC) # append to the vector of results
  }
  AUROC_matrix[j,] = n_AUROC # append to the matrix of results
}

# Plot the results
AUROC_matrix_means = colMeans(AUROC_matrix) # calculate the means of the experiment
AUROC_matrix_std = sqrt(apply(AUROC_matrix,2,var)) # calculate the std
library(ggplot2)
plot(AUROC_matrix_means, type = 'o') # starndard plot
AUROC_df = as.data.frame(AUROC_matrix_means) # to be able to plot with ggplot
mean_minus_std = AUROC_matrix_means - AUROC_matrix_std
mean_plus_std = AUROC_matrix_means + AUROC_matrix_std
p = ggplot(data=AUROC_df, aes(x=seq(20,100,10), y=AUROC_df$AUROC_matrix_means)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin= mean_minus_std, ymax = mean_plus_std)) +
  labs(title = "Mean AUROC vs Sample Size") + 
  xlab('Sample Size') +
  ylab('Mean AUROC') +
  geom_point(size=3, shape=21, fill="white")  # 21 is filled circle
p