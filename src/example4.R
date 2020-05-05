rm(list=ls()) # clean R environment
cat("\f") # clear the console

source('./utils.R') # source the utils functions
source('./bs_greedy_search.R') # source the bs GS
source('./MCMC.R') # source the MCMC script

set.seed(47) # for reproducibility

####### Data Generation
m_obs = c(20, 50, 100) # number of observations
var_noise = 1 # noise of the data
dataList = list() # generate a list of data
data.idx = 0
for (idx in 1:3) {
  for (jdx in 1:5) {
    data.idx = data.idx + 1
    dataList[[data.idx]] = make_test_Data(m_obs[idx], var_noise)  
  }
}

n_nodes = 11  #  no. of network nodes
start_incidence = matrix(0, n_nodes,n_nodes)
SCORES_undirected = matrix(0,11,11) # declare an empty edge score matrix
true_Net = make_true_Net() # get the true network structure
true_Net = true_Net + t(true_Net) # get also the simmetric true edges
SCORES_undirected = matrix(0,11,11) # declare an empty edge score matrix

# MCMC params
iterations  = 20000 # Total no. of MCMC iterations per run
step_save   = 100 # Thin-out by a factor of 100

# BS GS params
bsIterations = 40 # num BS iterations

mcmc_auroc = vector()
greedy_auroc = vector()
bs_gs_auroc = vector()
for (it in 1:15){
  # get the current data
  curr_data = dataList[[it]]
  
  # Run the MCMC algo
  sim = strMCMC(curr_data, start_incidence, iterations, step_save)
  n   = length(sim[[1]]) # Get length of the chain
  for (i in 1:n){
    # Sum the number of apperances on each score
    SCORES_undirected = (SCORES_undirected + sim[[1]][[i]] + t(sim[[1]][[i]]))
  }
  SCORES_undirected = SCORES_undirected/n # divivde by n to get the posterior edge score matrix
  AUROC = compute_AUROC(SCORES_undirected, true_Net) # calculate the AUC
  mcmc_auroc = c(mcmc_auroc, AUROC)
  cat('Finished MCMC!\n')
  
  # Run the GS algo
  greedyRes = greedy_search(curr_data, start_incidence)
  # Score the algo
  greedy_undirected = matrix(0,11,11)
  greedy_undirected = greedyRes + t(greedyRes)
  curr_greedy_auroc = compute_AUROC(greedy_undirected, true_Net)
  greedy_auroc = c(greedy_auroc, curr_greedy_auroc)
  cat('Finished Vanilla GS!\n\n')
  
  # Run the BS GS algo
  bsResults = bs_greedy_search(curr_data, bsIterations) # get the bs output list
  SCORES_undirected <- matrix(0,11,11)
  for (i in 1:bsIterations) {
    SCORES_undirected <- (SCORES_undirected + bsResults[[i]] + t(bsResults[[i]]))
  }
  
  # divide by the number of iterations to get the score for a certain edge
  SCORES_undirected = SCORES_undirected/bsIterations
  
  # Compute the AUROC
  bsAUROC    = compute_AUROC(SCORES_undirected,true_Net)
  bs_gs_auroc = c(bs_gs_auroc, bsAUROC) # append to the vector of results
}
