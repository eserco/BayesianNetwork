source('./greedy_search/Greedysearch.R', chdir = TRUE) # source the vanilla greedy search 

library(magrittr) # for the %>% operator 
library(dplyr)
library(purrr)
require(svMisc) # for the progress bar

bs_greedy_search = function(data, numIter) {
  '
    Will run the Greedy Search Bootstrap algorithm
    
    @param data matrix with the rows: num_nodes cols: num obs
    @param numIter number of iterations that the bootstrap sample will run
    
    @return bsResults list containing the best incidence matrix on every 
    iteration of the GS algorithm
  '
  n_nodes = nrow(data)  #  no. of network nodes
  start_incidence = matrix(0, n_nodes, n_nodes) # Generate an empty incidence matrix
  
  bootTrials = numIter # the amount of bs trials
  pbTick = (1 / bootTrials) * 100
  bsResults = list()
  progress(0) # inital tick of the pb
  for(i in 1:bootTrials) {
    # Sample from your data with replacement
    index       = sample(ncol(data), replace = TRUE)
    bsData      = data[1:n_nodes, index] # get the data from the bs sampled indexes
    
    currBsResults   = greedy_search(bsData, start_incidence) # run a greedy search on the bs sample
    bsResults[[i]] = currBsResults # append to the bsResults list
    progress(i * pbTick) # progress bar
  }
  
  return(bsResults) # return the bs list
}

