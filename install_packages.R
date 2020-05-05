install.packages("devtools") # for installing stuff
install.packages("svMisc") # Progress bar
install.packages('expm') # package for linear algebra stuff
install.packages('magrittr') # package for the %>% operator
install.packages('dplyr')
install.packages('purrr') # more dependencies
install.packages('ggplot2') # for graphs
# Commands to install the LIMMA package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma") # install the limma package
BiocManager::install("ROC") # install the ROC package
