#######################################################################################################################################################################
## MS_index_threshold_optimization.R: ##
########################################

## - defines the final MS index threshold

# Parameters: MSindex, threshold, MS_status
#       MSindex: the named vector (according to the sample IDs) of the difference of two cosine correlations (sample from two centroids), increasing sorted
#       threshold: the vector of thresholds that need to be tested
#       MS_status: named factor of MSI status (two levels: MSI/MSS); ordering of samples is the same as in tested dataset

#######################################################################################################################################################################


MS_index_thres_optim <- function(MSindex, threshold, MS_status)
{
  
  if(sum(is.na(MSindex))==0)
  {
    TPs <- sapply(threshold,FUN=function(x) sum(MS_status[match(names(MSindex),names(MS_status))][MSindex > x]=="MSI"))
    FPs <- sapply(threshold,FUN=function(x) sum(MS_status[match(names(MSindex),names(MS_status))][MSindex > x]=="MSS"))
    
    FNs <- sapply(threshold,FUN=function(x) sum(MS_status[match(names(MSindex),names(MS_status))][MSindex <= x]=="MSI"))
    TNs <- sapply(threshold,FUN=function(x) sum(MS_status[match(names(MSindex),names(MS_status))][MSindex <= x]=="MSS"))
    
    n <- length(MSindex)
    sensitivity_vect <- TPs/(TPs+FNs)
    specificity_vect <- TNs/(TNs+FPs)
    
    optimization_res <- list()
    optim_index <- which.max(sensitivity_vect + specificity_vect)
    optimization_res$optim_threshold <- threshold[optim_index]

    optimization_res
    
  } else 
  {
    NA
  }
  
}

