#######################################################################################################################################################################
## nearest_centroid_cosine.R: ##
################################

## - this function performs the neares centroid classification according to the cosine correlateion

# Parameters: data_train/data_test, signature, MS_status_train/MS_status_test, class_of_interest
#     data_train/data_test - matrix of expression values; samples in columns, features in rows
#     signature - character vector with existing signature of entrez gene IDs (e.g.: "entrezXXXX")
#     MS_status_train/MS_status_test - named factor of MSI status (two levels: MSI/MSS)
#     class_of_interest: a character specifying the positive group

#######################################################################################################################################################################


nearest_centroid_cosine <- function(data_train, data_test, signature, MS_status_train, MS_status_test, class_of_interest)

{

  library(lsa)
  possibilities <- levels(MS_status_train)
  negative_class <- possibilities[possibilities!=class_of_interest]
  
  ## centroids of two groups definition
  centroid_mat <- matrix(NA,nrow=length(intersect(signature, rownames(data_train))),ncol=length(possibilities))
  rownames(centroid_mat) <- intersect(signature, rownames(data_train))
  colnames(centroid_mat) <- possibilities
  
  for(i in c(1:length(possibilities)))
  {
    sam <- names(MS_status_train)[MS_status_train==possibilities[i]]
    tam <- intersect(signature, rownames(data_train))
    centroid_mat[,i] <- apply(data_train[tam,sam],1,mean)
  }
  
  ## cosine correlation of each sample to centroids and assignment  
  data_sig_train <- data_train[intersect(signature, rownames(data_train)),]
  data_sig_test <- data_test[intersect(signature, rownames(data_test)),]
  
  correlations_train <- rep(NA, length(possibilities))
  assignments_train <- rep(NA, dim(data_sig_train)[2])
  MSindex_train <- rep(NA, dim(data_sig_train)[2])
  names(MSindex_train) <- colnames(data_train)
  
  correlations_test <- rep(NA, length(possibilities))
  assignments_test <- rep(NA, dim(data_sig_test)[2])
  MSindex_test <- rep(NA, dim(data_sig_test)[2])
  names(MSindex_test) <- colnames(data_test)
  
  ## MS index definition
  for(i in 1:(dim(data_sig_train)[2])){
    for(j in 1:(length(possibilities))){
      correlations_train[j] <- cosine(data_sig_train[,i],centroid_mat[,j])
    }
    MSindex_train[i] <- correlations_train[1] - correlations_train[2] # GROUP OF INTEREST - NEGATIVE GROUP
    assignments_train[i] <- which.max(correlations_train)
  }
  
  for(i in c(1:length(possibilities)))
  {
    assignments_train[assignments_train==i] <- possibilities[i]
  }
  
  MSindex_train_sorted <- sort(MSindex_train)
  
  ## optimal threshold identification
  threshold <- seq(1,-1,by = -0.001)
  source("MSI_signature_discovery/MS_index_threshold_optimization.R")
  index_optim <- MS_index_thres_optim(MSindex=MSindex_train_sorted,
                                     threshold=threshold,
                                     MS_status=MS_status_train)

  ## test samples classification
  for(i in 1:ncol(data_sig_test)){
    for(j in 1:(length(possibilities))){
      correlations_test[j] <- cosine(data_sig_test[,i],centroid_mat[,j])
    }
    MSindex_test[i] <- correlations_test[1] - correlations_test[2] # GROUP OF INTEREST - NEGATIVE GROUP
    if(MSindex_test[i]  > index_optim$optim_threshold) {assignments_test[i] <- "MSI"} else {assignments_test[i] <- "MSS"}
  }
  names(assignments_test) <- names(MSindex_test)
  roc_curve <- roc(response = MS_status_test, predictor = MSindex_test, ci = TRUE)
  
  
  
  return(list(Centroids = centroid_mat,
              MSindex_optim = index_optim,
              Assignments_test = data.frame(MS_status_test=as.character(MS_status_test),class_result=assignments_test,MS_index=MSindex_test),
              # AUC_test = c(AUC_low, AUC, AUC_up)),
              roc_curve = roc_curve))
  
}

