#######################################################################################################################################################################
## MSI_signature_validation.R: ##
#################################

## - validates the MSI signature using an independent  dataset (the nearest centroid classification according to the cosine correlation)

# Parameters: expr_matrix, MS_status, signature_discovery_output, platform
#       expr_matrix: normalized microarray/RNAseq expression matrix, features in rows, samples in columns
#       MS_status: named factor of MSI status (two levels: MSI/MSS); ordering of samples is the same as in microarray_expr
#       signature_discovery_output: list, output of signature_discovery.R function
#       platform: character, either "microarray" or "RNAseq", input expression data specification

#######################################################################################################################################################################


MSI_siganture_validation <- function(expr_matrix, MS_status, signature_discovery_output, platform)

{
  
  library(lsa)

  if(platform!="microarray" && platform!="RNAseq"){ message("platform must be 'microarray' of 'RNAseq' !!! ") }
  
  if(platform=="microarray")
  {
    
    # cosine correlation of each sample to centroids and assignment  
    data_sig <- expr_matrix[intersect(rownames(signature_discovery_output$centroids_microarray), rownames(expr_matrix)),]
    centroid_mat_final_model <- signature_discovery_output$centroids_microarray[intersect(rownames(signature_discovery_output$centroids_microarray),rownames(expr_matrix)),]
    correlations <- rep(NA, 2)
    assignments <- rep(NA, dim(data_sig)[2])
    MSindex <- rep(NA, dim(data_sig)[2])
    names(MSindex) <- colnames(expr_matrix)
    
    for(i in 1:(dim(data_sig)[2])){
      for(j in 1:2){
        correlations[j] <- cosine(data_sig[,i],centroid_mat_final_model[,j])
      }
      MSindex[i] <- correlations[1] - correlations[2] # GROUP OF INTEREST (MSI) - NEGATIVE GROUP (MSS)
      if(MSindex[i]  > signature_discovery_output$optimal_threshold_microarray) {assignments[i] <- "MSI"} else {assignments[i] <- "MSS"}
    }
    names(assignments) <- colnames(data_sig)
    
  } else {
    
    # cosine correlation of each sample to centroids and assignment  
    data_sig <- expr_matrix[intersect(rownames(signature_discovery_output$centroids_RNAseq), rownames(expr_matrix)),]
    centroid_mat_final_model <- signature_discovery_output$centroids_RNAseq[intersect(rownames(signature_discovery_output$centroids_RNAseq), rownames(expr_matrix)),]
    correlations <- rep(NA, 2)
    assignments <- rep(NA, dim(data_sig)[2])
    MSindex <- rep(NA, dim(data_sig)[2])
    names(MSindex) <- colnames(expr_matrix)
    
    for(i in 1:(dim(data_sig)[2])){
      for(j in 1:2){
        correlations[j] <- cosine(data_sig[,i],centroid_mat_final_model[,j])
      }
      MSindex[i] <- correlations[1] - correlations[2] # GROUP OF INTEREST (MSI) - NEGATIVE GROUP (MSS)
      if(MSindex[i]  > signature_discovery_output$optimal_threshold_RNAseq) {assignments[i] <- "MSI"} else {assignments[i] <- "MSS"}
    }
    names(assignments) <- colnames(data_sig)

  }
  
  roc_curve <- roc(response = MS_status, predictor = MSindex, ci = TRUE)
  
  
  list(Assignments = assignments,
       ROC_pROC=roc_curve)
  
}
