#######################################################################################################################################################################
## signature_discovery.R: ##
############################

## - defines the final MSI gene expression signature from RNAseq experiment and microarray experiment
## - identifies the centroids for two groups: MSI and MSS

# Parameters: DEGs_RNAseq_expr, microarray_expr, MS_status_RNAseq, MS_status_array, Tian/Lanza/Kruhoffer/Giacomini
#    DEGs_RNAseq_expr - matrix of expression values; differentially expressed genes between phenotypes of interest in RNA-seq experiment; features in rows, samples in columns
#    microarray_expr - matrix of expression values; features in rows, samples in columns
#    MS_status_RNAseq - named factor of MSI status (two levels: MSI/MSS); ordering of samples is the same as in DEGs_RNAseq_expr
#    MS_status_microarray - named factor of MSI status (two levels: MSI/MSS); ordering of samples is the same as in microarray_expr
#    Tian/Lanza/Kruhoffer/Giacomini - character vectors of given signatures ("entrezXXXX")

#######################################################################################################################################################################


signature_discovery <- function(DEGs_RNAseq_expr, microarray_expr, MS_status_RNAseq, MS_status_microarray, Tian, Lanza, Kruhoffer, Giacomini) 
{
  
  library(lsa)
  library(gplots)

  source("MSI_signature_discovery/MS_index_threshold_optimization.R")
  
  input <- list(Tian = Tian,
                Lanza = Lanza,
                RNA_seq_A1 = rownames(DEGs_RNAseq_expr),
                Kruhoffer = Kruhoffer,
                Giacomini = Giacomini)
  
  tmp <- venn(input, show.plot=FALSE)
  int_all <- attr(tmp, "intersections")
  int_all <- int_all[regexpr(":",names(int_all))>0]
  signature <- unlist(int_all[regexpr("RNA_seq_A1",names(int_all))>0])
  print(paste("Missing genes in microarray:",setdiff(signature,rownames(microarray_expr))))
  signature <- intersect(signature,rownames(microarray_expr))
  
  print(paste("Number of genes in signature:",length(signature)))

  ### multicollinearity testing in microarray Quantile Normalised matrix
  int_mat <- abs(cor(t(microarray_expr[signature,])))
  int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
  
  while(sum(int_mat>0.75)>0)
  {
    print("removing multicollinearity in microarray")
    vtr <- as.numeric(names(sort(table(as.vector(which(int_mat>0.75,arr.ind = TRUE))),decreasing = TRUE)[1]))
    print(signature[vtr])
    signature <- signature[-vtr]
    int_mat <- abs(cor(t(microarray_expr[signature,])))
    int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
  }
  
  ### multicollinearity testing in RNAseq normalised matrix 
  int_mat <- abs(cor(t(DEGs_A1development[signature,])))
  int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
  
  while(sum(int_mat>0.75)>0)
  {
    print("removing multicollinearity in rnaseq")
    vtr <- as.numeric(names(sort(table(as.vector(which(int_mat>0.75,arr.ind = TRUE))),decreasing = TRUE)[1]))
    print(signature[vtr])
    signature <- signature[-vtr]
    int_mat <- abs(cor(t(DEGs_A1development[signature,])))
    int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
  }
  
  print(paste("Number of genes in signature:",length(signature)))
  
  
  ## microarray experiment
  signautre_expr <- microarray_expr[signature,]
  
  ### centroids of MSI/MSS
  centroids_microarray <- matrix(NA,nrow=length(rownames(signautre_expr)),ncol=2)
  rownames(centroids_microarray) <- signature
  
  
  if(is.factor(MS_status_microarray))
  {
    if(all(names(MS_status_microarray)==colnames(microarray_expr)))
    {
      colnames(centroids_microarray) <- levels(MS_status_microarray)
      
      for(i in c(1:length(levels(MS_status_microarray))))
      {
        sub_samples <- colnames(signautre_expr)[MS_status_microarray==levels(MS_status_microarray)[i]]
        centroids_microarray[,i] <- apply(signautre_expr[signature,sub_samples],1,mean)
      }
    } else {message("names(MS_status_microarray) differs from colnames(microarray_expr)")}
    
    
  } else {message("MS_status_microarray is NOT a factor")}
  

  ### cosine correlation of each sample to centroids and assignment
  correlations_microarray <- rep(NA, 2)
  assignments_microarray <- rep(NA, ncol(signautre_expr))
  MSindex_microarray <- rep(NA, ncol(signautre_expr))
  names(MSindex_microarray) <- colnames(signautre_expr)
  
  for(i in 1:ncol(signautre_expr)){
    for(j in 1:2){
      correlations_microarray[j] <- cosine(signautre_expr[,i],centroids_microarray[,j])
    }
    MSindex_microarray[i] <- correlations_microarray[1] - correlations_microarray[2] # GROUP OF INTEREST - NEGATIVE GROUP
  }
  
  MSindex_microarray_sorted <- sort(MSindex_microarray)
  optim_thres_microarray <- MS_index_thres_optim(MSindex = MSindex_microarray_sorted, threshold = seq(1,-1,by = -0.001), MS_status = MS_status_microarray[names(MSindex_microarray_sorted)])
  


  ## RNAseq experiment
  signautre_expr <- DEGs_RNAseq_expr[signature,]
  
  ### centroids of MSI/MSS
  centroids_RNAseq <- matrix(NA,nrow=length(rownames(signautre_expr)),ncol=2)
  rownames(centroids_RNAseq) <- signature
  
  if(is.factor(MS_status_microarray))
  {
    if(all(names(MS_status_RNAseq)==colnames(DEGs_RNAseq_expr)))
    {
      colnames(centroids_RNAseq) <- levels(MS_status_microarray)
      
      for(i in c(1:length(levels(MS_status_microarray))))
      {
        sub_samples <- colnames(signautre_expr)[MS_status_RNAseq==levels(MS_status_microarray)[i]]
        centroids_RNAseq[,i] <- apply(signautre_expr[signature,sub_samples],1,mean)
      }
    } else {message("names(MS_status_RNAseq) differs from colnames(DEGs_RNAseq_expr)")}
    
    
  } else {message("MS_status_RNAseq is NOT a factor")}
  
  
  ### cosine correlation of each sample to centroids and assignment
  correlations_RNAseq <- rep(NA, 2)
  assignments_RNAseq <- rep(NA, ncol(signautre_expr))
  MSindex_RNAseq <- rep(NA, ncol(signautre_expr))
  names(MSindex_RNAseq) <- colnames(signautre_expr)
  
  for(i in 1:ncol(signautre_expr)){
    for(j in 1:2){
      correlations_RNAseq[j] <- cosine(signautre_expr[,i],centroids_RNAseq[,j])
    }
    MSindex_RNAseq[i] <- correlations_RNAseq[1] - correlations_RNAseq[2] # GROUP OF INTEREST - NEGATIVE GROUP
  }
  
  MSindex_RNAseq_sorted <- sort(MSindex_RNAseq)
  optim_thres_RNAseq <- MS_index_thres_optim(MSindex = MSindex_RNAseq_sorted, threshold = seq(1,-1,by = -0.001), MS_status = MS_status_RNAseq[names(MSindex_RNAseq_sorted)])
  
  
  return (list(signature=signature, 
               centroids_microarray=centroids_microarray, 
               centroids_RNAseq=centroids_RNAseq, 
               optimal_threshold_microarray=optim_thres_microarray$optim_threshold,
               optimal_threshold_RNAseq=optim_thres_RNAseq$optim_threshold))
  
}


