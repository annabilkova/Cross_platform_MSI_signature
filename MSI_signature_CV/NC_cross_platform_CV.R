#######################################################################################################################################################################
## NC_cross_platform_CV.R: ##
#############################

## - computes the performance of the classifier using k-fold cross validation

# Parameters: microarray_expr, RNAseq_expr, MS_status_microarray, MS_status_RNAseq, lfc, adj_p_val_DEGs, K, batch_RNAseq, Tian/Lanza/Kruhoffer/Giacomini
#       microarray_expr: microarray expression matrix, features in rows, samples in columns, before quantile normalization!
#       RNAseq_expr: DGEList object with (at least) elements counts (table of unadjusted counts) and samples (data frame containing information about experimental group, library size and normalization factor for the library size), before normalization
#       MS_status_microarray: named factor of MSI status (two levels: MSI/MSS); ordering of samples is the same as in microarray_expr
#       MS_status_RNAseq: named factor of MSI status (two levels: MSI/MSS); ordering of samples is the same as in RNAseq_expr
#       lfc: LFC used in limma/edgeR to filter DEGs
#       adj_p_val_DEGs: adjusted p-value used in limma/edgeR to filter DEGs
#       K: integer, K-fold CV
#       batch_RNAseq: named character vector including batches of RNAseq dataset
#       Tian/Lanza/Kruhoffer/Giacomini - character vectors of given signatures ("entrezXXXX")

#######################################################################################################################################################################


NC_cross_platform_CV <- function(microarray_expr, RNAseq_expr, MS_status_microarray, MS_status_RNAseq, lfc, adj_p_val_DEGs, K, batch_RNAseq, Tian, Lanza, Kruhoffer, Giacomini)
{
  source("MSI_signature_CV/RNAseq_edgeR.R")
  source("MSI_signature_CV/nearest_centroid_cosine.R")
  split_text <- function(t,s,ch){unlist(lapply(strsplit(t,s,fixed = T),FUN=function(x) x[ch]))}
  
  library(preprocessCore)
  library(gplots)
  library(lsa)
  library(pROC)

  
  Assignments_rnaseq <- c()
  Assignments_microarray <- c()
  
  AUCMatrix_test <- matrix(NA,K,2,dimnames = list(1:K,c("microarray","RNAseq")))
  
  ## microarray 
  MSI_idx_microarray <- which(MS_status_microarray=="MSI")
  MSS_idx_microarray <- which(MS_status_microarray!="MSI")
  set.seed(5670)
  MSI_idx_random_samples_microarray <- sample(rep(1:K,length=length(MSI_idx_microarray)))
  set.seed(942)
  MSS_idx_random_samples_microarray <- sample(rep(1:K,length=length(MSS_idx_microarray)))
  classification_microarray <- c()
  
  ## rnaseq
  MSI_idx_RNAseq <- which(MS_status_RNAseq=="MSI")
  MSS_idx_RNAseq <- which(MS_status_RNAseq!="MSI")
  set.seed(2617)
  MSI_idx_random_samples_RNAseq <- sample(rep(1:K,length=length(MSI_idx_RNAseq)))
  set.seed(122)
  MSS_idx_random_samples_RNAseq <- sample(rep(1:K,length=length(MSS_idx_RNAseq)))
  classification_RNAseq <- c()
  

  # K-fold CV
  for(k in c(1:K))
  {
    
    print(paste("iteration",k))
    
    ## array train dataset
    trainIndex_microarray <- c(MSI_idx_microarray[MSI_idx_random_samples_microarray!=k],MSS_idx_microarray[MSS_idx_random_samples_microarray!=k])
    testIndex_microarray <- c(MSI_idx_microarray[MSI_idx_random_samples_microarray==k],MSS_idx_microarray[MSS_idx_random_samples_microarray==k])
    
    train_sample_ID_microarray <- colnames(microarray_expr)[trainIndex_microarray]
    test_sample_ID_microarray <- colnames(microarray_expr)[testIndex_microarray]
    train_MS_status_microarray <- MS_status_microarray[trainIndex_microarray]
    test_MS_status_microarray <- MS_status_microarray[testIndex_microarray]
    
    microarray_train <- microarray_expr[,train_sample_ID_microarray]
    print(paste("The ratio of MSI samples in train group: ",round(sum(train_MS_status_microarray=="MSI")/length(train_MS_status_microarray),3)))
    microarray_train_qn <- normalize.quantiles(microarray_train)
    colnames(microarray_train_qn) <- colnames(microarray_train)
    rownames(microarray_train_qn) <- rownames(microarray_train)
    
    microarray_test <- microarray_expr[,test_sample_ID_microarray]
    print(paste("The ratio of MSI samples in test group: ",round(sum(test_MS_status_microarray=="MSI")/length(test_MS_status_microarray),3)))
    microarray_test_qn <- normalize.quantiles(microarray_test)
    colnames(microarray_test_qn) <- colnames(microarray_test)
    rownames(microarray_test_qn) <- rownames(microarray_test)
    
    
    ## RNAseq train dataset (edgeR)
    trainIndex_RNAseq <- c(MSI_idx_RNAseq[MSI_idx_random_samples_RNAseq!=k],MSS_idx_RNAseq[MSS_idx_random_samples_RNAseq!=k])
    testIndex_RNAseq <- c(MSI_idx_RNAseq[MSI_idx_random_samples_RNAseq==k],MSS_idx_RNAseq[MSS_idx_random_samples_RNAseq==k])
    
    train_sample_ID_RNAseq <- colnames(RNAseq_expr)[trainIndex_RNAseq]
    train_MS_status_RNAseq <- MS_status_RNAseq[trainIndex_RNAseq]
    test_sample_ID_RNAseq <- colnames(RNAseq_expr)[testIndex_RNAseq]
    test_MS_status_RNAseq <- MS_status_RNAseq[testIndex_RNAseq]

    batch_train <- factor(batch_RNAseq[train_sample_ID_RNAseq])
    batch_test <- factor(batch_RNAseq[test_sample_ID_RNAseq])
    
    RNAseq_train <- RNAseq_expr[,train_sample_ID_RNAseq]
    print(paste("The ratio of MSI samples in train group: ",round(sum(train_MS_status_RNAseq=="MSI")/length(train_MS_status_RNAseq),3)))
    RNAseq_train_norm <- cpm(RNAseq_train, normalized.lib.sizes=TRUE, log = TRUE,prior.count = 2)
    RNAseq_train_norm <- removeBatchEffect(RNAseq_train_norm, batch_train)
    
    RNAseq_test <- RNAseq_expr[,test_sample_ID_RNAseq]
    print(paste("The ratio of MSI samples in test group: ",round(sum(test_MS_status_RNAseq=="MSI")/length(test_MS_status_RNAseq),3)))
    RNAseq_test_norm <- cpm(RNAseq_test, normalized.lib.sizes=TRUE,log = TRUE,prior.count = 2)
    RNAseq_test_norm <- removeBatchEffect(RNAseq_test_norm, batch_test)
    
    RNAseq_dataset_allsig <-RNAseq_edgeR(train_MS_status_RNAseq = train_MS_status_RNAseq, batch_train = batch_train, RNAseq_train = RNAseq_train, lfc = lfc, adj_p_val_DEGs = adj_p_val_DEGs)
    
    
    ### Signature
    input <- list(Tian = Tian,
                  Lanza = Lanza,
                  RNA_seq_A1 = rownames(RNAseq_dataset_allsig),
                  Kruhoffer = Kruhoffer,
                  Giacomini = Giacomini)
    tmp <- venn(input, show.plot=FALSE)
    int_all <- attr(tmp, "intersections")
    int_all <- int_all[regexpr(":",names(int_all))>0]
    EntrezIDs_signature <- unlist(int_all[regexpr("RNA_seq_A1",names(int_all))>0])
    EntrezIDs_signature <- intersect(EntrezIDs_signature,rownames(microarray_expr))
    print(paste("Number of genes in signature:",length(EntrezIDs_signature)))
    
    
    # multicollinearity testing in microarray train QN matrix
    int_mat <- abs(cor(t(microarray_train_qn[EntrezIDs_signature,])))
    int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
    
    while(sum(int_mat>0.75)>0)
    {
      print("removing multicollinearity in microarray")
      vtr <- as.numeric(names(sort(table(as.vector(which(int_mat>0.75,arr.ind = TRUE))),decreasing = TRUE)[1]))
      print(EntrezIDs_signature[vtr])
      EntrezIDs_signature <- EntrezIDs_signature[-vtr]
      int_mat <- abs(cor(t(microarray_train_qn[EntrezIDs_signature,])))
      int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
    }
    
    # multicollinearity testing in rnaseq train norm matrix    
    int_mat <- abs(cor(t(RNAseq_train_norm[EntrezIDs_signature,])))
    int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
    
    while(sum(int_mat>0.75)>0)
      {
        print("removing multicollinearity in rnaseq")
        vtr <- as.numeric(names(sort(table(as.vector(which(int_mat>0.75,arr.ind = TRUE))),decreasing = TRUE)[1]))
        print(EntrezIDs_signature[vtr])
        EntrezIDs_signature <- EntrezIDs_signature[-vtr]
        int_mat <- abs(cor(t(RNAseq_train_norm[EntrezIDs_signature,])))
        int_mat[lower.tri(int_mat,diag=TRUE)] <- 0
      }
    
    print(paste("Number of genes in signature:",length(EntrezIDs_signature)))
    
    nc_result_microarray <- nearest_centroid_cosine(data_train = microarray_train_qn, 
                                                    data_test = microarray_test_qn, 
                                                    signature = EntrezIDs_signature, 
                                                    MS_status_train = MS_status_microarray[trainIndex_microarray], 
                                                    MS_status_test = MS_status_microarray[testIndex_microarray], 
                                                    class_of_interest="MSI")
    
    AUCMatrix_test[k,which(colnames(AUCMatrix_test)=="microarray")] <- as.numeric(nc_result_microarray$roc_curve$auc)
    Assignments_microarray <- rbind(Assignments_microarray,nc_result_microarray$Assignments_test)
    
    nc_result_RNAseq <- nearest_centroid_cosine(data_train = RNAseq_train_norm, 
                                                data_test = RNAseq_test_norm, 
                                                signature = EntrezIDs_signature, 
                                                MS_status_train = MS_status_RNAseq[trainIndex_RNAseq], 
                                                MS_status_test = MS_status_RNAseq[testIndex_RNAseq], 
                                                class_of_interest="MSI")
    
    AUCMatrix_test[k,which(colnames(AUCMatrix_test)=="RNAseq")] <- as.numeric(nc_result_RNAseq$roc_curve$auc)
    Assignments_rnaseq <- rbind(Assignments_rnaseq,nc_result_RNAseq$Assignments_test)
      
  }
  
  
  AUC_res_microarray <- roc(response = Assignments_microarray$MS_status_test, predictor = Assignments_microarray$MS_index, ci = TRUE)
  AUC_res_RNAseq <- roc(response = Assignments_rnaseq$MS_status_test, predictor = Assignments_rnaseq$MS_index, ci = TRUE)
  
  
  return(list(p.value=adj_p_val_DEGs,
              AUC_microarray=AUC_res_microarray,
              AUC_RNAseq=AUC_res_RNAseq,
              lfc=lfc,
              KfoldCV=K,
              Assignments_rnaseq=Assignments_rnaseq,
              Assignments_microarray=Assignments_microarray))
  
}
