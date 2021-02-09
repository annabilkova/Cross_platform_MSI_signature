RNAseq_edgeR <- function(train_MS_status_RNAseq, batch_train, RNAseq_train, lfc, adj_p_val_DEGs)
{
  train_MS_status_RNAseq <- relevel(train_MS_status_RNAseq,ref="MSI")
  
  design <- model.matrix(~ 0 + train_MS_status_RNAseq + batch_train)
  yc_disp <- estimateDisp(RNAseq_train,design)
  fit <- glmQLFit(yc_disp,design)
  contrast.matrix <- makeContrasts(MSIvsMSS=train_MS_status_RNAseqMSI-train_MS_status_RNAseqMSS, levels=design)
  tr <- glmTreat(fit, contrast=contrast.matrix, lfc=lfc)
  RNAseq_dataset_allsig <- topTags(tr, n=Inf, adjust.method="BH", p.value=adj_p_val_DEGs)
  
  return(RNAseq_dataset_allsig)
}