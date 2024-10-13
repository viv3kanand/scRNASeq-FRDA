filterCell <- function(obj){
  meta <- obj@meta.data[, c("nCount_RNA", "nFeature_RNA")] |> 
    map(log10) |> 
    map(~c(10^(mean(.x) + 3*sd(.x)), 10^(mean(.x) - 3*sd(.x))))
  
  combined <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt") |> 
    subset(subset = nFeature_RNA > 200 & 
             nFeature_RNA < meta[[2]][1] & 
             nCount_RNA < meta[[1]][1] &
             percent.mt < 10)
  
  combined$log10GenesPerUMI <- log10(combined$nFeature_RNA) / log10(combined$nCount_RNA)
  
  return(combined)
}
