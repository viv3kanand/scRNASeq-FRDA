emp_fun <- function(obj){
  obj@meta.data[, c("nCount_RNA", "nFeature_RNA")] |> 
    map(log10) |> 
    map(~c(10^(mean(.x) + 3*sd(.x)), 10^(mean(.x) - 3*sd(.x))))
}

filt_fun <- function(obj, emp){
  PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt") |> 
    subset(subset = nFeature_RNA > 200 &
             nFeature_RNA < emp[[2]][1] & 
             nCount_RNA < emp[[1]][1] &
             percent.mt < 10)
}

filterCell <- function(obj){
  emp <- emp_fun(obj)
  combined <- filt_fun(obj, emp)
  
  combined$log10GenesPerUMI <- log10(combined$nFeature_RNA) / log10(combined$nCount_RNA)
  
  return(combined)
}