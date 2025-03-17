run_dcats <- function(obj, Condition = "Condition") {
  require(DCATS)
  require(tidyverse)
  knn_mat = knn_simMat(obj@graphs$SCT_snn, obj$cell_type)
  
  count_mat = table(obj$orig.ident, obj$cell_type)
  
  design_mat = obj@meta.data |> 
    group_by(orig.ident, Condition) |> 
    summarize() |> 
    tibble::column_to_rownames("orig.ident")
  
  res <- dcats_GLM(count_mat, design_mat, similarity_mat = knn_mat)
  df <- as.data.frame(res)
  colnames(df) <- names(res)
  
  return(df)
}