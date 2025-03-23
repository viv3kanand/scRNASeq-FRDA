applyGSVA = function(gset_df, group_col, gene_col, expr_mat_list, kcdf = c("Gaussian", "Poisson")) {
  require(GSVA)
  require(BiocParallel)
  require(dplyr)
  
  stopifnot(inherits(geneset, "tbl_df") &
              inherits(group_col, "character") &
              inherits(gene_col, "character") &
              inherits(expr_mat_list, "list"))
  
  kcdf = match.arg(kcdf)
  
  gset_list <- split(gset_df[gene_col], gset_df[group_col])
  
  for (i in seq_along(expr_mat_list)) {
    expr_mat = expr_mat_list[[i]]
    
    expr_mat = as.data.frame(expr_mat)
    rownames(expr_mat) = expr_mat[[1]]
    expr_mat = expr_mat[,-1, drop = FALSE] |> as.matrix()
    
    gsvapar = gsvaParam(exprData = expr_mat, geneSets = geneset, kcdf = kcdf)
    res = gsva(gsvapar, verbose = TRUE, BPPARAM = MulticoreParam(workers = 40))
    
    resList[[i]] = as.data.frame(t(res))
  }
  return(resList)
}