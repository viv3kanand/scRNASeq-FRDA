applyGSVA = function(gset_df, group_col, gene_col, expr_mat, kcdf = c("Gaussian", "Poisson")) {
  require(GSVA)
  require(BiocParallel)
  require(dplyr)

  stopifnot(inherits(gset_df, "tbl_df") &
              inherits(group_col, "character") &
              inherits(gene_col, "character") &
              inherits(expr_mat_list, "tbl_df"))
  
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required but not installed!")
  }
  
  kcdf = match.arg(kcdf)
  
  gset_list <- split(as.character(gset_df[[gene_col]]), as.character(gset_df[[group_col]]))
  
  expr_mat = vst_gsva |> 
    tibble::column_to_rownames("symbol") |>
    as.matrix()

  gsvapar = gsvaParam(exprData = expr_mat, geneSets = gset_list, kcdf = kcdf)
  res = gsva(gsvapar, verbose = TRUE, BPPARAM = MulticoreParam(workers = 40))
  
  return(res)
}