applyGSVA = function(geneset, group_col, gene_col, expr_matList, kcdf = c("Gaussian", "Poisson")) {
  require(GSVA)
  require(BiocParallel)
  require(dplyr)
  
  stopifnot(inherits(geneset, "tbl_df") &
              inherits(group_col, "character") &
              inherits(gene_col, "character") &
              inherits(expr_matList, "list"))
  
  kcdf = match.arg(kcdf)
  
  for (i in seq_along(expr_matList)) {
    expr_mat = expr_matList[[i]]
    
    gsvapar = gsvaParam(exprData = expr_mat, geneSets = geneset, kcdf = kcdf)
    res = gsva(gsvapar, verbose = TRUE, BPPARAM = MulticoreParam(workers = 40))
  }
}