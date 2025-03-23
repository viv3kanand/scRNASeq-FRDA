applyGSVA = function(geneset, group_col, gene_col, exprmatList, kcdf = c("Gaussian", "Poisson")) {
  require(GSVA)
  require(BiocParallel)
  require(dplyr)
  
  stopifnot(inherits(geneset, "tbl_df") &
              inherits(group_col, "character") &
              inherits(gene_col, "character") &
              inherits(exprmaListt, "matrix"))
  
  kcdf = match.arg(kcdf)
  
  for (i in seq_along(exprmatList)) {
    exprmat = exprmatList[[i]]
    
    gsvapar = gsvaParam(exprData = exprmat, geneSets = geneset, kcdf = kcdf)
    res = gsva(gsvapar, verbose = TRUE, BPPARAM = MulticoreParam(workers = 40))
  }
}