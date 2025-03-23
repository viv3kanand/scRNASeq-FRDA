require(GSVA)
require(BiocParallel)
require(dplyr)


applyGSVA = function(geneset, group_col, gene_col, exprmat, kcdf = c("Gaussian", "Poisson")) {
  
  stopifnot(inherits(geneset, "tbl_df") &
              inherits(group_col, "character") &
              inherits(gene_col, "character") &
              inherits(exprmat, "matrix"))
  
  kcdf = match.arg(kcdf)
  
  gsvapar = gsvaParam(exprData = exprmat, geneSets = geneset, kcdf = kcdf)
  res = gsva(gsvapar, verbose = TRUE, BPPARAM = MulticoreParam(workers = 40))
}