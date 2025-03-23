require(GSVA)
require(BiocParallel)
require(dplyr)


applyGSVA = function(geneset, group_col, gene_col, exprmat, kcdf = c("Gaussian", "Poisson")) {
  
  
  gsvapar = gsvaParam(exprData = exprmat, geneSets = geneset, kcdf = "Gaussian")
  res = gsva(gsvapar, verbose = TRUE, BPPARAM = MulticoreParam(workers = 40))
}