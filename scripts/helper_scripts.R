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

knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10(guide = "axis_logticks") + 
    scale_y_log10(guide = "axis_logticks") +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}