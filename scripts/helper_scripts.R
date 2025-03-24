# empirical filter functions
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

# zscore
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


# plot functions
knee_plot <- function(bc_rank) {
  require(DropletUtils)
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) |> 
    distinct() |> 
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

# https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p) {
  library(gtable)
  library(lemon)
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

font <- "Helvetica"
text_size <- 6
line_size <- 0.5

custom_theme <- function(){
  theme_classic() +
    theme(text = element_text(family = font,
                              size = text_size),
          line = element_line(linewidth = line_size),
          legend.title = element_text(size = rel(1.25)),
          axis.title = element_text(size = rel(1.25)),
          axis.text = element_text(size = rel(1.15)),
          plot.title = element_text(size = rel(1.5),
                                    hjust = 0.5),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          legend.position = "right",
          legend.text=element_text(size = rel(1.15)))
}

empty_axis <- function(){
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
}