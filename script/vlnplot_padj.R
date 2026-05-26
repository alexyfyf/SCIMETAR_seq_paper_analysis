library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

# hacky way to add p adjust for ggpubr stat_compare_means when using comparison
# https://github.com/kassambara/ggpubr/issues/65

vlnplot_padj <- function(seuratobj, 
                         vars = c("dCt_adj", "SingleR.DevCell", "Treatment", "genotype_srsf2"),
                         facetorder = names(col_ref), 
                         test = "t.test", padj = 'fdr',
                         fill = col_group[1:2],
                         col = col_genotype,
                         lims = c(-1.5,12)) {
  
  # unpack your inputs
  yvar     <- vars[1]   # measurement on y‐axis
  facetvar <- vars[2]   # faceting variable
  xvar     <- vars[3]   # x‐axis / fill / colour
  colvar  <- vars[4]
  
  # 1. Fetch & prep data
  df <- FetchData(seuratobj, vars = vars) %>%
    # ensure grouping is a factor
    mutate(!!xvar := factor(.data[[xvar]]),
           !!facetvar := factor(.data[[facetvar]], levels = facetorder))
  
  # stat.test <- compare_means(
  #   formula         = as.formula(paste(yvar, "~", xvar)),
  #   group.by        = facetvar,
  #   data            = df,
  #   method          = test,
  #   p.adjust.method = "BH"
  # )
  
  # print(stat.test)
  
  # 2. Build a matching palette
  
  # 3. Build the ggplot with tidy‐eval
  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_violin(aes(fill = .data[[xvar]]),
                trim = FALSE, scale = "width", colour = NA) +
    # use facet_grid with your variable in the columns
    facet_grid(cols = vars(.data[[facetvar]]), scales = "free_y") +
    scale_fill_manual(values = fill) +
    scale_color_manual(values = col) +
    stat_compare_means(
      method      = test,
      label       = "p.format",
      hide.ns     = TRUE,
      tip.length  = 0.02
    ) +
    scale_y_continuous(limits = lims) +
    theme_classic() +
    theme(
      strip.text      = element_text(face = "bold"),
      axis.title.x    = element_blank(),
      # legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
  
  if (!is.na(colvar)) {
    p <- p + 
      geom_jitter(aes(col = .data[[colvar]]) ,
                  position = position_jitterdodge(
                    jitter.width = 0.1,    # controls vertical jitter
                    dodge.width = 0.9     # controls horizontal separation between groups
                  ),
                  size     = 0.8, 
                  alpha    = 0.8,
                  shape    = 16) 
  } else  {
    p <- p + 
      geom_jitter(position = position_jitter(width = 0.15),
                size     = 0.8, col = 'black',
                alpha    = 0.6,
                shape    = 16) 
  }
  
  # 4. Grab the raw p‐values from that stat layer and do a global FDR adjust
  pb      <- ggplot_build(p)
  layer_i <- which(vapply(pb$data, function(d) "p" %in% names(d), logical(1)))[1]
  pb$data[[layer_i]]$label <- 
    p.adjust(pb$data[[layer_i]]$p, method = padj) %>% signif(2)
  
  # 5. Draw the updated plot
  gt <- ggplot_gtable(pb)
  grid::grid.newpage()
  grid::grid.draw(gt)
  invisible(gt)
}

