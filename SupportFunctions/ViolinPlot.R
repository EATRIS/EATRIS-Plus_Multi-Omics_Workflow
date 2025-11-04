plot_violin_unified <- function(data, x_var, y_var, facet_var = NULL, palette="jco") {
  
  plt <- ggviolin(data,
                  x = x_var,
                  y = y_var,
                  color = x_var,
                  add = c("jitter", "boxplot"),
                  palette = palette,
                  draw_quantiles = 0.5)
  
  # Add statistics (ANOVA by default)
  plt <- plt + stat_compare_means(method = "anova")
  
  # Facet if requested
  if (!is.null(facet_var)) {
    plt <- plt + facet_wrap(as.formula(paste("~", facet_var)))
  }
  
  # Unified theme
  plt <- plt + theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    text = element_text(size = 10)
  )
  
  return(plt)
}