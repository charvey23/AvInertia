# This script contains all applicable plotting functions

plot_CGloc <- function(clean_pts,mass_properties,mass_properties_skin,mass_properties_bone,mass_properties_feathers,mass_properties_muscle){


  #--- Predefine the main theme ----
  th <- ggplot2::theme_classic() +
    ggplot2::theme(
      # Text
      axis.title  = ggplot2::element_text(size = 10),
      axis.text   = ggplot2::element_text(size = 10, colour = "black"),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, unit = "pt")),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      # Axis line
      axis.line   = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(-5,"pt"),
      # Legend
      legend.position  = 'none',
      # Background transparency
      # Background of panel
      panel.background = ggplot2::element_rect(fill = "transparent"),
      # Background behind actual data points
      plot.background  = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  CGplot <- ggplot2::ggplot() +
    # Add in data
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_skin, object == "CGy")$value, y=subset(mass_properties_skin, object == "CGx")$value), col = "#05B805") +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_bone, object == "CGy")$value, y=subset(mass_properties_bone, object == "CGx")$value), col = "#BB0000") +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_feathers, object == "CGy")$value, y=subset(mass_properties_feathers, object == "CGx")$value), col = "#249494") +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_muscle, object == "CGy")$value, y=subset(mass_properties_muscle, object == "CGx")$value), col = "#BB0000") +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties, component == "wing" & object == "CGy")$value, y=subset(mass_properties, component == "wing" & object == "CGx")$value), col = "black", size  = 3) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[,2], y=clean_pts[,1]), col = "gray") +
    # Theme
    th +
    # Axis control
    ggplot2::scale_y_continuous(name = "y (m)",limits = c(-0.4,0.01)) +
    ggplot2::scale_x_continuous(name = "x (m)", limits = c(0,0.6), position = "top") +
    ggthemes::geom_rangeframe() +
    ggplot2::coord_fixed() +
    ggplot2::annotate(geom = "segment", x = log(0), xend = log(0), y = -0.4, yend = 0) +
    ggplot2::annotate(geom = "segment", x = 0, xend = 0.6, y = -log(0), yend = -log(0))

  return(CGplot)
}
