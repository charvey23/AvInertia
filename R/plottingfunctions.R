# This script contains all applicable plotting functions

plot_CGloc <- function(clean_pts,mass_properties,mass_properties_skin,mass_properties_bone,mass_properties_feathers,mass_properties_muscle){


  #--- Predfine the main theme ----
  th <- ggplot2::theme_classic() +
    ggplot2::theme(
      # Text
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10, colour = "black"),
      axis.text.x = element_text(margin = margin(t = 10, unit = "pt")),
      axis.text.y = element_text(margin = margin(r = 10)),
      # Axis line
      axis.line = element_blank(),
      axis.ticks.length = unit(-5,"pt"),
      # Legend
      legend.position = 'none',
      # Background transparency
      # Background of panel
      panel.background = element_rect(fill = "transparent"),
      # Background behind actual data points
      plot.background = element_rect(fill = "transparent", color = NA)
    )

  CGplot <- ggplot2::ggplot() +
    #add in data
    geom_point(aes(x=subset(mass_properties_skin, object == "CGy")$value, y=subset(mass_properties_skin, object == "CGx")$value), col = "#05B805") +
    geom_point(aes(x=subset(mass_properties_bone, object == "CGy")$value, y=subset(mass_properties_bone, object == "CGx")$value), col = "#BB0000") +
    geom_point(aes(x=subset(mass_properties_feathers, object == "CGy")$value, y=subset(mass_properties_feathers, object == "CGx")$value), col = "#249494") +
    geom_point(aes(x=subset(mass_properties_muscle, object == "CGy")$value, y=subset(mass_properties_muscle, object == "CGx")$value), col = "#BB0000") +
    geom_point(aes(x=subset(mass_properties, component == "wing" & object == "CGy")$value, y=subset(mass_properties, component == "wing" & object == "CGx")$value), col = "black", size  = 3) +
    geom_point(aes(x=clean_pts[,2], y=clean_pts[,1]), col = "gray") +
    #theme
    th +
    #axis control
    scale_y_continuous(name = "y (m)",limits = c(-0.4,0.01)) +
    scale_x_continuous(name = "x (m)", limits = c(0,0.6), position = "top") +
    geom_rangeframe() +
    coord_fixed() +
    annotate(geom = "segment", x = log(0), xend = log(0), y = -0.4, yend = 0) +
    annotate(geom = "segment", x = 0, xend = 0.6, y = -log(0), yend = -log(0))

  return(CGplot)
}
