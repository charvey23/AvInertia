# This script contains all applicable plotting functions

plot_CGloc <- function(clean_pts,mass_properties,mass_properties_skin,mass_properties_bone,mass_properties_feathers,mass_properties_muscle, prop_tertiary1, prop_tertiary2, plot_var){
  primaries   = mass_properties_feathers[grep("P",mass_properties_feathers$component),]
  secondaries = mass_properties_feathers[grep("S",mass_properties_feathers$component),]

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

  CGplot_x <- ggplot2::ggplot() +
    # Add in data
    ggplot2::geom_point(ggplot2::aes(x=prop_tertiary1$CG[2], y=prop_tertiary1$CG[1]), col = "black",fill = "#BEA4BD", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=prop_tertiary2$CG[2], y=prop_tertiary2$CG[1]), col = "black",fill = "#BEA4BD", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_skin, object == "CGy")$value, y=subset(mass_properties_skin, object == "CGx")$value), col = "black", fill = "#BEA4BD", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_bone, object == "CGy")$value, y=subset(mass_properties_bone, object == "CGx")$value), col = "black", fill = "#FAC5C6", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(primaries, object == "CGy")$value, y=subset(primaries, object == "CGx")$value), col = "black", fill = "#A0B3DC", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(secondaries, object == "CGy")$value, y=subset(secondaries, object == "CGx")$value), col = "black", fill = "#9AD09B", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_muscle, object == "CGy")$value, y=subset(mass_properties_muscle, object == "CGx")$value), col = "black", fill = "#FAC5C6", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties, component == "wing" & object == "CGy")$value, y=subset(mass_properties, component == "wing" & object == "CGx")$value), col = "black", size  = 3) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[1:4,2], y=clean_pts[1:4,1]), col = "black", fill = "white", pch = 22) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[5:10,2], y=clean_pts[5:10,1]), col = "gray50", pch = 15) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[5:10,2], y=clean_pts[5:10,1]), col = "black", pch = 0) +
    ggplot2::geom_point(ggplot2::aes(x=0.11005137, y=0.00354631), col = "gray50", pch = 15) +
    ggplot2::geom_point(ggplot2::aes(x=0.11005137, y=0.00354631), col = "black", pch = 0) +
    # Theme
    th +
    # Axis control
    ggplot2::scale_y_continuous(name = "x (m)",limits = c(-0.4,0.01)) +
    ggplot2::scale_x_continuous(name = "y (m)", limits = c(0,0.6), position = "top") +
    ggthemes::geom_rangeframe() +
    ggplot2::coord_fixed() +
    ggplot2::annotate(geom = "segment", x = log(0), xend = log(0), y = -0.4, yend = 0) +
    ggplot2::annotate(geom = "segment", x = 0, xend = 0.6, y = -log(0), yend = -log(0))

  if(plot_var == "yx"){
    plot(CGplot_x)
  }

  CGplot_z <- ggplot2::ggplot() +
    # Add in data
    ggplot2::geom_point(ggplot2::aes(x=prop_tertiary1$CG[2], y=-prop_tertiary1$CG[3]), col = "black",fill = "#BEA4BD", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=prop_tertiary2$CG[2], y=-prop_tertiary2$CG[3]), col = "black",fill = "#BEA4BD", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_skin, object == "CGy")$value, y=-subset(mass_properties_skin, object == "CGz")$value), col = "black",fill = "#BEA4BD", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_bone, object == "CGy")$value, y=-subset(mass_properties_bone, object == "CGz")$value),col = "black", fill = "#FAC5C6", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(primaries, object == "CGy")$value, y=-subset(primaries, object == "CGz")$value), col = "black",fill = "#A0B3DC", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(secondaries, object == "CGy")$value, y=-subset(secondaries, object == "CGz")$value), col = "black",fill = "#9AD09B", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties_muscle, object == "CGy")$value, y=-subset(mass_properties_muscle, object == "CGz")$value), fill = "black",col = "#FAC5C6", pch = 21) +
    ggplot2::geom_point(ggplot2::aes(x=subset(mass_properties, component == "wing" & object == "CGy")$value, y=-subset(mass_properties, component == "wing" & object == "CGz")$value), col = "black", size  = 3) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[1:4,2], y=-clean_pts[1:4,3]), col = "black", fill = "white", pch = 22) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[5:10,2], y=-clean_pts[5:10,3]), col = "gray50", pch = 15) +
    ggplot2::geom_point(ggplot2::aes(x=clean_pts[5:10,2], y=-clean_pts[5:10,3]), col = "black", pch = 0) +
    ggplot2::geom_point(ggplot2::aes(x=0.11005137, y=-0.01009165), col = "gray50", pch = 15) +
    ggplot2::geom_point(ggplot2::aes(x=0.11005137, y=-0.01009165), col = "black", pch = 0) +
    # Theme
    th +
    # Axis control
    ggplot2::scale_y_continuous(name = "z (m)",limits = c(-0.2,0.2)) +
    ggplot2::scale_x_continuous(name = "y (m)", limits = c(0,0.6), breaks= c(0,0.1), position = "top") +
    ggthemes::geom_rangeframe() +
    ggplot2::coord_fixed() +
    ggplot2::annotate(geom = "segment", x = log(0), xend = log(0), y = -0.2, yend = 0.2) +
    ggplot2::annotate(geom = "segment", x = 0, xend = 0.6, y = -log(0), yend = -log(0))

  if(plot_var == "yz"){
    plot(CGplot_z)
    }

  return(CGplot_x)
}
