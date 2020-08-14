# This script contains all applicable plotting functions

plot_CGloc <- function(clean_pts,mass_properties,mass_properties_skin,mass_properties_bones,mass_properties_feathers,mass_properties_muscles){


  #--- Predfine the main theme ----
  th <- theme_classic() +
    theme(
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

  label_elb = "Elbow angle (?)"
  label_man = "Wrist angle (?)"

  CGplot <- ggplot() +
    #add in data
    geom_point(aes(x=subset(mass_properties_skin, object == "CGy")$value, y=subset(mass_properties_skin, object == "CGz")$value), col = "gray") +
    geom_point(aes(x=subset(mass_properties_bone, object == "CGy")$value, y=subset(mass_properties_bone, object == "CGz")$value), col = "yellow") +
    geom_point(aes(x=subset(mass_properties_feathers, object == "CGy")$value, y=subset(mass_properties_feathers, object == "CGz")$value), col = "blue") +
    geom_point(aes(x=subset(mass_properties_muscle, object == "CGy")$value, y=subset(mass_properties_muscle, object == "CGz")$value), col = "red") +
    geom_point(aes(x=subset(mass_properties, object == "CGy")$value, y=subset(mass_properties, object == "CGz")$value), col = "green") +
    geom_point(aes(x=clean_pts[,2], y=clean_pts[,3]), col = "black")
    #theme
    #th +
    #axis control
    #scale_x_continuous(name = label_elb,limits = c(30,170), breaks = c(30,50,70,90,110,130,150,170)) +
    #scale_y_continuous(name = label_man, limits = c(100,180), breaks = c(100,120,140,160,180)) +
    #geom_rangeframe() +
    #annotate(geom = "segment", x = 0, xend = 0, y = 100, yend = 180) +
    #annotate(geom = "segment", x = 30, xend = 170, y = 0, yend = 0)

  return(CGplot)
}
