test_max = aggregate(list(prop_q_dot = dat_final$prop_q_dot, full_m = dat_final$full_m),  by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), max)
test_min = aggregate(list(prop_q_dot = dat_final$prop_q_dot, full_m = dat_final$full_m),  by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), min)

ROM_plot <- ggplot()+
  geom_polygon(data = chulls_xz, aes(x = elbow, y = manus), col = "black", fill = NA, alpha = 0.5) +
  geom_point(data = dat_final[which(dat_final$species %in% test_max$species & dat_final$prop_q_dot %in% test_max$prop_q_dot),],
             aes(x = elbow, y = manus, col = log(prop_q_dot)), alpha = 0.9) +
  geom_point(data = dat_final[which(dat_final$species %in% test_min$species & dat_final$prop_q_dot %in% test_min$prop_q_dot),],
             aes(x = elbow, y = manus, col = log(prop_q_dot)), alpha = 0.9, pch = 17) +
  facet_wrap(~species_order, nrow = 2) +
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.line =  element_line(colour = "black", linetype=1))+
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)') +
  scale_y_continuous(name='Wrist angle (°)')

## to compute the body extremes

extremes = aggregate(list(BeakTipx_orgShoulder = dat_final$BeakTipx_orgShoulder,
                          Centrez_orgShoulder  = dat_final$Centrez_orgShoulder,
                          TailTipx_orgShoulder = dat_final$TailTipx_orgShoulder,
                          MaxWidthx_orgShoulder = dat_final$MaxWidthx_orgShoulder,
                          Dorsalz_orgShoulder  = dat_final$Dorsalz_orgShoulder,
                          Ventralz_orgShoulder = dat_final$Ventralz_orgShoulder),  by=list(species_order = dat_final$species_order), mean)
x_shape1 = melt(extremes[,c("species_order","BeakTipx_orgShoulder","MaxWidthx_orgShoulder","TailTipx_orgShoulder")], id= "species_order")
x_shape2 = melt(extremes[,c("species_order","MaxWidthx_orgShoulder","BeakTipx_orgShoulder")], id= "species_order")
x_shape = rbind(x_shape1,x_shape2)
z_shape1 = melt(extremes[,c("species_order","Centrez_orgShoulder","Dorsalz_orgShoulder")], id= "species_order")
z_shape2 = melt(extremes[,c("species_order","Centrez_orgShoulder","Ventralz_orgShoulder")], id= "species_order")
z_shape3 = melt(extremes[,c("species_order","Centrez_orgShoulder")], id= "species_order")
z_shape = rbind(z_shape1,z_shape2,z_shape3)

### Ixx vs Izz plots

chulls_Ixxzz <- ddply(dat_final[,c("species_order","BirdID","full_Ixx_specific","full_Izz_specific")], .(species_order, BirdID),
                      function(df) df[chull(df$full_Ixx_specific, df$full_Izz_specific), ])


fullIxxIzz_plot <- ggplot()+
  geom_polygon(data = chulls_Ixxzz,aes(x = full_Ixx_specific, y = full_Izz_specific, col = species_order, group = interaction(species_order,BirdID)), fill = NA, alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0) +
  #facet_wrap(~species_order, nrow = 3) +
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  th +
  coord_fixed()

## other I parameters

for (i in 1:nrow(dat_bird)){
  tmp            = subset(dat_final[,c("species","species_order","BirdID","elbow","manus",
                                       "full_CGx","wing_CGy","full_CGz",
                                       "full_Ixx","full_Iyy","full_Izz","full_Ixz",
                                       "pitch_div", "yaw_div")], species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])

  tmp$scaled_CGx = (tmp$full_CGx - min(tmp$full_CGx))/(max(tmp$full_CGx) - min(tmp$full_CGx))
  tmp$scaled_CGy = (tmp$wing_CGy - min(tmp$wing_CGy))/(max(tmp$wing_CGy) - min(tmp$wing_CGy))
  tmp$scaled_CGz = (tmp$full_CGz - min(tmp$full_CGz))/(max(tmp$full_CGz) - min(tmp$full_CGz))

  tmp$scaled_Ixx = (tmp$full_Ixx - min(tmp$full_Ixx))/(max(tmp$full_Ixx) - min(tmp$full_Ixx))
  tmp$scaled_Iyy = (tmp$full_Iyy - min(tmp$full_Iyy))/(max(tmp$full_Iyy) - min(tmp$full_Iyy))
  tmp$scaled_Izz = (tmp$full_Izz - min(tmp$full_Izz))/(max(tmp$full_Izz) - min(tmp$full_Izz))
  tmp$scaled_Ixz = (tmp$full_Ixz - min(tmp$full_Ixz))/(max(tmp$full_Ixz) - min(tmp$full_Ixz))

  tmp$scaled_pitch_div = (tmp$pitch_div - min(tmp$pitch_div))/(max(tmp$pitch_div) - min(tmp$pitch_div))
  tmp$scaled_yaw_div   = (tmp$yaw_div - min(tmp$yaw_div))/(max(tmp$yaw_div) - min(tmp$yaw_div))
  if (i == 1){
    dat_scaled = tmp
  }else{
    dat_scaled = rbind(dat_scaled,tmp)
  }
}

ROM_plot <- ggplot()+
  geom_point(data = dat_scaled, aes(x = elbow, y = manus, col = scaled_pitch_div)) +
  facet_wrap(~species_order) + th


ROM_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = elbow, y = manus, col = yaw_div)) +
  facet_wrap(~species_order) + th





## ----------------- Plot of linear model outputs ------------------
# order species to match the phylogeny
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
col_xz = "#af0069"

Ixx_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars

  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Ixx_sp_man_lb1,min_Ixx_sp_man_lb2), xmax = pmax(max_Ixx_sp_man_ub1,max_Ixx_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Ixx_sp_elb_lb1,min_Ixx_sp_elb_lb2), xmax = pmax(max_Ixx_sp_elb_ub1,max_Ixx_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_man+mean_Ixx_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_man+mean_Ixx_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_man+mean_Ixx_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_elb+mean_Ixx_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_elb+mean_Ixx_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_elb+mean_Ixx_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-0.05,0.15), breaks = c(-0.05,0,0.05,0.1,0.15), labels = c(-0.005,0,0.005,0.01,0.015)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.05, xend = 0.15, y = log(0), yend = log(0))


Iyy_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars

  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Iyy_sp_man_lb1,min_Iyy_sp_man_lb2), xmax = pmax(max_Iyy_sp_man_ub1,max_Iyy_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Iyy_sp_elb_lb1,min_Iyy_sp_elb_lb2), xmax = pmax(max_Iyy_sp_elb_ub1,max_Iyy_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_man+mean_Iyy_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_man+mean_Iyy_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_man+mean_Iyy_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_elb+mean_Iyy_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_elb+mean_Iyy_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_elb+mean_Iyy_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-0.05,0.15), breaks = c(-0.05,0,0.05,0.1,0.15), labels = c(-0.005,0,0.005,0.01,0.015)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.05, xend = 0.15, y = log(0), yend = log(0))



Izz_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Izz_sp_man_lb1,min_Izz_sp_man_lb2), xmax = pmax(max_Izz_sp_man_ub1,max_Izz_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Izz_sp_elb_lb1,min_Izz_sp_elb_lb2), xmax = pmax(max_Izz_sp_elb_ub1,max_Izz_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_man+mean_Izz_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_man+mean_Izz_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_man+mean_Izz_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_elb+mean_Izz_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_elb+mean_Izz_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_elb+mean_Izz_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  # need to adjust
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-0.05,0.15), breaks = c(-0.05,0,0.05,0.1,0.15), labels = c(-0.005,0,0.005,0.01,0.015)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.05, xend = 0.15, y = log(0), yend = log(0))


# -------- Combine panels into figure ------------
#exported as 6x14
bottomrow <- plot_grid(phylo_plot,Ixx_specific,Iyy_specific, Izz_specific,Ixz_specific,
                       #arrangement data
                       ncol = 5,
                       rel_widths = c(1.5,1,1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")

bottomrow <- plot_grid(phylo_plot,Ixx_angleeffect,Iyy_angleeffect, Izz_angleeffect,
                       #arrangement data
                       ncol = 4,
                       rel_widths = c(1.5,1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")

#exported as 7.5x9
figure_final <- plot_grid(toprow,bottomrow,
                          #arrangement data
                          ncol = 1, nrow = 2, rel_heights = c(1,1))
# Visualize the wings as required - For each calibration verify that the axis is RH
m = 1:nrow(dat_wing_curr)
max = 0.4
plot(dat_wing$pt2_Y[m],dat_wing$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_wing$pt3_Y[m],dat_wing$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_wing$pt4_Y[m], dat_wing$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_wing$pt1_Y[m],dat_wing$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_wing$pt8_Y[m],dat_wing$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_wing$pt10_Y[m],dat_wing$pt10_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

plot(dat_clean$pt2_Y[m],dat_clean$pt2_Z[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_clean$pt3_Y[m],dat_clean$pt3_Z[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_clean$pt4_Y[m], dat_clean$pt4_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_clean$pt8_Y[m],dat_clean$pt8_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_clean$pt11_Y[m],dat_clean$pt11_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")
points(dat_clean$pt1_Y[m],dat_clean$pt1_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")

plot(dat_clean$pt2_Z[m],dat_clean$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_clean$pt3_Z[m],dat_clean$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_clean$pt4_Z[m], dat_clean$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_clean$pt1_Z[m],dat_clean$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_clean$pt8_Z[m],dat_clean$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_clean$pt10_Z[m],dat_clean$pt10_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

x = c(dat_clean$pt1_X[m],dat_clean$pt2_X[m],dat_clean$pt3_X[m],dat_clean$pt4_X[m],
      dat_clean$pt6_X[m],dat_clean$pt7_X[m],dat_clean$pt8_X[m],dat_clean$pt9_X[m],
      dat_clean$pt10_X[m],dat_clean$pt11_X[m],dat_clean$pt12_X[m])
y = c(dat_clean$pt1_Y[m],dat_clean$pt2_Y[m],dat_clean$pt3_Y[m],dat_clean$pt4_Y[m],
      dat_clean$pt6_Y[m],dat_clean$pt7_Y[m],dat_clean$pt8_Y[m],dat_clean$pt9_Y[m],
      dat_clean$pt10_Y[m],dat_clean$pt11_Y[m],dat_clean$pt12_Y[m])
z = c(dat_clean$pt1_Z[m],dat_clean$pt2_Z[m],dat_clean$pt3_Z[m],dat_clean$pt4_Z[m],
      dat_clean$pt6_Z[m],dat_clean$pt7_Z[m],dat_clean$pt8_Z[m],dat_clean$pt9_Z[m],
      dat_clean$pt10_Z[m],dat_clean$pt11_Z[m],dat_clean$pt12_Z[m])
plot3d(x, y, z, col = c("black","gray20","gray50","gray70","red", "orange", "yellow", "green", "blue", "navy", "purple"))

plot(dat_clean$elbow,dat_clean$manus)

plot(test$pt2_Y[m],test$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(test$pt3_Y[m],test$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(test$pt4_Y[m], test$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(test$pt1_Y[m],test$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(test$pt8_Y[m],test$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(test$pt9_Y[m],test$pt9_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

plot(test$pt2_Y[m],test$pt2_Z[m], xlim =c(-max,max), ylim = c(-max,max))
points(test$pt3_Y[m],test$pt3_Z[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(test$pt4_Y[m], test$pt4_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(test$pt1_Y[m],test$pt1_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(test$pt8_Y[m],test$pt8_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(test$pt9_Y[m],test$pt9_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")


points(tmp$pt2_Y[m],tmp$pt2_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(tmp$pt3_Y[m],tmp$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(tmp$pt4_Y[m], tmp$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(tmp$pt1_Y[m],tmp$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(tmp$pt8_Y[m],tmp$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(tmp$pt9_Y[m],tmp$pt9_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")


# ## TO CHECK DATA
# x_vertices = c(dat_wing$pt1_X[i], dat_wing$pt6_X[i],dat_wing$pt7_X[i],dat_wing$pt8_X[i],dat_wing$pt9_X[i],dat_wing$pt10_X[i],dat_wing$pt11_X[i],dat_wing$pt12_X[i])
# y_vertices = c(dat_wing$pt1_Y[i], dat_wing$pt6_Y[i],dat_wing$pt7_Y[i],dat_wing$pt8_Y[i],dat_wing$pt9_Y[i],dat_wing$pt10_Y[i],dat_wing$pt11_Y[i],dat_wing$pt12_Y[i])
# z_vertices = c(dat_wing$pt1_Z[i], dat_wing$pt6_Z[i],dat_wing$pt7_Z[i],dat_wing$pt8_Z[i],dat_wing$pt9_Z[i],dat_wing$pt10_Z[i],dat_wing$pt11_Z[i],dat_wing$pt12_Z[i])
#
#
x = c(dat_wing$pt1_X[m],dat_wing$pt2_X[m],dat_wing$pt3_X[m],dat_wing$pt4_X[m],
      dat_wing$pt6_X[m],dat_wing$pt7_X[m],dat_wing$pt8_X[m],dat_wing$pt9_X[m],
      dat_wing$pt10_X[m],dat_wing$pt11_X[m],dat_wing$pt12_X[m])
y = c(dat_wing$pt1_Y[m],dat_wing$pt2_Y[m],dat_wing$pt3_Y[m],dat_wing$pt4_Y[m],
      dat_wing$pt6_Y[m],dat_wing$pt7_Y[m],dat_wing$pt8_Y[m],dat_wing$pt9_Y[m],
      dat_wing$pt10_Y[m],dat_wing$pt11_Y[m],dat_wing$pt12_Y[m])
z = c(dat_wing$pt1_Z[m],dat_wing$pt2_Z[m],dat_wing$pt3_Z[m],dat_wing$pt4_Z[m],
      dat_wing$pt6_Z[m],dat_wing$pt7_Z[m],dat_wing$pt8_Z[m],dat_wing$pt9_Z[m],
      dat_wing$pt10_Z[m],dat_wing$pt11_Z[m],dat_wing$pt12_Z[m])
# plot3d(x, y, z, col = c("black","gray20","gray50","gray70","red", "orange", "yellow", "green", "blue", "navy", "purple"))


min_lb     <- aggregate(list(min_CGx_sp_int_lb = dat_model_out$CGx_sp_int_lb,
                             min_CGy_sp_int_lb = dat_model_out$CGy_sp_int_lb,
                             min_CGz_sp_int_lb = dat_model_out$CGz_sp_int_lb,
                             min_CGx_sp_elb_lb1 = dat_model_out$CGx_sp_elb_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[1],
                             min_CGx_sp_elb_lb2 = dat_model_out$CGx_sp_elb_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[3],
                             min_CGy_sp_elb_lb1 = dat_model_out$CGy_sp_elb_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[1],
                             min_CGy_sp_elb_lb2 = dat_model_out$CGy_sp_elb_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[3],
                             min_CGz_sp_elb_lb1 = dat_model_out$CGz_sp_elb_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[1],
                             min_CGz_sp_elb_lb2 = dat_model_out$CGz_sp_elb_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[3],
                             min_CGx_sp_man_lb1 = dat_model_out$CGx_sp_man_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[1],
                             min_CGx_sp_man_lb2 = dat_model_out$CGx_sp_man_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[3],
                             min_CGy_sp_man_lb1 = dat_model_out$CGy_sp_man_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[1],
                             min_CGy_sp_man_lb2 = dat_model_out$CGy_sp_man_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[3],
                             min_CGz_sp_man_lb1 = dat_model_out$CGz_sp_man_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[1],
                             min_CGz_sp_man_lb2 = dat_model_out$CGz_sp_man_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[3]),
                        by=list(species = dat_model_out$species), min)
max_ub     <- aggregate(list(max_CGx_sp_int_ub = dat_model_out$CGx_sp_int_ub,
                             max_CGy_sp_int_ub = dat_model_out$CGy_sp_int_ub,
                             max_CGz_sp_int_ub = dat_model_out$CGz_sp_int_ub,
                             max_CGx_sp_elb_ub1 = dat_model_out$CGx_sp_elb_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[1],
                             max_CGx_sp_elb_ub2 = dat_model_out$CGx_sp_elb_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[3],
                             max_CGy_sp_elb_ub1 = dat_model_out$CGy_sp_elb_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[1],
                             max_CGy_sp_elb_ub2 = dat_model_out$CGy_sp_elb_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[3],
                             max_CGz_sp_elb_ub1 = dat_model_out$CGz_sp_elb_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[1],
                             max_CGz_sp_elb_ub2 = dat_model_out$CGz_sp_elb_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[3],
                             max_CGx_sp_man_ub1 = dat_model_out$CGx_sp_man_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[1],
                             max_CGx_sp_man_ub2 = dat_model_out$CGx_sp_man_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[3],
                             max_CGy_sp_man_ub1 = dat_model_out$CGy_sp_man_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[1],
                             max_CGy_sp_man_ub2 = dat_model_out$CGy_sp_man_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[3],
                             max_CGz_sp_man_ub1 = dat_model_out$CGz_sp_man_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[1],
                             max_CGz_sp_man_ub2 = dat_model_out$CGz_sp_man_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[3]),
                        by=list(species = dat_model_out$species), max)
ci_bounds <- merge(min_lb,max_ub, by = c("species"))

dat_CG_plot     <- aggregate(list(mean_CGx_sp_int = dat_model_out$CGx_sp_int,
                                  mean_CGy_sp_int = dat_model_out$CGy_sp_int,
                                  mean_CGz_sp_int = dat_model_out$CGz_sp_int,
                                  mean_CGx_sp_elb = dat_model_out$CGx_sp_elb,
                                  mean_CGy_sp_elb = dat_model_out$CGy_sp_elb,
                                  mean_CGz_sp_elb = dat_model_out$CGz_sp_elb,
                                  mean_CGx_sp_man = dat_model_out$CGx_sp_man,
                                  mean_CGy_sp_man = dat_model_out$CGy_sp_man,
                                  mean_CGz_sp_man = dat_model_out$CGz_sp_man,
                                  mean_CGx_sp_elbman = dat_model_out$CGx_sp_elbman,
                                  mean_CGy_sp_elbman = dat_model_out$CGy_sp_elbman,
                                  mean_CGz_sp_elbman = dat_model_out$CGz_sp_elbman),  by=list(species = dat_model_out$species), mean)

test <- merge(aggregate(list(min_CGx_specific = dat_comp$min_CGx_specific,
                             min_CGy_specific = dat_comp$min_wing_CGy_specific,
                             min_CGz_specific = dat_comp$min_CGz_specific),  by=list(species = dat_comp$species), min),
              aggregate(list(max_CGx_specific = dat_comp$max_CGx_specific,
                             max_CGy_specific = dat_comp$max_wing_CGy_specific,
                             max_CGz_specific = dat_comp$max_CGz_specific),  by=list(species = dat_comp$species), max), by = "species")
test2 <- cbind(stack(test$min_CGx_specific,test$max_CGx_specific),rbind(test$min_CGz_specific,test$max_CGz_specific),rbind(test$species,test$species))
colnames(test2) <- c("x","z","species")
CGx_specific <- ggplot() +
  # add background info
  #geom_polygon(data = test2, aes(x = x, y = y, group = species), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  #geom_errorbarh(data = test, aes(xmin = -min_CGx_specific, xmax = -max_CGx_specific, y = species), height = 0, alpha = 0.5)+
  # geom_point(data = aggregate(list(max_CGx_specific = dat_comp$max_CGx_specific,
  #                                  max_CGz_specific = dat_comp$max_CGz_specific),  by=list(species = dat_comp$species), max),
  #            aes(x = -max_CGx_specific, y = max_CGz_specific), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  # geom_point(data = aggregate(list(min_CGx_specific = dat_comp$min_CGx_specific,
  #                                  min_CGz_specific = dat_comp$min_CGz_specific),  by=list(species = dat_comp$species), min),
  #            aes(x = -min_CGx_specific, y = min_CGz_specific), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGx_specific = dat_comp$mean_CGx_specific,
                                   mean_CGz_specific = dat_comp$mean_CGz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = -mean_CGx_specific, y = mean_CGz_specific, col = species), fill = "white", pch = 21, alpha = 1, size = 1.5) +
  geom_point(data = aggregate(list(shoulderx_specific = dat_bird$shoulderx_specific,
                                   shoulderz_specific = dat_bird$shoulderz_specific),  by=list(species = dat_bird$species), mean),
             aes(x = -shoulderx_specific, y = shoulderz_specific, col = species), pch = 2, alpha = 1) +
  #theme control
  facet_wrap(~species, ncol = 4)+
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "x (% of full body length)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))


CGx_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = test, aes(xmin = -min_CGx_specific, xmax = -max_CGx_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_CGx_specific = dat_comp$max_CGx_specific),  by=list(species = dat_comp$species), max),
             aes(x = -max_CGx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_CGx_specific = dat_comp$min_CGx_specific),  by=list(species = dat_comp$species), min),
             aes(x = -min_CGx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGx_specific = dat_comp$mean_CGx_specific),  by=list(species = dat_comp$species), mean),
             aes(x = -mean_CGx_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  geom_point(data = aggregate(list(shoulderx_specific = dat_bird$shoulderx_specific),  by=list(species = dat_bird$species), mean),
             aes(x = -shoulderx_specific, y = species), col = col_x, pch = 2, alpha = 1) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "x (% of full body length)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

CGy_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = test, aes(xmin = min_CGy_specific, xmax = max_CGy_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_CGy_specific = dat_comp$max_wing_CGy_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_CGy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_CGy_specific = dat_comp$min_wing_CGy_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_CGy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGy_specific = dat_comp$mean_wing_CGy_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_CGy_specific, y = species), fill = "white", col = col_x,pch = 21, alpha = 1) +
  geom_point(data = aggregate(list(shouldery = dat_final$shouldery_specific),  by=list(species = dat_final$species), mean),
             aes(x = shouldery, y = species), fill = col_x, col = "black", pch = 2, alpha = 1) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "y (% of half span)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

CGz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_point(data = aggregate(list(shoulderz_specific = dat_bird$shoulderz_specific),  by=list(species = dat_bird$species), mean),
             aes(x = shoulderz_specific, y = species), fill = col_x, col = "black", pch = 2, alpha = 1) +
  geom_errorbarh(data = test, aes(xmin = min_CGz_specific, xmax = max_CGz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_CGz_specific = dat_comp$max_CGz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_CGz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_CGz_specific = dat_comp$min_CGz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_CGz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGz_specific = dat_comp$mean_CGz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_CGz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "z (% of max. body height)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))



CGx_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1, xmax = 2.5), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.5, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = -pmin(min_CGx_sp_man_lb1,min_CGx_sp_man_lb2), xmax = -pmax(max_CGx_sp_man_ub1,max_CGx_sp_man_ub2), y = species),
                 col = col_z, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = -pmin(min_CGx_sp_elb_lb1,min_CGx_sp_elb_lb2), xmax = -pmax(max_CGx_sp_elb_ub1,max_CGx_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_man+mean_CGx_sp_elbman*elb_fixed[1]), y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_man+mean_CGx_sp_elbman*elb_fixed[2]), y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_man+mean_CGx_sp_elbman*elb_fixed[3]), y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_elb+mean_CGx_sp_elbman*man_fixed[1]), y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_elb+mean_CGx_sp_elbman*man_fixed[2]), y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_elb+mean_CGx_sp_elbman*man_fixed[3]), y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  # to get %/deg we need to divide the output slopes by 10 instead I will just adjust the labels accordingly
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-1,2.5) ,breaks = c(-1,0,1,2), labels = c(-0.1,0,0.1,0.2))+
  geom_rangeframe() +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))


CGy_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1, xmax = 2.5), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.50, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGy_sp_man_lb1,min_CGy_sp_man_lb2), xmax = pmax(max_CGy_sp_man_ub1,max_CGy_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGy_sp_elb_lb1,min_CGy_sp_elb_lb2), xmax = pmax(max_CGy_sp_elb_ub1,max_CGy_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  # add data
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_man+mean_CGy_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_man+mean_CGy_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_man+mean_CGy_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_elb+mean_CGy_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_elb+mean_CGy_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_elb+mean_CGy_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  # to get %/deg we need to divide the output slopes by 10 instead I will just adjust the labels accordingly
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-1,2.5) ,breaks = c(-1,0,1,2), labels = c(-0.1,0,0.1,0.2))+
  geom_rangeframe() +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))


CGz_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1, xmax = 2.4), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.50, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGz_sp_man_lb1,min_CGz_sp_man_lb2), xmax = pmax(max_CGz_sp_man_ub1,max_CGz_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGz_sp_elb_lb1,min_CGz_sp_elb_lb2), xmax = pmax(max_CGz_sp_elb_ub1,max_CGz_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  # add data
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_man+mean_CGz_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_man+mean_CGz_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_man+mean_CGz_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_elb+mean_CGz_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_elb+mean_CGz_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_elb+mean_CGz_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  # to get %/deg we need to divide the output slopes by 10 instead I will just adjust the labels accordingly
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-1,2.4) ,breaks = c(-1,0,1,2), labels = c(-0.1,0,0.1,0.2))+
  geom_rangeframe() +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))


### To plot the maximum wing shapes
# adjust data
max_wing$species_order = factor(max_wing$species, levels = phylo_order$species)
tmpx = melt(max_wing[,c("species_order","pt12_X","pt6_X","pt7_X","pt8_X","pt9_X","pt10_X","pt11_X","edge_X")], id = "species_order")
tmpy = melt(max_wing[,c("species_order","pt12_Y","pt6_Y","pt7_Y","pt8_Y","pt9_Y","pt10_Y","pt11_Y","edge_Y")], id = "species_order")
wing_shape = cbind(tmpx,tmpy[,c(2,3)])
colnames(wing_shape) <- c("species_order","variable_x","value_x","variable_y","value_y")










### ----------- Gull only data ---------------
gull_model_Ixx <- lm(full_Ixx~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_model_Iyy <- lm(full_Iyy~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_model_Izz <- lm(full_Izz~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_model_Ixz <- lm(full_Ixz~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_Ixx <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Ixx)) + th

gull_Iyy <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Iyy)) + th

gull_Izz <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Izz)) + th

gull_Ixz <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Ixz)) + th


## ----------------------------------

validation_Izz_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = wing_Izz, col = species_order)) +
  geom_line(data = unique(dat_final[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_final[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Ixx_body_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = full_Ixx, col = species_order)) +
  geom_point(aes(x = 0.0875, y = 928/(100^2*1000)), col = "black") + # Hedrick and Biewener
  geom_point(aes(x = 0.2893, y = 12889/(100^2*1000)), col = "black") + # Hedrick and Biewener
  th +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')


test <- lm(prop_q_dot_nd ~ elbow*manus + species_order, data = dat_final)
xgrid <-  seq(floor(min(dat_final$elbow)), ceiling(max(dat_final$elbow)), 1)
ygrid <-  seq(floor(min(dat_final$manus)), ceiling(max(dat_final$manus)), 1)
zgrid <-  phylo_order$species
data.fit       <- expand.grid(elbow = xgrid, manus = ygrid, species_order = zgrid)
data.fit$prop_q  <-  predict(test, newdata = data.fit)
