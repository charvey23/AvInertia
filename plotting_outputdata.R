library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe
library(gridExtra) # for using grid arrange


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
    # Background transparency
    # Background of panel
    panel.background = ggplot2::element_rect(fill = "transparent"),
    # Background behind actual data points
    plot.background  = ggplot2::element_rect(fill = "transparent", color = NA))
#
# ## TO CHECK DATA
# x_vertices = c(dat_wing$pt1_X[i], dat_wing$pt6_X[i],dat_wing$pt7_X[i],dat_wing$pt8_X[i],dat_wing$pt9_X[i],dat_wing$pt10_X[i],dat_wing$pt11_X[i],dat_wing$pt12_X[i])
# y_vertices = c(dat_wing$pt1_Y[i], dat_wing$pt6_Y[i],dat_wing$pt7_Y[i],dat_wing$pt8_Y[i],dat_wing$pt9_Y[i],dat_wing$pt10_Y[i],dat_wing$pt11_Y[i],dat_wing$pt12_Y[i])
# z_vertices = c(dat_wing$pt1_Z[i], dat_wing$pt6_Z[i],dat_wing$pt7_Z[i],dat_wing$pt8_Z[i],dat_wing$pt9_Z[i],dat_wing$pt10_Z[i],dat_wing$pt11_Z[i],dat_wing$pt12_Z[i])
#
#
# x = c(dat_wing$pt1_X[m],dat_wing$pt2_X[m],dat_wing$pt3_X[m],dat_wing$pt4_X[m],
#       dat_wing$pt6_X[m],dat_wing$pt7_X[m],dat_wing$pt8_X[m],dat_wing$pt9_X[m],
#       dat_wing$pt10_X[m],dat_wing$pt11_X[m],dat_wing$pt12_X[m])
# y = c(dat_wing$pt1_Y[m],dat_wing$pt2_Y[m],dat_wing$pt3_Y[m],dat_wing$pt4_Y[m],
#       dat_wing$pt6_Y[m],dat_wing$pt7_Y[m],dat_wing$pt8_Y[m],dat_wing$pt9_Y[m],
#       dat_wing$pt10_Y[m],dat_wing$pt11_Y[m],dat_wing$pt12_Y[m])
# z = c(dat_wing$pt1_Z[m],dat_wing$pt2_Z[m],dat_wing$pt3_Z[m],dat_wing$pt4_Z[m],
#       dat_wing$pt6_Z[m],dat_wing$pt7_Z[m],dat_wing$pt8_Z[m],dat_wing$pt9_Z[m],
#       dat_wing$pt10_Z[m],dat_wing$pt11_Z[m],dat_wing$pt12_Z[m])
# plot3d(x, y, z, col = c("black","gray20","gray50","gray70","red", "orange", "yellow", "green", "blue", "navy", "purple"))

CG_side_view <- ggplot()+
  geom_point(data = dat_bird, aes(x = -shoulderx_specific, y = -shoulderz_specific, col = log(full_m)), pch = 2) +
  geom_point(data = dat_final, aes(x = -full_CGx_specific, y = -full_CGz_specific, col = log(full_m)), pch = 10) + th +
  scale_y_continuous(name = "z (m)",limits = c(-0.15,0.15)) +
  scale_x_continuous(name = "x (m)",limits = c(0,1)) +
  coord_fixed() + facet_wrap(~species, ncol = 2)

CG_top_view <- ggplot()+
  geom_point(data = dat_final, aes(y = wing_CGx_specific, x = wing_CGy_specific, col = elbow), pch = 10) + th +
  scale_y_continuous(name = "x (m)",limits = c(-0.15,0.15)) +
  scale_x_continuous(name = "y (m)",limits = c(0,1)) +
  coord_fixed() + facet_wrap(~species, ncol = 2)

CG_front_view <- ggplot()+
  geom_point(data = dat_final, aes(y = wing_CGz_specific, x = wing_CGy_specific, col = log(full_m)), pch = 10) + th +
  scale_y_continuous(name = "x (m)",limits = c(-0.15,0.15)) +
  scale_x_continuous(name = "y (m)",limits = c(0,1)) +
  coord_fixed() + facet_wrap(~species, ncol = 2)


cohens_effect <- ggplot()+
  geom_density(data = dat_model_out, aes(x = CGx_elb_eff), fill = "blue")


ROM_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = elbow, y = manus, col = abs(full_Ixz/full_Izz))) + th + facet_wrap(~species, ncol = 2)

q_dot_plot <- ggplot()+
  geom_point(data = subset(dat_final, species != "oce_leu"), aes(x = full_m, y = prop_q_dot_nd, col = clade), pch = 15) + th

del_M_plot <- ggplot()+
  geom_point(data = subset(dat_all_e, BirdID != "21_0203" &species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade)) +
  geom_point(data = subset(dat_all_t, BirdID != "21_0203" &species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade), pch = 15) +
  th +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

cg_plot <- ggplot()+
  geom_point(data = subset(dat_final, species == "col_liv"), aes(x = manus, y = -full_CGx_specific, col = elbow)) + th

CG_range_plot <- ggplot()+
  geom_point(data = dat_bird, aes(x = total_bird_mass, y = range_CGx_specific, col = clade)) + th +
  scale_x_continuous(trans='log10')

Ixz_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = full_m, y = wing_Ixz/(full_m*torso_length^2), col = species)) + th


## ----------------- Plot of linear model outputs ------------------
# order species to match the phylogeny
phylo_order = c("bra_can","ana_pla", "lop_imp", "lop_nyc","chr_amh","col_liv","cho_min","cyp_nig","pel_ery","ard_her",
                "col_aur","tyt_alb","acc_str","acc_coo","fal_per","fal_col","meg_alc","cya_ste","cor_cor")

shading <- data.frame(col1 = levels(dat_int_model$species)[seq(from = 2, to = max(as.numeric(as.factor(dat_int_model$species)))-1, by = 2)],
                      col2 = levels(dat_int_model$species)[seq(from = 3, to = max(as.numeric(as.factor(dat_int_model$species))), by = 2)])
man_fixed = c(80,100,120)*0.001
elb_fixed = c(80,100,120)*0.001
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
# incase af0069

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

CG_int_effect <- ggplot() +
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = min_CGx_sp_int_lb, xmax = max_CGx_sp_int_ub, y = species),
                 col = col_x, position = position_nudge(y = 0.7), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = min_CGy_sp_int_lb, xmax = max_CGy_sp_int_ub, y = species),
                 col = col_y, position = position_nudge(y = 0.5), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = -min_CGz_sp_int_lb, xmax = -max_CGz_sp_int_ub, y = species),
                 col = col_z, position = position_nudge(y = 0.3), alpha = 0.2, height = 0.3, size = 1) +
  # add data
  geom_point(data = dat_model_out, aes(x = CGx_sp_int, y = species), col = col_x, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_int, y = species), col = col_y, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_int, y = species), col = col_z, position = position_nudge(y = 0.3)) +
  #theme control
  th +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order), name = "") +
  scale_x_continuous(name = expression(paste("Intercept x 10"^"3"," (%)", sep ="")), limits= c(-0.5,0.3), breaks = c(-0.5,-0.25,0,0.25), labels = c(-50,-25,0,25)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 1, yend = 19) +
  annotate(geom = "segment", x = -0.5, xend = 0.25, y = log(0), yend = log(0))

CG_elb_effect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1) +
  geom_segment(aes(y = 1, yend = 19, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGx_sp_elb_lb1,min_CGx_sp_elb_lb2), xmax = pmax(max_CGx_sp_elb_ub1,max_CGx_sp_elb_ub2), y = species),
                 col = col_x, position = position_nudge(y = 0.7), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGy_sp_elb_lb1,min_CGy_sp_elb_lb2), xmax = pmax(max_CGy_sp_elb_ub1,max_CGy_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = 0.5), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = -pmin(min_CGz_sp_elb_lb1,min_CGz_sp_elb_lb2), xmax = -pmax(max_CGz_sp_elb_ub1,max_CGz_sp_elb_ub2), y = species),
                 col = col_z, position = position_nudge(y = 0.3), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_model_out, aes(x = CGx_sp_elb+CGx_sp_elbman*man_fixed[1], y = species), col = col_x, alpha = 0.4, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGx_sp_elb+CGx_sp_elbman*man_fixed[2], y = species), col = col_x, alpha = 0.7, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGx_sp_elb+CGx_sp_elbman*man_fixed[3], y = species), col = col_x, alpha = 1, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_elb+CGy_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_elb+CGy_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_elb+CGy_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_elb-CGz_sp_elbman*man_fixed[1], y = species), col = col_z, alpha = 0.4, position = position_nudge(y = 0.3)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_elb-CGz_sp_elbman*man_fixed[2], y = species), col = col_z, alpha = 0.7, position = position_nudge(y = 0.3)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_elb-CGz_sp_elbman*man_fixed[3], y = species), col = col_z, alpha = 1, position = position_nudge(y = 0.3)) +
  #theme control
  th +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order), name = "")+
  scale_x_continuous(name = expression(paste("Elbow Coefficient of Regression x 10"^"3"," (%/°)", sep ="")), limits= c(-2.1,2.1)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 1, yend = 19) +
  annotate(geom = "segment", x = -2, xend = 2, y = log(0), yend = log(0))

CG_man_effect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1) +
  geom_segment(aes(y = 1, yend = 19, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGx_sp_man_lb1,min_CGx_sp_man_lb2), xmax = pmax(max_CGx_sp_man_ub1,max_CGx_sp_man_ub2), y = species),
                 col = col_x, position = position_nudge(y = 0.7), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGy_sp_man_lb1,min_CGy_sp_man_lb2), xmax = pmax(max_CGy_sp_man_ub1,max_CGy_sp_man_ub2), y = species),
                 col = col_y, position = position_nudge(y = 0.5), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = -pmin(min_CGz_sp_man_lb1,min_CGz_sp_man_lb2), xmax = -pmax(max_CGz_sp_man_ub1,max_CGz_sp_man_ub2), y = species),
                 col = col_z, position = position_nudge(y = 0.3), alpha = 0.2, height = 0.3, size = 1) +
  # add data
  geom_point(data = dat_model_out, aes(x = CGx_sp_man+CGx_sp_elbman*elb_fixed[1], y = species), col = col_x, alpha = 0.2, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGx_sp_man+CGx_sp_elbman*elb_fixed[2], y = species), col = col_x, alpha = 0.5, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGx_sp_man+CGx_sp_elbman*elb_fixed[3], y = species), col = col_x, alpha = 1, position = position_nudge(y = 0.7)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_man+CGy_sp_elbman*elb_fixed[1], y = species), col = col_y, alpha = 0.2, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_man+CGy_sp_elbman*elb_fixed[2], y = species), col = col_y, alpha = 0.5, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = CGy_sp_man+CGy_sp_elbman*elb_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = 0.5)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_man-CGz_sp_elbman*elb_fixed[1], y = species), col = col_z, alpha = 0.2, position = position_nudge(y = 0.3)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_man-CGz_sp_elbman*elb_fixed[2], y = species), col = col_z, alpha = 0.5, position = position_nudge(y = 0.3)) +
  geom_point(data = dat_model_out, aes(x = -CGz_sp_man-CGz_sp_elbman*elb_fixed[3], y = species), col = col_z, alpha = 1, position = position_nudge(y = 0.3)) +
  # theme control
  th +
  # axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order), name = "") +
  scale_x_continuous(name = expression(paste("Wrist Coefficient of Regression x 10"^"3"," (%/°)", sep ="")), limits= c(-2.1,2.1))+
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 1, yend = 19) +
  annotate(geom = "segment", x = -2, xend = 2, y = log(0), yend = log(0))

# -------- Combine panels into figure ------------
bottomrow <- plot_grid(CG_int_effect,CG_elb_effect, CG_man_effect,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1,1,1),
                       #labels
                       labels = c("A","B","C"),
                       label_size = 10,
                       label_fontfamily = "sans")

bottomrow <- plot_grid(CG_man_effect,CG_man_effect1,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(1,1,1),
                       #labels
                       labels = c("A","B","C"),
                       label_size = 10,
                       label_fontfamily = "sans")
#exported as 7.5x9
figure_final <- plot_grid(toprow,bottomrow,
                           #arrangement data
                           ncol = 1, nrow = 2, rel_heights = c(1,1))



## ----------------- Plots for sharing with others -----------------

validation_Ixx_plot <- ggplot()+
  geom_point(data = dat_all_e, aes(x = full_m, y = wing_Ixx, col = clade)) +
  geom_line(data = unique(dat_all_e[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_all_e[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_all_e[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Izz_plot <- ggplot()+
  geom_point(data = dat_all, aes(x = full_m, y = wing_Izz, col = clade)) +
  geom_line(data = unique(dat_all[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_all[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_all[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Ixx_body_plot <- ggplot()+
  geom_point(data = dat_all, aes(x = full_m, y = full_Ixx, col = clade)) +
  geom_point(aes(x = 0.0875, y = 928/(100^2*1000)), col = "black") + # Hedrick and Biewener
  geom_point(aes(x = 0.2893, y = 12889/(100^2*1000)), col = "black") + # Hedrick and Biewener
  th +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')
