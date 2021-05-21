library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe
library(gridExtra) # for using grid arrange
library(ggtree)    # for plotting phylogeny
library(alphahull) # for computing convex hulls
library(plyr)      # for computing convex hulls
library(reshape)   # for melt
library(sensiPhy)
library(geomorph)
source("plotting_info.R")
source("plotting_dataprep.R")

## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 1 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------
#exported as 5x5
validation_Ixx_plot <- ggplot()+
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = 3.76*10^-3*full_m^1.89, ymax = 3.76*10^-3*full_m^2.22),
              col = NA, fill = "gray70", alpha = 0.15) + # Kirkpatrick, 1994, 99% fiducity limits
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = 1.94*10^-3*full_m^1.787, ymax = 1.94*10^-3*full_m^2.134),
              col = NA, fill = "gray50", alpha = 0.2) + # Berg and Rayner, 1995, 95% fiducity limits
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Ixx_val_mcmc_output$solutions[1,1])*full_m^Ixx_val_mcmc_output$solutions[2,2], ymax = exp(Ixx_val_mcmc_output$solutions[1,1])*full_m^Ixx_val_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.2) + # AvInertia
  #geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_wing_Ixx, ymin = sachs_pred_Ixx, col = species_order), alpha = 0.5) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "gray70", linetype = 3) + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray50", linetype = 2) + # Berg and Rayner, 1995
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Ixx_val_mcmc_output$solutions[1,1])*full_m^Ixx_val_mcmc_output$solutions[2,1]),
            col = "black") + # AvInertia
  #geom_point(data = dat_comp, aes(x = full_m, y = sachs_pred_Ixx, fill = species_order, col = species_order), pch = 22, size = 1.5) + # Sachs, 2005
  geom_point(data = dat_comp, aes(x = full_m, y = max_wing_Ixx, col = species_order), pch = 19, size = 2, alpha = 0.6)+
  geom_point(data = dat_comp, aes(x = full_m, y = max_wing_Ixx, col = species_order), pch = 1, size = 2)+
  # theme control
  th  +
  theme(legend.position="none") +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #axis control
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Maximum wing I"["xx"]," (kg-m"^2,")", sep = "")),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

#exported as 3.5x10
ROM_plot <- ggplot()+
  geom_polygon(data = chulls_elbman, aes(x = elbow, y = manus, fill = species_order), col = NA) +
  geom_point(data = subset(dat_final, species == "col_liv" & BirdID == "20_0300" & TestID == "20_0300_4" & FrameID == 317), aes(x = elbow, y = manus, fill = species_order),size = 0.5) +
  geom_point(data = subset(dat_final, species == "col_liv" & BirdID == "20_0300" & TestID == "20_0300_2" & FrameID == 94), aes(x = elbow, y = manus, fill = species_order),size = 0.5) +
  facet_wrap(~species_order, nrow = 3) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA))+
  theme(legend.position="none") +
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)', breaks = c(50,150)) +
  scale_y_continuous(name='Wrist angle (°)', breaks = c(50,150)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = 50, yend = 150) +
  annotate(geom = "segment", x = 50, xend = 150, y = -Inf, yend = -Inf) +
  annotate(geom = "segment", x = 50, xend = 50, y = -Inf, yend = 50) +
  annotate(geom = "segment", x = 150, xend = 150, y = -Inf, yend = 50) +
  annotate(geom = "segment", x = -Inf, xend = 30, y = 50, yend = 50) +
  annotate(geom = "segment", x = -Inf, xend = 30, y = 150, yend = 150)

## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 2 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------

# Panel A plot
CGxz_plot <- ggplot()+
  # add data
  geom_rect(data = vertices_cgxz, aes(xmin =0, ymin = -0.1, xmax = -0.2, ymax = 0, group = species_order), linetype="dashed", color = "gray", fill = NA, alpha = 0.25) +
  # shoulder motion polygons
  geom_rect(data = cg_x_shape[which(cg_x_shape$variable == "mean_CGxback"),], aes(xmin = value, xmax = cg_x_shape$value[which(cg_x_shape$variable == "mean_CGxfore")],
                ymin = -cg_z_shape$value[which(cg_z_shape$variable == "mean_CGzup")],  ymax = -cg_z_shape$value[which(cg_z_shape$variable == "mean_CGzdn")],
                fill = species_order), col = NA, alpha = 0.4) +
  # elbow and wrist polygons
  geom_polygon(data = vertices_cgxz, aes(x = full_CGx_specific_orgShoulder, y = -full_CGz_specific_orgShoulder, col = species_order, fill = species_order, group = BirdID)) +
  facet_wrap(~species_order, ncol = 1) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  th +
  #theme control
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(legend.position="none") +
  coord_fixed() +
  scale_x_continuous(name='x (%)', limits = c(-0.2,0.05), breaks = c(-0.2,-0.1,0), labels = c(-20,-10,0)) +
  scale_y_continuous(name='z (%)', limits = c(-0.12,0.08), breaks = c(-0.1,0), labels = c(10,0)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.1, yend = 0) +
  annotate(geom = "segment", x = 0, xend = -0.2, y = log(0), yend = log(0))

## -------------------- Panel B ------------------------

# plot
wingCGxy_plot <- ggplot() +
  geom_rect(data = vertices_cgxy, aes(xmin =0, ymin = -0.14, xmax = 0.3, ymax = 0, group = species_order), linetype="dashed", color = "gray", fill = NA, alpha = 0.25) +
  #geom_path(data = wing_shape, aes(x = value_y, y = value_x), col = "black", alpha = 0.2) +
  geom_polygon(data = vertices_cgxy, aes(x = wing_CGy_specific_orgShoulder, y = wing_CGx_specific_orgShoulder, col = species_order, fill = species_order, group = BirdID)) +
  facet_wrap(~species_order, ncol = 1) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(legend.position="none") +
  coord_fixed() +
  scale_x_continuous(name='y (%)', limits = c(-0.01,0.31), breaks = c(0,0.1,0.2,0.3), labels = c(0,10,20,30)) +
  scale_y_continuous(name='x (%)', limits = c(-0.15,0.1), breaks = c(-0.1,0), labels = c(-10,0)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.1, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 0.3, y = log(0), yend = log(0))


effect_CGx <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(CGx_elb_etap = dat_model_out$CGx_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGx_man_etap = dat_model_out$CGx_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGx_elbman_etap = dat_model_out$CGx_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +  # axis control
  scale_y_continuous(limits = c(0,12), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))

effect_CGz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(CGz_elb_etap = dat_model_out$CGz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGz_man_etap = dat_model_out$CGz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGz_elbman_etap = dat_model_out$CGz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))

effect_CGy <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(CGy_elb_etap = dat_model_out$CGy_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGy_man_etap = dat_model_out$CGy_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGy_elbman_etap = dat_model_out$CGy_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))


CG_isometry_fullbird <- ggplot() +
  # add specific line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,1]),
            col = "gray") +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray", alpha = 0.3) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,1]),
            col = "gray") +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray", alpha = 0.3) +
  # add line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^CGx_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^CGx_model_mcmc_output$solutions[2,2], ymax = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^CGx_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.3) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^CGz_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^CGz_model_mcmc_output$solutions[2,2], ymax = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^CGz_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.3) +
  # add specific data
  #geom_errorbar(data = dat_comp, aes(x = full_m, ymin = -min_CGx_specific, ymax = -max_CGx_specific), width = 0.08, col = "gray", alpha = 0.4) +
  geom_point(data = dat_comp, aes (x = full_m, y = -mean_CGx_specific_orgBeak), col = "gray", pch = 0) +
  #geom_errorbar(data = dat_comp, aes(x = full_m, ymin = min_CGz_specific, ymax = max_CGz_specific), width = 0.08, col = "gray", alpha = 0.4) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_CGz_specific_orgDorsal), col = "gray", pch = 2) +
  # add exact data
  geom_point(data = dat_comp, aes (x = full_m, y = -mean_CGx_orgBeak), col = "black", pch = 15) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_CGz_orgDorsal), col = "black", pch = 17) +
  # add iso lines
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^(1/3)), col = "black", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^(1/3)), col = "black", linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean CG position (m)", sep = "")), limits = c(0.005,1),
                     breaks = c(1E-2,1E-1,1E0), labels = c(expression(10^-2),expression(10^-1),expression(10^0)),
                     sec.axis = sec_axis(~.*100, name= expression(paste("Mean ", bar("CG"), " (%)", sep = "")), breaks = c(1,10,30,100)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-2, yend = 1E0) +
  annotate(geom = "segment", x = Inf, xend = Inf, y = 1E-2, yend = 1E0, col = "gray") +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

CG_isometry_singlewing <- ggplot() +
  # add specific line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,1]),
            col = "gray") +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray", alpha = 0.3) +
  # add line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^CGy_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^CGy_model_mcmc_output$solutions[2,2], ymax = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^CGy_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.3) +
  # add data
  geom_point(data = dat_comp, aes (x = full_m, y = mean_wing_CGy), col = "black", pch = 16) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_wing_CGy_specific), col = "gray", pch = 1) +
  # add iso lines
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^(1/3)), col = "black", linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean CGy"["wing"], " position (m)", sep = "")), limits = c(0.005,1),
                     breaks = c(1E-2,1E-1,1E0), labels = c(expression(10^-2),expression(10^-1),expression(10^0)),
                     sec.axis = sec_axis(~.*100, name= expression(paste("Mean ", bar("CGy")["wing"], " (%)", sep = "")), breaks = c(1,10,30,100)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-2, yend = 1E0) +
  annotate(geom = "segment", x = Inf, xend = Inf, y = 1E-2, yend = 1E0, col = "gray") +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)


# -------- Combine into panel ------------

leftcol <- plot_grid(phylo_plot_complete,CGxz_plot, wingCGxy_plot,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1,0.6,0.6),
                       #labels
                       labels = c("","",""),
                       label_size = 10,
                       label_fontfamily = "sans")

#export as 9x3.5 if inputting independently
rightcol_ind <- plot_grid(effect_CGx,effect_CGz,CG_isometry_fullbird,effect_CGy,CG_isometry_singlewing,
                      #arrangement data
                      ncol = 1,
                      rel_heights = c(0.5,0.5,1.3,0.5,1.3),
                      #labels
                      labels = c("",""),
                      label_size = 10,
                      label_fontfamily = "sans")

#exported as 12x12
Figure2_final <- plot_grid(leftcol, rightcol_ind,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(1,0.8),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")

## ------------------------------------------------------------------------------------
## --------------------------------- Figure 3 -----------------------------------------
## ------------------------------------------------------------------------------------

I_ratio <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  #geom_errorbarh(data = dat_I_sp_plot, aes(xmin = -0.5*(max_Ixx-min_Iyy), xmax = 0.5*(max_Ixx-min_Iyy), y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(Ixx_range =  dat_comp$max_Ixx/dat_comp$min_Ixx),  by=list(species = dat_comp$species), max),
             aes(x = Ixx_range, y = species), col = "black",alpha = 0.85, size = 1, pch = 0) +
  geom_point(data = aggregate(list(Iyy_range = dat_comp$max_Iyy/dat_comp$min_Iyy),  by=list(species = dat_comp$species), max),
             aes(x = Iyy_range, y = species), col = "gray40",alpha = 0.85, size = 1, pch = 1) +
  geom_point(data = aggregate(list(Izz_range =  dat_comp$max_Izz/dat_comp$min_Izz),  by=list(species = dat_comp$species), max),
             aes(x = Izz_range, y = species), col = "gray60",alpha = 0.85, size = 1, pch = 2) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  scale_x_continuous(name = "Maximum I/Minimum I", limits= c(1,9), breaks = c(1,3,5,7,9)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 1, xend = 9, y = log(0), yend = log(0))

effect_Ixx <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Ixx_elb_etap = dat_model_out$Ixx_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixx_man_etap = dat_model_out$Ixx_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixx_elbman_etap = dat_model_out$Ixx_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,10), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))

effect_Iyy <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Iyy_elb_etap = dat_model_out$Iyy_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Iyy_man_etap = dat_model_out$Iyy_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Iyy_elbman_etap = dat_model_out$Iyy_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,10), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))
effect_Izz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Izz_elb_etap = dat_model_out$Izz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Izz_man_etap = dat_model_out$Izz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Izz_elbman_etap = dat_model_out$Izz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +  # axis control
  scale_y_continuous(limits = c(0,10), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))
effect_Ixz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Ixz_elb_etap = dat_model_out$Ixz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixz_man_etap = dat_model_out$Ixz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixz_elbman_etap = dat_model_out$Ixz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,10), name = "", breaks = c(0,5,10)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 10, x = log(0), xend = log(0))

I_isometry <- ggplot() +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^Ixx_model_mcmc_output$solutions[2,2], ymax = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^Ixx_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.15) +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^Iyy_model_mcmc_output$solutions[2,2], ymax = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^Iyy_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray40", alpha = 0.15) +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Izz_model_mcmc_output$solutions[2,2], ymax = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Izz_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray60", alpha = 0.15) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^Ixx_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^Iyy_model_mcmc_output$solutions[2,1]),
            col = "gray40") +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Izz_model_mcmc_output$solutions[2,1]),
            col = "gray60") +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "black", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "gray40", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "gray60", linetype = 2) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_Ixx), col = "black", pch = 0, size = 0.8) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_Iyy), col = "gray40", pch = 1, size = 0.8)+
  geom_point(data = dat_comp, aes (x = full_m, y = mean_Izz), col = "gray60", pch = 2, alpha = 0.8, size = 0.8)+
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean moment of inertia (kg-m"^2,")", sep = "")), limits = c(1E-6,2E-1),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)


q_plot <- ggplot() +
  # add model fit
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(max_q_model_mcmc_output$solutions[1,1])*full_m^max_q_model_mcmc_output$solutions[2,1]),
            col = "gray60") +
  #add data
  geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_q, ymin = min_q, col = species_order), alpha = 0.5) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q, col = species_order), alpha = 0.6, pch = 19,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q, col = species_order), pch = 1,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_q, col = species_order), alpha = 0.6, pch = 19,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_q, col = species_order), pch = 1,size = 0.8) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10), name = "Mass (kg)",
                     labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', limits = c(0.1,18), breaks = c(0.1,1,10),
                     labels = c(expression(10^-1),expression(10^0),expression(10^1)), name = expression(paste("Pitch agility (s"^{-2},")")))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-1, yend = 1E1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

q_nd_plot <- ggplot() +
  #add data
  geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_q_nd, ymin = min_q_nd, col = species_order), alpha = 0.5) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q_nd, col = species_order), alpha = 0.6, pch = 19,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q_nd, col = species_order), pch = 1,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_q_nd, col = species_order), alpha = 0.6, pch = 19,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_q_nd, col = species_order), pch = 1,size = 0.8) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10), name = "Mass (kg)",
                     labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', limits = c(0.1,1.4), breaks = c(0.1,1),
                     labels = c(expression(10^-1),expression(10^0)), name = expression(paste("Normalized pitch agility")))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-1, yend = 1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

# -------- Combine panels into figure ------------
leftcol <- plot_grid(I_isometry,q_plot,q_nd_plot,
                       #arrangement data
                       ncol = 1,
                       rel_heights = c(1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")
rightcol <- plot_grid(I_ratio,blank_plot,
                       #arrangement data
                       ncol = 1,
                       rel_heights = c(2,1.1),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")
#exported as 7.5x6
toprow <- plot_grid(leftcol,blank_plot, rightcol,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1,0.25,0.75),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")
#exported as 4x4
effectcol <- plot_grid(effect_Ixx,effect_Iyy,effect_Izz,effect_Ixz,
                      #arrangement data
                      ncol = 1,
                      rel_heights = c(1,1,1,1),
                      #labels
                      labels = c("",""),
                      label_size = 10,
                      label_fontfamily = "sans")
## ------------------------------------------------------------------------------------
## --------------------------------- Figure 4 -----------------------------------------
## ------------------------------------------------------------------------------------

var_test_plot <- ggplot()+
  geom_point(data = subset(dat_var_tot, component != "full_I"),
             aes(x = mean_sig_sq, y = 1,fill = mean_sig_sq), pch = 21, col = "black") +
  geom_point(data = dat_var_tot, aes(x = ind_var_CGx$opt$sigsq, y = 1,fill = ind_var_CGx$opt$sigsq), col = "black", size = 2, pch = 22) +
  geom_point(data = dat_var_tot, aes(x = ind_var_CGz$opt$sigsq, y = 1,fill = ind_var_CGz$opt$sigsq), col = "black", size = 2, pch = 22) +
  geom_point(data = dat_var_tot, aes(x = ind_var_CGx_sh$opt$sigsq, y = 1,fill = ind_var_CGx_sh$opt$sigsq), col = "black", size = 2, pch = 22) +
  geom_point(data = dat_var_tot, aes(x = ind_var_CGz_sh$opt$sigsq, y = 1,fill = ind_var_CGz_sh$opt$sigsq), col = "black", size = 2, pch = 22) +
  scale_fill_gradient(low = "white", high = "#2B2530", name = "rate of variance evolution") +
  th

# Extract the correct colours
merge(subset(ggplot_build(var_test_plot)$data[[1]], shape ==21),subset(dat_var_tot, component != "full_I"), by.x = "x", by.y = "mean_sig_sq")

# exported as 5x5 - 5.5x4
var_sens_plot <- ggplot()+
  geom_dotplot(data = dat_var_range, aes(x = component, y = mean_sig_sq), col = NA, fill = "black", binaxis='y', stackdir='center', dotsize = 8, binwidth = 1/400, alpha = 0.6) +
  geom_point(data = dat_var_tot, aes(x = "full_CG", y = ind_var_CGx$opt$sigsq),fill = "gray",  col = "black", size = 1.5, pch = 22, alpha = 0.3) +
  geom_point(data = dat_var_tot, aes(x = "full_CG", y = ind_var_CGz$opt$sigsq),fill = "gray",  col = "black", size = 1.5, pch = 24, alpha = 0.3) +
  geom_point(data = dat_var_tot, aes(x = "full_CG_sh", y = ind_var_CGx_sh$opt$sigsq),fill = "gray",  col = "black", size = 1.5, pch = 22, alpha = 0.3) +
  geom_point(data = dat_var_tot, aes(x = "full_CG_sh", y = ind_var_CGz_sh$opt$sigsq),fill = "gray",  col = "black", size = 1.5, pch = 24, alpha = 0.3) +
  geom_point(data = dat_var_tot, aes(x = component, y = mean_sig_sq), fill = "gray", col = "black", size = 2.5, pch = 23) +
  th+
  scale_y_continuous(trans = "log10", limits = c(0.005,0.12), breaks = c(0.005,0.01,0.05,0.1), name = expression(paste(sigma^2,""[mult]))) +
  geom_rangeframe() +
  annotate(geom = "segment", y = 0.005, yend = 0.1, x = log(0), xend = log(0)) +
  coord_flip()

#exported as 11x3.5
chord_length_plot <- ggplot()+
  geom_point(data = stab_check, aes(x = 1, y = 1, col = log(c_l_max), size = max_q_nd*1.5), alpha = 0.5) +
  geom_point(data = stab_check, aes(x = 1, y = 1, col = log(c_l_max), size = max_q_nd*1.5), pch = 1, alpha = 1) +
  geom_point(data = stab_check, aes(x = 0.98, y = 1, col = log(c_l_min), size = min_q_nd*1.5), alpha = 0.5) +
  geom_point(data = stab_check, aes(x = 0.98, y = 1, col = log(c_l_min), size = min_q_nd*1.5), pch = 1, alpha = 1) +
  th +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA))+
  facet_wrap(~species_order, ncol = 1) +
  scale_size(range = c(0.05, 10), name="Normalized pitch agility") + #scales based on area
  scale_colour_gradient2(midpoint=0, low="#1D0747", mid="#FFF9D3", high="#00510A", space ="Lab", name = "Estimated static margin ratio") +
  scale_y_continuous(limits = c(0.955,1.025)) +
  scale_x_continuous(limits = c(0.955,1.025))
