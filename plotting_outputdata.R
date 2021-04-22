library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe
library(gridExtra) # for using grid arrange
library(ggtree)    # for plotting phylogeny
library(alphahull) # for computing convex hulls
library(plyr)      # for computing convex hulls
library(reshape)   # for melt

source("plotting_info.R")
source("plotting_dataprep.R")

## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 1 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------

validation_Ixx_plot <- ggplot()+
  geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_wing_Ixx, ymin = sachs_pred_Ixx, col = species_order), alpha = 0.5) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "gray") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray", linetype = 2) + # Berg and Rayner, 1995
  geom_point(data = dat_comp, aes(x = full_m, y = sachs_pred_Ixx, fill = species_order, col = species_order), pch = 22, size = 1.5) + # Sachs, 2005
  geom_point(data = dat_comp, aes(x = full_m, y = max_wing_Ixx, fill = species_order), pch = 21, col = "black", size = 2)+
  th  +
  theme(legend.position="none") +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Maximum wing I"["xx"]," (kg-m"^2,")", sep = "")),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

ROM_plot <- ggplot()+
  geom_polygon(data = chulls_elbman, aes(x = elbow, y = manus, fill = species_order), col = NA) +
  geom_point(data = subset(dat_final, species == "col_liv" & BirdID == "20_0300" & TestID == "20_0300_4" & FrameID == 317), aes(x = elbow, y = manus, fill = species_order),size = 0.5) +
  geom_point(data = subset(dat_final, species == "col_liv" & BirdID == "20_0300" & TestID == "20_0300_2" & FrameID == 94), aes(x = elbow, y = manus, fill = species_order),size = 0.5) +
  facet_wrap(~species_order, nrow = 3) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  #theme control
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.line =  element_line(colour = "black", linetype=1))+
  theme(legend.position="none") +
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)', breaks = c(50,150)) +
  scale_y_continuous(name='Wrist angle (°)', breaks = c(50,150))

## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 2 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------

# Panel A plot
CGxz_plot <- ggplot()+
  # add data
  geom_rect(data = vertices_cgxz, aes(xmin =0, ymin = -0.1, xmax = 0.2, ymax = 0, group = species_order), linetype="dashed", color = "gray", fill = NA, alpha = 0.25) +
  # geom_path(data = x_shape, aes(x = -value, y = -z_shape$value), col = "black", alpha = 0.2) +
  geom_polygon(data = cg_x_shape, aes(x = -value, y = -cg_z_shape$value, fill = species_order), col = NA, alpha = 0.5) +
  geom_polygon(data = vertices_cgxz, aes(x = -full_CGx_specific_orgShoulder, y = -full_CGz_specific_orgShoulder, col = species_order, fill = species_order, group = BirdID)) +
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
  scale_x_continuous(name='x (%)', limits = c(-0.04,0.2), breaks = c(0,0.1,0.2), labels = c(0,-10,-20)) +
  scale_y_continuous(name='z (%)', limits = c(-0.12,0.08), breaks = c(-0.1,0), labels = c(10,0)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.1, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 0.2, y = log(0), yend = log(0))

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
  scale_y_continuous(name='x (%)', limits = c(-0.15,0.1), breaks = c(-0.1,0,0.1), labels = c(-10,0,10)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.1, yend = 0.1) +
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
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

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
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

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
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12)) +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

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
  geom_point(data = dat_comp, aes (x = full_m, y = -mean_CGx_specific), col = "gray", pch = 19) +
  #geom_errorbar(data = dat_comp, aes(x = full_m, ymin = min_CGz_specific, ymax = max_CGz_specific), width = 0.08, col = "gray", alpha = 0.4) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_CGz_specific), col = "gray", pch = 1) +
  # add exact data
  geom_point(data = dat_comp, aes (x = full_m, y = -mean_CGx_orgBeak), col = "black", pch = 19) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_CGz_orgDorsal), col = "black", pch = 1) +
  # add iso lines
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^(1/3)), col = "black", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^(1/3)), col = "black", linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean center of gravity (m)", sep = "")), limits = c(0.006,1),
                     breaks = c(1E-2,1E-1,1E0), labels = c(expression(10^-2),expression(10^-1),expression(10^0)),
                     sec.axis = sec_axis(~.*100, name= "Normalized mean center of gravity (%)", breaks = c(1,10,30,100)))+
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
  geom_point(data = dat_comp, aes (x = full_m, y = mean_wing_CGy), col = "black") +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_wing_CGy_specific), col = "gray") +
  # add iso lines
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^(1/3)), col = "black", linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean center of gravity (m)", sep = "")), limits = c(0.006,1),
                     breaks = c(1E-2,1E-1,1E0), labels = c(expression(10^-2),expression(10^-1),expression(10^0)),
                     sec.axis = sec_axis(~.*100, name= "Normalized mean center of gravity (%)", breaks = c(1,10,30,100)))+
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

rightcol <- plot_grid(blank_plot,effect_CGx,effect_CGz,CG_isometry_fullbird,effect_CGy,CG_isometry_singlewing,blank_plot,
                    #arrangement data
                    ncol = 1,
                    rel_heights = c(0.5,0.6,0.6,1,0.6,1,0.5),
                    #labels
                    labels = c("",""),
                    label_size = 10,
                    label_fontfamily = "sans")

#exported as 12x12
Figure2_final <- plot_grid(leftcol, rightcol,
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
Ixx_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Ixx_specific, xmax = max_Ixx_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Ixx_specific = dat_comp$max_Ixx_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Ixx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(min_Ixx_specific = dat_comp$min_Ixx_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Ixx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(mean_Ixx_specific = dat_comp$mean_Ixx_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Ixx_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 0.9) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  scale_x_continuous(name = expression(paste(bar("I"["xx"]))), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Iyy_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Iyy_specific, xmax = max_Iyy_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Iyy_specific = dat_comp$max_Iyy_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Iyy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(min_Iyy_specific = dat_comp$min_Iyy_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Iyy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(mean_Iyy_specific = dat_comp$mean_Iyy_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Iyy_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 0.9) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  scale_x_continuous(name = expression(paste(bar("I"["yy"]))), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Izz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Izz_specific, xmax = max_Izz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Izz_specific = dat_comp$max_Izz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Izz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(min_Izz_specific = dat_comp$min_Izz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Izz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(mean_Izz_specific = dat_comp$mean_Izz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Izz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 0.9) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  scale_x_continuous(name = expression(paste(bar("I"["zz"]))), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Ixz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "gray")+
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Ixz_specific, xmax = max_Ixz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Ixz_specific = dat_comp$max_Ixz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Ixz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(min_Ixz_specific = dat_comp$min_Ixz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Ixz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85, size = 0.8) +
  geom_point(data = aggregate(list(mean_Ixz_specific = dat_comp$mean_Ixz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Ixz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 0.9) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  scale_x_continuous(name = expression(paste(bar("I"["xz"]))), limits= c(-0.015,0.015), breaks = c(-0.015,0,0.015), labels = c(-1.5,0,1.5)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.015, xend = 0.015, y = log(0), yend = log(0))

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
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

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
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_Izz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Izz_elb_etap = dat_model_out$Izz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Izz_man_etap = dat_model_out$Izz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Izz_elbman_etap = dat_model_out$Izz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

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
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = lab_eta) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

I_isometry <- ggplot() +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^Ixx_model_mcmc_output$solutions[2,2], ymax = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^Ixx_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.15) +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Iyy_model_mcmc_output$solutions[2,2], ymax = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^Iyy_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray40", alpha = 0.15) +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Izz_model_mcmc_output$solutions[2,2], ymax = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Izz_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray60", alpha = 0.15) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^Ixx_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^Iyy_model_mcmc_output$solutions[2,1]),
            col = "gray40") +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^Izz_model_mcmc_output$solutions[2,1]),
            col = "gray60") +
  # geom_line(data = dat_comp, aes(x = full_m, y = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "black", linetype = 2) +
  # geom_line(data = dat_comp, aes(x = full_m, y = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "gray40", linetype = 2) +
  # geom_line(data = dat_comp, aes(x = full_m, y = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "gray60", linetype = 2) +
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
  #add data
  geom_point(data = dat_comp, aes(x = full_m, y = max_q, col = species_order), alpha = 0.6, pch = 19) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q, col = species_order), pch = 1) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10), name = "Mass (kg",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', limits = c(0.1,10), breaks = c(0.1,1,10), labels = c(expression(10^-1),expression(10^0),expression(10^1)), name = expression(paste("Longitudinal agility metric (s"^{-2},")")))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-1, yend = 1E1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

# -------- Combine panels into figure ------------
#exported as 6x14
midrow <- plot_grid(effect_Ixx,effect_Iyy, effect_Izz,effect_Ixz,
                       #arrangement data
                       nrow = 4,
                       rel_widths = c(1,1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")
toprow <- plot_grid(phylo_plot_complete,Ixx_specific,Iyy_specific, Izz_specific,Ixz_specific,
                       #arrangement data
                       ncol = 5,
                       rel_widths = c(0.5,1,1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")
bottomrow <- plot_grid(midrow,I_isometry,q_plot,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")
#exported as 7x8
figure_final <- plot_grid(toprow,bottomrow,
                           #arrangement data
                           ncol = 1, nrow = 2, rel_heights = c(1,0.7))

## ------------------------------------------------------------------------------------
## --------------------------------- Figure 4 -----------------------------------------
## ------------------------------------------------------------------------------------

var_test_plot <- ggplot()+
  geom_point(data = dat_var_tot, aes(x = log(min_val), y = log(mean_sig_sq), col = mean_sig_sq)) +
  scale_colour_gradientn(colors = cc_var, name = "rate of variance evolution") +
  th
# Extract the correct colours
ggplot_build(var_test_plot)$data
