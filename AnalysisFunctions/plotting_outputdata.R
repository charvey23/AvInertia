source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/process_outputdata.R")
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/plotting_setup.R")
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/calc_evoparams.R")
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/analyse_CGsens.R")
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/calc_CGparams.R")
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/calc_MOIparams.R")
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/calc_manparams.R")


## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 1 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------
#----------- Panel h -------------
validation_Ixx_plot <- ggplot()+#exported as 4x5
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = 3.76*10^-3*full_m^1.89, ymax = 3.76*10^-3*full_m^2.22),
              col = NA, fill = "gray70", alpha = 0.15) + # Kirkpatrick, 1994, 99% fiducity limits
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = 1.94*10^-3*full_m^1.787, ymax = 1.94*10^-3*full_m^2.134),
              col = NA, fill = "gray50", alpha = 0.2) + # Berg and Rayner, 1995, 95% fiducity limits
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(Ixx_val_mcmc_output$solutions[1,1])*full_m^Ixx_val_mcmc_output$solutions[2,2], ymax = exp(Ixx_val_mcmc_output$solutions[1,1])*full_m^Ixx_val_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.2) + # AvInertia
  #geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_wing_Ixx, ymin = sachs_pred_Ixx, col = species_order), alpha = 0.5) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black", linetype = 3) + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "black", linetype = 2) + # Berg and Rayner, 1995
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(Ixx_val_mcmc_output$solutions[1,1])*full_m^Ixx_val_mcmc_output$solutions[2,1]),
            col = "black") + # AvInertia
  #geom_point(data = dat_comp, aes(x = full_m, y = sachs_pred_Ixx, fill = species_order, col = species_order), pch = 22, size = 1.5) + # Sachs, 2005
  geom_point(data = Ixx_max, aes(x = full_m, y = wing_hum_Ixx, col = species_order), pch = 19, size = 2, alpha = 0.6)+
  geom_point(data = Ixx_max, aes(x = full_m, y = wing_hum_Ixx, col = species_order), pch = 1, size = 2)+
  # theme control
  th  +
  theme(legend.position="none") +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #axis control
  scale_x_continuous(trans='log10', name = "Body mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Maximum wing I"["xx"]," (kg-m"^2,")", sep = "")),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

#----------- Panel g -------------
# create the convex hulls
chulls_elbman <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                       function(df) df[chull(df$elbow, df$manus), ])

ROM_plot <- ggplot()+ #exported as 3.5x10
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
               aes(x = CGx_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(CGx_man_etap = dat_model_out$CGx_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(CGx_elbman_etap = dat_model_out$CGx_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
               aes(x = CGz_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(CGz_man_etap = dat_model_out$CGz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(CGz_elbman_etap = dat_model_out$CGz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
               aes(x = CGy_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(CGy_man_etap = dat_model_out$CGy_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(CGy_elbman_etap = dat_model_out$CGy_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.3) +
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray60", alpha = 0.3) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,1]),
            col = "gray60") +
  # add specific data
  geom_point(data = dat_comp, aes (x = full_m, y = -mean_CGx_specific_orgShoulder), col = "black", pch = 0) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_CGz_specific_orgDorsal), col = "gray60", pch = 2) +
  # add iso lines
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "black", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray60", linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Body mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean ", bar("CG"), " (% of full length)", sep = "")), limits = c(0.02,0.2),
                     breaks = c(0.01,0.02,0.1,0.2), labels = c(1,2,10,20)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 0.02, yend = 0.2) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

CG_isometry_wingonly <- ggplot() +
  # add specific line fits
  geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray40", alpha = 0.3) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,1]),
            col = "gray40") +
  # add specific data
  geom_point(data = dat_comp, aes (x = full_m, y = mean_wing_CGy_specific), col = "gray40", pch = 1) +
  # add iso lines
  geom_line(data = dat_comp, aes(x = full_m, y = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^(0)), col = "gray40", linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Body mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean ", bar("CG"), " (% of halfspan)", sep = "")), limits = c(0.02,0.2),
                     breaks = c(0.01,0.02,0.1,0.2), labels = c(1,2,10,20)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 0.02, yend = 0.2) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)


CGrange_plot <- ggplot() +
  # add model fit
  geom_ribbon(data = shoulder_motion, aes(x = span_ratio,
                                          ymin = (CG_range_model_mcmc_output$solutions[1,1] + CG_range_model_mcmc_output$solutions[2,2]*span_ratio),
                                          ymax = (CG_range_model_mcmc_output$solutions[1,1] + CG_range_model_mcmc_output$solutions[2,3]*span_ratio)),
              col = NA, fill = "gray50", alpha = 0.2) +
  geom_line(data = shoulder_motion, aes(x = span_ratio, y = (CG_range_model_mcmc_output$solutions[1,1] + CG_range_model_mcmc_output$solutions[2,1]*span_ratio)),
            col = "black") +
  #add data
  geom_point(data = shoulder_motion, aes(x = span_ratio, y = range_CGx_specific, col = species_order), alpha = 0.6, pch = 15,size = 2) +
  #geom_point(data = shoulder_motion, aes(x = span_ratio, y = range_CGx_specific, col = species_order), pch = 0,size = 0.8) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(limits = c(0,3.25), breaks = c(0,1,2,3,4), labels = c(0,100,200,300,400), name = "Maximum wingspan/body length (%)") +
  scale_y_continuous(limits = c(-0.025,0.25), breaks = c(0,0.05,0.1,0.15,0.2,0.25), labels = c(0,5,10,15,20,25), name = "Range of the CG position (% of body length)")+
  geom_rangeframe() +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = 0, yend = 0.25) +
  annotate(geom = "segment", x = 0, xend = 3, y = -Inf, yend = -Inf)


# -------- Combine into panel ------------

leftcol <- plot_grid(phylo_plot_complete,CGxz_plot, wingCGxy_plot,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1,0.6,0.6),
                       #labels
                       labels = c("","",""),
                       label_size = 10,
                       label_fontfamily = "sans")

rightcol_ind <- plot_grid(effect_CGx,effect_CGz,effect_CGy,CG_isometry_fullbird,CG_isometry_wingonly,CGrange_plot,
                      #arrangement data
                      ncol = 1,
                      rel_heights = c(0.5,0.5,0.5,1.3,1.3,1.4),
                      #labels
                      labels = c("",""),
                      label_size = 10,
                      label_fontfamily = "sans")

#exported as 12x12
Figure2_final <- plot_grid(leftcol, rightcol_ind,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(1,0.7),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")

## ------------------------------------------------------------------------------------
## --------------------------------- Figure 3 -----------------------------------------
## ------------------------------------------------------------------------------------

# Panel A

set_alp = 1
col_tail  = "gray30"
col_torso = "gray40"
col_neck  = "gray60"
col_head  = "gray80"
col_skin  = "#BEA4BD"
col_musc  = "#FAC5C6"
col_feat  = "#A0B3DC"
col_bone  = "#FCDFC7"
Ixx_cont = ggplot() +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx+neck_con_Ixx+torso_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx+neck_con_Ixx+torso_con_Ixx+tail_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_tail, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx+neck_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx+neck_con_Ixx+torso_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_torso, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx+neck_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_neck, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx+head_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_head, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx+skin_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_skin, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx+feat_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx+musc_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_musc, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Ixx,
                                                         xmax=bone_con_Ixx+feat_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_feat, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=0,
                                                         xmax=bone_con_Ixx,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_bone, alpha = set_alp) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  scale_y_continuous(breaks = seq(1,22,1), labels = phylo_order$species, trans="reverse", name = "") +
  scale_x_continuous(limits = c(0,1.01), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = Inf, yend = Inf)

Iyy_cont = ggplot() +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy+torso_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy+torso_con_Iyy+tail_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_tail, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy+torso_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_torso, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_neck, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_head, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_skin, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy+feat_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_musc, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_feat, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=0,
                                                         xmax=bone_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_bone, alpha = set_alp) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  scale_y_continuous(breaks = seq(1,22,1), labels = phylo_order$species, trans="reverse", name = "") +
  scale_x_continuous(limits = c(0,1.01), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = Inf, yend = Inf)

Izz_cont = ggplot() +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz+neck_con_Izz+torso_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz+neck_con_Izz+torso_con_Izz+tail_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_tail, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz+neck_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz+neck_con_Izz+torso_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_torso, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz+neck_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_neck, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz+head_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_head, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz+feat_con_Izz+musc_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz+musc_con_Izz+skin_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_skin, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz+feat_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz+musc_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_musc, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=bone_con_Izz,
                                                         xmax=bone_con_Izz+feat_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_feat, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="max"), aes(xmin=0,
                                                         xmax=bone_con_Izz,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_bone, alpha = set_alp) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  scale_y_continuous(breaks = seq(1,22,1), labels = phylo_order$species, trans="reverse", name = "") +
  scale_x_continuous(limits = c(0,1.01), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = Inf, yend = Inf)


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
  scale_x_continuous(name = "Maximum I/Minimum I", limits= c(1,11.8), breaks = c(1,3,5,7,9,11)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 1, xend = 11, y = log(0), yend = log(0))

effect_Ixx <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Ixx_elb_etap = dat_model_out$Ixx_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(Ixx_man_etap = dat_model_out$Ixx_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(Ixx_elbman_etap = dat_model_out$Ixx_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
               aes(x = Iyy_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(Iyy_man_etap = dat_model_out$Iyy_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(Iyy_elbman_etap = dat_model_out$Iyy_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
               aes(x = Izz_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(Izz_man_etap = dat_model_out$Izz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(Izz_elbman_etap = dat_model_out$Izz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
               aes(x = Ixz_elb_etap), fill = col_elb, alpha = 0.6)  +
  geom_density(data = aggregate(list(Ixz_man_etap = dat_model_out$Ixz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_man_etap), fill = col_man, alpha = 0.6)  +
  geom_density(data = aggregate(list(Ixz_elbman_etap = dat_model_out$Ixz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_elbman_etap), fill = col_elbman, alpha = 0.6)  +
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
  #geom_line(data = dat_comp, aes(x = full_m, y = exp(Ixx_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "black", linetype = 2) +
  #geom_line(data = dat_comp, aes(x = full_m, y = exp(Iyy_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "gray40", linetype = 2) +
  #geom_line(data = dat_comp, aes(x = full_m, y = exp(Izz_model_mcmc_output$solutions[1,1])*full_m^(5/3)), col = "gray60", linetype = 2) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_Ixx), col = "black", pch = 0, size = 0.8) +
  geom_point(data = dat_comp, aes (x = full_m, y = mean_Iyy), col = "gray40", pch = 1, size = 0.8)+
  geom_point(data = dat_comp, aes (x = full_m, y = mean_Izz), col = "gray60", pch = 2, alpha = 0.8, size = 0.8)+
  th +
  scale_x_continuous(trans='log10', name = "Body mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Mean moment of inertia (kg-m"^2,")", sep = "")), limits = c(1E-6,2E-1),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)


# -------- Combine panels into figure ------------
bottomrow <- plot_grid(I_ratio,Ixx_cont,Iyy_cont,Izz_cont,
                           #arrangement data
                           nrow = 1,
                           rel_heights = c(1,1,1,1),
                           #labels
                           labels = c("c","d","e","f"),
                           label_size = 10,
                           label_fontfamily = "sans")
toprow <- plot_grid(I_isometry,blank_plot,blank_plot,blank_plot,
                     #arrangement data
                     nrow = 1,
                     rel_widths = c(1.8,1,1,1),
                     #labels
                     labels = c("a","b"),
                     label_size = 10,
                     label_fontfamily = "sans")
Figure3_final <- plot_grid(toprow,bottomrow, # exported as 6x6.8
                           #arrangement data
                           ncol = 1,
                           rel_heights = c(0.7,1),
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
q_nd_plot <- ggplot() +
  # geom_ribbon(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, ymin = exp(max_q_nd_model_mcmc_output$solutions[1,1])*full_m^max_q_nd_model_mcmc_output$solutions[2,2], ymax = exp(max_q_nd_model_mcmc_output$solutions[1,1])*full_m^max_q_nd_model_mcmc_output$solutions[2,3]),
  #             col = NA, fill = "gray60", alpha = 0.15) +
  # # add model fit
  # geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(max_q_nd_model_mcmc_output$solutions[1,1])*full_m^max_q_nd_model_mcmc_output$solutions[2,1]),
  #           col = "gray60") +
  geom_rect(aes(ymin = 0, ymax = 1.5, xmin = 0, xmax = Inf), alpha = 0.1) +
  #add data
  geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_q_nd, ymin = min_q_nd, col = species_order), alpha = 0.5) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q_nd, col = species_order), alpha = 0.6, pch = 15,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_q_nd, col = species_order), pch = 0,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_q_nd, col = species_order), alpha = 0.6, pch = 15,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_q_nd, col = species_order), pch = 0,size = 0.8) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10), name = "Body mass (kg)",
                     labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(limits = c(-1.5,1.5), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), name = expression(paste("Normalized pitch agility")))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = -1.5, yend = 1.5) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = -Inf, yend = -Inf)

sm_nd_plot <- ggplot() +
  geom_rect(aes(ymin = -0.6, ymax = 0, xmin = 0, xmax = Inf), alpha = 0.1) +
  #add data
  annotate(geom="text", x=5, y=0.4, label="max static margin \n optimal phenotype",color="black") +
  annotate(geom="text", x=5, y=-0.29, label="min static margin \n optimal phenotype",color="black") +
  annotate("segment", x = 5, xend = 5, y = 0.34, yend = 0.28, colour = "black", size=0.4,
           arrow=arrow(length = unit(0.15, 'cm'), type = 'closed')) +
  annotate("segment", x = 5, xend = 5, y = -0.23, yend = -0.17, colour = "black", size=0.4,
           arrow=arrow(length = unit(0.15, 'cm'), type = 'closed')) +
  geom_hline(yintercept = OU_maxsm$opt$z0, linetype = "dashed") +

  geom_hline(yintercept = OU_minsm$opt$z0, linetype = "dashed") +
  #geom_hline(yintercept = maxstab_tail_OU_outputs$opt$z0, linetype = "dashed", col = "gray70") +
  #geom_hline(yintercept = minstab_tail_OU_outputs$opt$z0, linetype = "dashed", col = "gray70") +
  geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_sm_nd, ymin = min_sm_nd, col = species_order), alpha = 0.5) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_sm_nd, col = species_order), alpha = 0.6, pch = 15,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = max_sm_nd, col = species_order), pch = 0,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_sm_nd, col = species_order), alpha = 0.6, pch = 15,size = 0.8) +
  geom_point(data = dat_comp, aes(x = full_m, y = min_sm_nd, col = species_order), pch = 0,size = 0.8) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10), name = "Body mass (kg)",
                     labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(limits = c(-0.6,0.6), name = expression(paste("Normalized static margin (% of c"["rmax"],")")))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = -0.6, yend = 0.6) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = -Inf, yend = -Inf)


#exported as 3x6.5
plot_agility <- plot_grid(q_nd_plot,sm_nd_plot,
                       #arrangement data
                       nrow = 1,
                       rel_heights = c(1,1),
                       #labels
                       labels = c("a","b"),
                       label_size = 10,
                       label_fontfamily = "sans")

#exported as 11x4
chord_length_plot <- ggplot()+
  geom_point(data = stab_check, aes(x = 1, y = 1, col = max_sm_nd, size = abs(max_q_nd)), alpha = 0.8) +
  geom_point(data = stab_check, aes(x = 1, y = 1, col = max_sm_nd, size = abs(max_q_nd)), pch = 1, alpha = 1) +
  geom_point(data = stab_check, aes(x = 0.98, y = 1, col = min_sm_nd, size = abs(min_q_nd)), alpha = 0.8) +
  geom_point(data = stab_check, aes(x = 0.98, y = 1, col = min_sm_nd, size = abs(min_q_nd)), pch = 1, alpha = 1) +
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
  #scale_size(range = c(0.05, 10), name="Normalized pitch agility", breaks = c(0.3,0.6,0.9,1.2,1.5), labels = c(0.2,0.4,0.6,0.8,1)) + #scales based on area
  scale_colour_gradient2(midpoint=0, low="#1D0747", mid="#FFF9D3", high="#00510A", space ="Lab", name = "static margin (% of max. root chord)") +
  scale_y_continuous(limits = c(0.955,1.025)) +
  scale_x_continuous(limits = c(0.955,1.025))


Iyy_cont = ggplot() +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy+torso_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy+torso_con_Iyy+tail_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_tail, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy+torso_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_torso, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy+neck_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_neck, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy+head_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_head, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy+skin_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_skin, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy+feat_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy+musc_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_musc, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=bone_con_Iyy,
                                                         xmax=bone_con_Iyy+feat_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_feat, alpha = set_alp) +
  geom_rect(data=subset(I_contr,type_metric=="min"), aes(xmin=0,
                                                         xmax=bone_con_Iyy,
                                                         ymin=as.numeric(as.factor(species_order))-0.5,
                                                         ymax=as.numeric(as.factor(species_order))+0.5), stat = "identity", fill = col_bone, alpha = set_alp) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  scale_y_continuous(breaks = seq(1,22,1), labels = phylo_order$species, trans="reverse", name = "") +
  scale_x_continuous(limits = c(0,1.01), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = Inf, yend = Inf)

### ------------- Extended Data Figure ----------

test_max = aggregate(list(prop_q_dot = dat_final$prop_q_dot,
                          full_m = dat_final$full_m,
                          sm_nd = dat_final$sm_nd),
                     by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), max)
test_min = aggregate(list(prop_q_dot = dat_final$prop_q_dot,
                          full_m = dat_final$full_m,
                          sm_nd = dat_final$sm_nd),
                     by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), min)

chulls_xz <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                   function(df) df[chull(df$elbow, df$manus), ])

ROM_plot_sm <- ggplot()+
  #geom_polygon(data = chulls_xz, aes(x = elbow, y = manus), col = "black", fill = NA, alpha = 0.5) +
  geom_point(data = dat_final, aes(x =elbow, y = manus, col = sm_nd)) +
  geom_point(data = dat_final[which(dat_final$species %in% test_max$species & dat_final$sm_nd %in% test_max$sm_nd),],
             aes(x = elbow, y = manus, fill = sm_nd), col = "black", alpha = 1, pch = 23) +
  geom_point(data = dat_final[which(dat_final$species %in% test_min$species & dat_final$sm_nd %in% test_min$sm_nd),],
             aes(x = elbow, y = manus, fill = sm_nd), col = "black", alpha = 1, pch = 23) +
  facet_wrap(~species_order, nrow = 2) +
  #colour control
  scale_colour_gradient2(midpoint=0, low="#1D0747", mid="#FFF9D3", high="#00510A", space ="Lab", name = "static margin (% of max. root chord)") +
  scale_fill_gradient2(midpoint=0, low="#1D0747", mid="#FFF9D3", high="#00510A", space ="Lab", name = "static margin (% of max. root chord)") +
  #theme control
  th +
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA))+
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


### ------------- Plot the tail volume coefficients for the supplemental materials ------------
tail_plot <- ggplot() +
  geom_point(data = aggregate(list(Vh = dat_final$Vh, full_m = dat_final$full_m),
                              by=list(species_order = dat_final$species_order, BirdID = dat_final$BirdID), max),
             aes(x = full_m, y = Vh, col = species_order)) +
  geom_point(data = aggregate(list(Vh = dat_final$Vh, full_m = dat_final$full_m),
                              by=list(species_order = dat_final$species_order, BirdID = dat_final$BirdID), min),
             aes(x = full_m, y = Vh, col = species_order)) +
  geom_segment(data = aggregate(list(Vh = dat_final$Vh, full_m = dat_final$full_m, Vh_max = -dat_final$Vh),
                              by=list(species_order = dat_final$species_order, BirdID = dat_final$BirdID), min),
             aes(x = full_m, xend = full_m, y = Vh, yend = -Vh_max, col = species_order)) +
  geom_hline(yintercept = 0.4) +
  annotate(geom="text", x=0.05, y=2, label="Approximate lower bound for \n traditional aircraft tail volume coefficients",
           color="black") +
  annotate("segment", x = 0.05, xend = 0.07, y = 1.3, yend = 0.45, colour = "black", size=0.4,
           arrow=arrow(length = unit(0.25, 'cm'), type = 'closed')) +
  #colour control
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  th +
  theme(legend.position="none") +
  # axis control
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10), name = "Body mass (kg)",
                     labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', limits = c(10^-3,10), breaks = c(0.001,0.01,0.1,1,10),
                     labels = c(expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)), name = expression(paste("Tail volume coefficient")))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-3, yend = 10) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)
