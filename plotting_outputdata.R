library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe
library(gridExtra) # for using grid arrange
library(ggtree)    # for plotting phylogeny
library(plyr)      # for computing convex hulls

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

#--- Blank plot to fill space ---
blank_plot <- ggplot() + theme_void()

## ----------------- Plot of linear model outputs ------------------
# order species to match the phylogeny

man_fixed = c(80,100,120)*0.001
elb_fixed = c(80,100,120)*0.001


# ---------- Pre-define colours ----------
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
col_xz = "#af0069"

col_elb = "#55b3b1"
col_man = "#af0069"
col_elbman = "#5F0995"

# Species rainbow colours
cc_rain1  <- scales::seq_gradient_pal("#DA9101", "#CA302F", "Lab")(seq(0,1,length.out=6))
cc_rain2  <- scales::seq_gradient_pal("#FCC201", "#DA9101", "Lab")(seq(0,1,length.out=3))
cc_rain3  <- scales::seq_gradient_pal("#138715", "#FCC201", "Lab")(seq(0,1,length.out=6))
cc_rain4  <- scales::seq_gradient_pal("#1FC3CD", "#138715", "Lab")(seq(0,1,length.out=4))
cc_rain5  <- scales::seq_gradient_pal("#304CC8", "#1FC3CD", "Lab")(seq(0,1,length.out=5))
cc_rain6  <- scales::seq_gradient_pal("#5F2CC8", "#304CC8", "Lab")(seq(0,1,length.out=3))
cc_rain   <- c(cc_rain6, cc_rain5[2:5], cc_rain4[2:4], cc_rain3[2:6],cc_rain2[2:3],cc_rain1[2:6])

# -------------- Set up common labels -----------
lab_eta = expression(paste(eta^2))


# --------------- Set up Phylogeny info -----------
phylo_label = c("Cooper's hawk","Sharp-shinned hawk", "Western grebe", "Mallard","Great blue heron","Canada goose","Common nighthawk",
                "Lady Amherst's pheasant","Northern flicker","Pigeon","Common raven","Steller's jay","Black swift",
                "Merlin","Peregrine falcon","Leach's storm petrel","Glaucous-winged gull","Himalayan monal","Silver pheasant",
                "Belted kingfisher","American white pelican","Barn owl")
dat_phylo_plot <- data.frame(label = pruned_mcc$tip.label, genus = phylo_label)
phylo_plot <- ggtree(pruned_mcc) %<+% dat_phylo_plot %>% rotate(23) +
  geom_tiplab(aes(label=genus)) + xlim(0,150) + theme(plot.margin = unit(c(14,8,14,8), "mm"))

phylo_data <- ggtree::ggtree(pruned_mcc)[["data"]][,c("label","y")]
phylo_order <- merge(phylo_data[complete.cases(phylo_data),], unique(dat_bird[,c("binomial","species")]), by.x = "label", by.y = "binomial")
phylo_order <- arrange(phylo_order,y)

dat_comp$species_order = factor(dat_comp$species, levels = phylo_order$species)
dat_final$species_order = factor(dat_final$species, levels = phylo_order$species)
shoulder_motion$species_order = factor(shoulder_motion$species, levels = phylo_order$species)

dat_model_out$species <- factor(dat_model_out$species, levels = phylo_order$species)
shading <- data.frame(col1 = levels(dat_model_out$species)[seq(from = 1, to = max(as.numeric(as.factor(dat_model_out$species)))-1, by = 2)],
                      col2 = levels(dat_model_out$species)[seq(from = 2, to = max(as.numeric(as.factor(dat_model_out$species))), by = 2)])

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

chulls_elbman <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                   function(df) df[chull(df$elbow, df$manus), ])

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
# Panel A adjust data
chulls_CGxz <- ddply(dat_final[,c("species_order","BirdID","full_CGx_specific_orgShoulder","full_CGz_specific_orgShoulder")], .(species_order, BirdID),
                   function(df) df[chull(df$full_CGx_specific_orgShoulder, df$full_CGz_specific_orgShoulder), ])

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

CG_extremes = aggregate(list(mean_CGz     = dat_final$full_CGz_specific_orgShoulder,
                             mean_CGx     = dat_final$full_CGx_specific_orgShoulder),  by=list(species_order = dat_final$species_order), mean)
tmp = aggregate(list(mean_CGxback  = shoulder_motion$backwards_CGx_specific,
                     mean_CGzup = shoulder_motion$upwards_CGz_specific),  by=list(species_order = shoulder_motion$species_order), min)
CG_extremes = merge(CG_extremes, tmp, by = "species_order")
tmp = aggregate(list(mean_CGxfore  = shoulder_motion$forward_CGx_specific,
                     mean_CGzdn = shoulder_motion$downwards_CGz_specific),  by=list(species_order = shoulder_motion$species_order), max)
CG_extremes = merge(CG_extremes, tmp, by = "species_order")

cg_x_shape1 = melt(CG_extremes[,c("species_order","mean_CGxfore","mean_CGx","mean_CGxback")], id= "species_order")
cg_x_shape2 = melt(CG_extremes[,c("species_order","mean_CGx")], id= "species_order")
cg_x_shape = rbind(cg_x_shape1,cg_x_shape2)
cg_z_shape1 = melt(CG_extremes[,c("species_order","mean_CGz","mean_CGzup")], id= "species_order")
cg_z_shape2 = melt(CG_extremes[,c("species_order","mean_CGz","mean_CGzdn")], id= "species_order")
cg_z_shape = rbind(cg_z_shape1,cg_z_shape2)

# need to adjust species order so that pigeon comes last

# Panel A plot
CGxz_plot <- ggplot()+
  # add data
  geom_rect(data = chulls_CGxz, aes(xmin =0, ymin = -0.1, xmax = 0.2, ymax = 0, group = species_order), linetype="dashed", color = "gray", fill = NA, alpha = 0.25) +
  # geom_path(data = x_shape, aes(x = -value, y = -z_shape$value), col = "black", alpha = 0.2) +
  geom_polygon(data = cg_x_shape, aes(x = -value, y = -cg_z_shape$value, fill = species_order), col = NA, alpha = 0.5) +
  geom_polygon(data = chulls_CGxz, aes(x = -full_CGx_specific_orgShoulder, y = -full_CGz_specific_orgShoulder, col = species_order, fill = species_order, group = BirdID)) +
  facet_wrap(~species_order, nrow = 3) +
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
# adjust data
max_wing$species_order = factor(max_wing$species, levels = phylo_order$species)
tmpx = melt(max_wing[,c("species_order","pt12_X","pt6_X","pt7_X","pt8_X","pt9_X","pt10_X","pt11_X","edge_X")], id = "species_order")
tmpy = melt(max_wing[,c("species_order","pt12_Y","pt6_Y","pt7_Y","pt8_Y","pt9_Y","pt10_Y","pt11_Y","edge_Y")], id = "species_order")
wing_shape = cbind(tmpx,tmpy[,c(2,3)])
colnames(wing_shape) <- c("species_order","variable_x","value_x","variable_y","value_y")
chulls_CGxy <- ddply(dat_final[,c("species_order","BirdID","wing_CGx_specific_orgShoulder","wing_CGy_specific_orgShoulder")], .(species_order, BirdID),
                     function(df) df[chull(df$wing_CGx_specific_orgShoulder, df$wing_CGy_specific_orgShoulder), ])
# plot
wingCGxy_plot <- ggplot() +
  geom_rect(data = chulls_CGxy, aes(xmin =0, ymin = -0.14, xmax = 0.3, ymax = 0, group = species_order), linetype="dashed", color = "gray", fill = NA, alpha = 0.25) +
  #geom_path(data = wing_shape, aes(x = value_y, y = value_x), col = "black", alpha = 0.2) +
  geom_polygon(data = chulls_CGxy, aes(x = wing_CGy_specific_orgShoulder, y = wing_CGx_specific_orgShoulder, col = species_order, fill = species_order, group = BirdID)) +
  facet_wrap(~species_order, nrow = 3) +
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




### ------------ Fit PGLSS model to the positions -----------
CGx_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(-mean_CGx_orgBeak) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = morpho_data_means,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
CGx_model_mcmc_output  = summary(CGx_model_mcmc)

CGx_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(-mean_CGx_specific) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = morpho_data_means,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
CGx_sp_model_mcmc_output  = summary(CGx_sp_model_mcmc)

CGy_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_wing_CGy) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = morpho_data_means,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
CGy_model_mcmc_output  = summary(CGy_model_mcmc)

CGy_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_wing_CGy_specific) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = morpho_data_means,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
CGy_sp_model_mcmc_output  = summary(CGy_sp_model_mcmc)

CGz_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_CGz_orgDorsal) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = morpho_data_means,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
CGz_model_mcmc_output  = summary(CGz_model_mcmc)

CGz_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_CGz_specific) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = morpho_data_means,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
CGz_sp_model_mcmc_output  = summary(CGz_sp_model_mcmc)


CG_isometry_fullbird <- ggplot() +
  # add specific line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,1]),
            col = "gray") +
  geom_ribbon(data = data.fit, aes(x = full_m, ymin = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGx_sp_model_mcmc_output$solutions[1,1])*full_m^CGx_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray", alpha = 0.3) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,1]),
            col = "gray") +
  geom_ribbon(data = data.fit, aes(x = full_m, ymin = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGz_sp_model_mcmc_output$solutions[1,1])*full_m^CGz_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray", alpha = 0.3) +
  # add line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^CGx_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_ribbon(data = data.fit, aes(x = full_m, ymin = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^CGx_model_mcmc_output$solutions[2,2], ymax = exp(CGx_model_mcmc_output$solutions[1,1])*full_m^CGx_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "black", alpha = 0.3) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^CGz_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_ribbon(data = data.fit, aes(x = full_m, ymin = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^CGz_model_mcmc_output$solutions[2,2], ymax = exp(CGz_model_mcmc_output$solutions[1,1])*full_m^CGz_model_mcmc_output$solutions[2,3]),
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
  geom_ribbon(data = data.fit, aes(x = full_m, ymin = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,2], ymax = exp(CGy_sp_model_mcmc_output$solutions[1,1])*full_m^CGy_sp_model_mcmc_output$solutions[2,3]),
              col = NA, fill = "gray", alpha = 0.3) +
  # add line fits
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^CGy_model_mcmc_output$solutions[2,1]),
            col = "black") +
  geom_ribbon(data = data.fit, aes(x = full_m, ymin = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^CGy_model_mcmc_output$solutions[2,2], ymax = exp(CGy_model_mcmc_output$solutions[1,1])*full_m^CGy_model_mcmc_output$solutions[2,3]),
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

effectsize_full <- plot_grid(effect_CGx,effect_CGz,
                                 #arrangement data
                                 nrow = 2,
                                 rel_heights = c(1,1),
                                 #labels
                                 labels = c("",""),
                                 label_size = 10,
                                 label_fontfamily = "sans")
toprow <- plot_grid(blank_plot, effectsize_full,
                    #arrangement data
                    ncol = 2,
                    rel_widths = c(1,0.4),
                    #labels
                    labels = c("",""),
                    label_size = 10,
                    label_fontfamily = "sans")

fullbird_panel <- plot_grid(CGxz_plot, CG_isometry_fullbird,
                              #arrangement data
                              ncol = 2,
                              rel_widths = c(1,0.4),
                              #labels
                              labels = c("","",""),
                              label_size = 10,
                              label_fontfamily = "sans")

effectsize_single <- plot_grid(blank_plot,effect_CGy,
                             #arrangement data
                             nrow = 2,
                             rel_heights = c(1,1),
                             #labels
                             labels = c("",""),
                             label_size = 10,
                             label_fontfamily = "sans")

midrow <- plot_grid(blank_plot, effectsize_single,
                    #arrangement data
                    ncol = 2,
                    rel_widths = c(1,0.4),
                    #labels
                    labels = c("",""),
                    label_size = 10,
                    label_fontfamily = "sans")

bottomrow <- plot_grid(wingCGxy_plot, CG_isometry_singlewing,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(1,0.4),
                       #labels
                       labels = c("","",""),
                       label_size = 10,
                       label_fontfamily = "sans")
#exported as 12x12
Figure2_final <- plot_grid(toprow,fullbird_panel,midrow,bottomrow,
                       #arrangement data
                       nrow = 4,
                       rel_heights = c(0.8,1,0.8,1),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")






chulls_Ixxzz <- ddply(dat_final[,c("species_order","BirdID","full_Ixx_specific","full_Izz_specific")], .(species_order, BirdID),
                     function(df) df[chull(df$full_Ixx_specific, df$full_Izz_specific), ])

# Panel C plot
wingIxxIzz_plot <- ggplot()+
  geom_polygon(data = chulls_Ixxzz, aes(x = full_Ixx_specific, y = full_Izz_specific, col = species_order, fill = species_order, group = BirdID)) +
  facet_wrap(~species_order, ncol = 5) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  #theme control
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        # Background of panel
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.line =  element_line(colour = "black", linetype=1))+
  theme(legend.position="none") +
  coord_fixed() +
  scale_x_continuous(name='y (%)', limits = c(-0.01,1.3)) +
  scale_y_continuous(name='x (%)', limits = c(-0.8,0.2))




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
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
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
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
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
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
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
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

# -------- Combine panels into figure ------------
bottomrow <- plot_grid(phylo_plot,CGx_specific,CGy_specific, CGz_specific,
                       #arrangement data
                       ncol = 4,
                       rel_widths = c(1.5,1,1,1),
                       #labels
                       labels = c("A","B","C","D"),
                       label_size = 10,
                       label_fontfamily = "sans")

#exported as 6x12
figure_final <- plot_grid(toprow,bottomrow,
                           #arrangement data
                           ncol = 1, nrow = 2, rel_heights = c(1,1))
#expot as 6x10
bottomrow <- plot_grid(phylo_plot,CGx_angleeffect,CGy_angleeffect,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1.5,1,1),
                       #labels
                       labels = c("A","B","C"),
                       label_size = 10,
                       label_fontfamily = "sans")

bottomrow <- plot_grid(cohens_effect_x,cohens_effect_y,cohens_effect_z,
                       #arrangement data
                       nrow = 3,
                       rel_heights = c(1,1,1),
                       #labels
                       labels = c("A","B","C"),
                       label_size = 10,
                       label_fontfamily = "sans")

effectsizes <- plot_grid(effect_CGx,effect_CGy,effect_CGz, effect_Ixx,effect_Iyy,effect_Izz,effect_Ixz,
                       #arrangement data
                       nrow = 7,
                       #labels
                       labels = c("CGx","CGy","CGz","Ixx","Iyy","Izz","Ixz"),
                       label_size = 10,
                       label_fontfamily = "sans")


## ----------------- Plot of linear model outputs ------------------
# order species to match the phylogeny
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
col_xz = "#af0069"

min_I_lb   <- aggregate(list(min_Ixx_sp_int_lb = dat_model_out$Ixx_sp_int_lb,
                             min_Iyy_sp_int_lb = dat_model_out$Iyy_sp_int_lb,
                             min_Izz_sp_int_lb = dat_model_out$Izz_sp_int_lb,
                             min_Ixz_sp_int_lb = dat_model_out$Ixz_sp_int_lb,
                             min_Ixx_sp_elb_lb1 = dat_model_out$Ixx_sp_elb_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[1],
                             min_Ixx_sp_elb_lb2 = dat_model_out$Ixx_sp_elb_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[3],
                             min_Iyy_sp_elb_lb1 = dat_model_out$Iyy_sp_elb_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[1],
                             min_Iyy_sp_elb_lb2 = dat_model_out$Iyy_sp_elb_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[3],
                             min_Izz_sp_elb_lb1 = dat_model_out$Izz_sp_elb_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[1],
                             min_Izz_sp_elb_lb2 = dat_model_out$Izz_sp_elb_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[3],
                             min_Ixz_sp_elb_lb1 = dat_model_out$Ixz_sp_elb_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[1],
                             min_Ixz_sp_elb_lb2 = dat_model_out$Ixz_sp_elb_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[3],
                             min_Ixx_sp_man_lb1 = dat_model_out$Ixx_sp_man_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[1],
                             min_Ixx_sp_man_lb2 = dat_model_out$Ixx_sp_man_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[3],
                             min_Iyy_sp_man_lb1 = dat_model_out$Iyy_sp_man_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[1],
                             min_Iyy_sp_man_lb2 = dat_model_out$Iyy_sp_man_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[3],
                             min_Izz_sp_man_lb1 = dat_model_out$Izz_sp_man_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[1],
                             min_Izz_sp_man_lb2 = dat_model_out$Izz_sp_man_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[3],
                             min_Ixz_sp_man_lb1 = dat_model_out$Ixz_sp_man_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[1],
                             min_Ixz_sp_man_lb2 = dat_model_out$Ixz_sp_man_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[3]),
                        by=list(species = dat_model_out$species), min)
max_I_ub   <- aggregate(list(max_Ixx_sp_int_ub = dat_model_out$Ixx_sp_int_ub,
                             max_Iyy_sp_int_ub = dat_model_out$Iyy_sp_int_ub,
                             max_Izz_sp_int_ub = dat_model_out$Izz_sp_int_ub,
                             max_Ixz_sp_int_ub = dat_model_out$Ixz_sp_int_ub,
                             max_Ixx_sp_elb_ub1 = dat_model_out$Ixx_sp_elb_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[1],
                             max_Ixx_sp_elb_ub2 = dat_model_out$Ixx_sp_elb_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[3],
                             max_Iyy_sp_elb_ub1 = dat_model_out$Iyy_sp_elb_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[1],
                             max_Iyy_sp_elb_ub2 = dat_model_out$Iyy_sp_elb_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[3],
                             max_Izz_sp_elb_ub1 = dat_model_out$Izz_sp_elb_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[1],
                             max_Izz_sp_elb_ub2 = dat_model_out$Izz_sp_elb_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[3],
                             max_Ixz_sp_elb_ub1 = dat_model_out$Ixz_sp_elb_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[1],
                             max_Ixz_sp_elb_ub2 = dat_model_out$Ixz_sp_elb_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[3],
                             max_Ixx_sp_man_ub1 = dat_model_out$Ixx_sp_man_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[1],
                             max_Ixx_sp_man_ub2 = dat_model_out$Ixx_sp_man_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[3],
                             max_Iyy_sp_man_ub1 = dat_model_out$Iyy_sp_man_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[1],
                             max_Iyy_sp_man_ub2 = dat_model_out$Iyy_sp_man_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[3],
                             max_Izz_sp_man_ub1 = dat_model_out$Izz_sp_man_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[1],
                             max_Izz_sp_man_ub2 = dat_model_out$Izz_sp_man_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[3],
                             max_Ixz_sp_man_ub1 = dat_model_out$Ixz_sp_man_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[1],
                             max_Ixz_sp_man_ub2 = dat_model_out$Ixz_sp_man_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[3]),
                        by=list(species = dat_model_out$species), max)
ci_I_bounds <- merge(min_I_lb,max_I_ub, by = c("species"))

dat_I_plot     <- aggregate(list(mean_Ixx_sp_int = dat_model_out$Ixx_sp_int,
                                 mean_Iyy_sp_int = dat_model_out$Iyy_sp_int,
                                 mean_Izz_sp_int = dat_model_out$Izz_sp_int,
                                 mean_Ixx_sp_elb = dat_model_out$Ixx_sp_elb,
                                 mean_Iyy_sp_elb = dat_model_out$Iyy_sp_elb,
                                 mean_Izz_sp_elb = dat_model_out$Izz_sp_elb,
                                 mean_Ixz_sp_elb = dat_model_out$Ixz_sp_elb,
                                 mean_Ixx_sp_man = dat_model_out$Ixx_sp_man,
                                 mean_Iyy_sp_man = dat_model_out$Iyy_sp_man,
                                 mean_Izz_sp_man = dat_model_out$Izz_sp_man,
                                 mean_Ixz_sp_man = dat_model_out$Ixz_sp_man,
                                 mean_Ixx_sp_elbman = dat_model_out$Ixx_sp_elbman,
                                 mean_Iyy_sp_elbman = dat_model_out$Iyy_sp_elbman,
                                 mean_Izz_sp_elbman = dat_model_out$Izz_sp_elbman,
                                 mean_Ixz_sp_elbman = dat_model_out$Ixz_sp_elbman),  by=list(species = dat_model_out$species), mean)


dat_I_sp_plot <- merge(aggregate(list(min_Ixx_specific = dat_comp$min_Ixx_specific,
                             min_Iyy_specific = dat_comp$min_Iyy_specific,
                             min_Izz_specific = dat_comp$min_Izz_specific,
                             min_Ixz_specific = dat_comp$min_Ixz_specific),  by=list(species = dat_comp$species), min),
              aggregate(list(max_Ixx_specific = dat_comp$max_Ixx_specific,
                             max_Iyy_specific = dat_comp$max_Iyy_specific,
                             max_Izz_specific = dat_comp$max_Izz_specific,
                             max_Ixz_specific = dat_comp$max_Ixz_specific),  by=list(species = dat_comp$species), max), by = "species")

Ixx_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Ixx_specific, xmax = max_Ixx_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Ixx_specific = dat_comp$max_Ixx_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Ixx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Ixx_specific = dat_comp$min_Ixx_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Ixx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Ixx_specific = dat_comp$mean_Ixx_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Ixx_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["xx"]," (% of maximum inertia)", sep = "")), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))


Iyy_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Iyy_specific, xmax = max_Iyy_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Iyy_specific = dat_comp$max_Iyy_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Iyy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Iyy_specific = dat_comp$min_Iyy_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Iyy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Iyy_specific = dat_comp$mean_Iyy_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Iyy_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["yy"]," (% of maximum inertia)", sep = "")), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Izz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Izz_specific, xmax = max_Izz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Izz_specific = dat_comp$max_Izz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Izz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Izz_specific = dat_comp$min_Izz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Izz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Izz_specific = dat_comp$mean_Izz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Izz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["zz"]," (% of maximum inertia)", sep = "")), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Ixz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Ixz_specific, xmax = max_Ixz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Ixz_specific = dat_comp$max_Ixz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Ixz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Ixz_specific = dat_comp$min_Ixz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Ixz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Ixz_specific = dat_comp$mean_Ixz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Ixz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["xz"]," (% of maximum inertia)", sep = "")), limits= c(-0.015,0.015), breaks = c(-0.015,0,0.015), labels = c(-1.5,0,1.5)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.015, xend = 0.015, y = log(0), yend = log(0))

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


Ixx_iso_model <- lm(log(max_Ixx) ~ 1 + offset((5/3)*log(full_m)), data = dat_comp)
Iyy_iso_model <- lm(log(max_Iyy) ~ 1 + offset((5/3)*log(full_m)), data = dat_comp)
Izz_iso_model <- lm(log(max_Izz) ~ 1 + offset((5/3)*log(full_m)), data = dat_comp)

I_isometry <- ggplot() +
  geom_point(data = dat_comp, aes (x = full_m, y = max_Ixx), col = col_x) +
  geom_point(data = dat_comp, aes (x = full_m, y = max_Iyy), col = "#479030")+
  geom_point(data = dat_comp, aes (x = full_m, y = max_Izz), col = col_z)+
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(Ixx_iso_model)[1])*full_m^(5/3)), col = col_x, linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(Iyy_iso_model)[1])*full_m^(5/3)), col = "#479030", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(Izz_iso_model)[1])*full_m^(5/3)), col = col_z, linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Maximum moment of inertia (kg-m"^2,")", sep = "")),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)



# export as 5x5
CGx_iso_model <- lm(log(-min_CGx) ~ 1 + offset((1/3)*log(full_m)), data = dat_comp)
CGy_iso_model <- lm(log(max_wing_CGy) ~ 1 + offset((1/3)*log(full_m)), data = dat_comp)
CGz_iso_model <- lm(log(max_CGz) ~ 1 + offset((1/3)*log(full_m)), data = dat_comp)

CG_isometry <- ggplot() +
  geom_point(data = dat_comp, aes (x = full_m, y = -min_CGx), col = col_x)+
  geom_point(data = dat_comp, aes (x = full_m, y = max_wing_CGy), col = "#479030")+
  geom_point(data = dat_comp, aes (x = full_m, y = max_CGz), col = col_z)+
  # verify the fits match PGLSS model
  #geom_line(data = test, aes(x = full_m, y = exp(fit)), col = col_x) +
  #geom_ribbon(data = test,
  #            aes(x = full_m, ymin = exp(lwr), ymax = exp(upr)),
  #            col = col_x, fill = col_x, alpha = 0.4) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(CGx_iso_model)[1])*full_m^(1/3)), col = col_x, linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(CGy_iso_model)[1])*full_m^(1/3)), col = "#479030", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(CGz_iso_model)[1])*full_m^(1/3)), col = col_z, linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Center of gravity (m)", sep = "")),
                     breaks = c(1E-2,1E-1,1E0), labels = c(expression(10^-2),expression(10^-1),expression(10^0)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-2, yend = 1E0) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

bottomrow <- plot_grid(CG_isometry,I_isometry,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(1,1),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")

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

## ----------------- Plots for sharing with others -----------------


validation_Izz_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = wing_Izz, col = clade)) +
  geom_line(data = unique(dat_final[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_final[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Ixx_body_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = full_Ixx, col = clade)) +
  geom_point(aes(x = 0.0875, y = 928/(100^2*1000)), col = "black") + # Hedrick and Biewener
  geom_point(aes(x = 0.2893, y = 12889/(100^2*1000)), col = "black") + # Hedrick and Biewener
  th +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

q_plot <- ggplot() +
  geom_point(data = aggregate(list(prop_q_dot = dat_final$prop_q_dot, full_m = dat_final$full_m),  by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), max),
             aes(x = full_m, y = prop_q_dot, fill = species), col = "black", alpha = 0.8, pch = 21) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 0.66*full_m^-0.35),col = "gray") + # Kirkpatrick, 1994
  #geom_point(data = aggregate(list(prop_q_dot = dat_final$prop_q_dot, full_m = dat_final$full_m),  by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), min),
  #           aes(x = full_m, y = prop_q_dot, col = species)) +
  th +
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10)) +
  scale_y_continuous(trans='log10', limits = c(0.1,10), breaks = c(0.1,1,10))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-1, yend = 1E1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

test_plot <- ggplot() +
  geom_point(aes(x = dat_final$elbow[which(dat_final$species %in% test$species & dat_final$prop_q_dot %in% test$prop_q_dot)],
                 y = dat_final$manus[which(dat_final$species %in% test$species & dat_final$prop_q_dot %in% test$prop_q_dot)],
                 col = log(dat_final$prop_q_dot[which(dat_final$species %in% test$species & dat_final$prop_q_dot %in% test$prop_q_dot)]))) + th

ROM_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = elbow, y = manus, col = species_order), alpha = 0.5) +
  facet_wrap(~species_order, nrow = 2) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = ggplot2::element_rect(fill = "transparent"),
    # Background behind actual data points
    plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
    axis.line =  element_line(colour = "black", linetype=1))+
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)') +
  scale_y_continuous(name='Wrist angle (°)')

# Species rainbow colours
cc_rain1  <- scales::seq_gradient_pal("#DA9101", "#CA302F", "Lab")(seq(0,1,length.out=6))
cc_rain2  <- scales::seq_gradient_pal("#FCC201", "#DA9101", "Lab")(seq(0,1,length.out=3))
cc_rain3  <- scales::seq_gradient_pal("#138715", "#FCC201", "Lab")(seq(0,1,length.out=6))
cc_rain4  <- scales::seq_gradient_pal("#304CC8", "#138715", "Lab")(seq(0,1,length.out=6))
cc_rain5  <- scales::seq_gradient_pal("#4E2388", "#304CC8", "Lab")(seq(0,1,length.out=5))
cc_rain   <- c(cc_rain5, cc_rain4[2:6], cc_rain3[2:6],cc_rain2[2:3],cc_rain1[2:6])

chulls_xz <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                   function(df) df[chull(df$elbow, df$manus), ])

ROM_plot <- ggplot()+
  geom_polygon(data = chulls_xz, aes(x = elbow, y = manus, fill = species_order), col = NA, alpha = 0.5) +
  #stat_contour_filled(data = data.fit, aes(x = elbow, y = manus, z = prop_q, colour = ..level..),
  #             breaks = quantile(data.fit$prop_q, seq(0, 1, 0.05)), size = 0.8) +
  facet_wrap(~species_order, nrow = 2) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.line =  element_line(colour = "black", linetype=1))+
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)') +
  scale_y_continuous(name='Wrist angle (°)') + theme(legend.position="none")


test <- lm(prop_q_dot_nd ~ elbow*manus + species_order, data = dat_final)
xgrid <-  seq(floor(min(dat_final$elbow)), ceiling(max(dat_final$elbow)), 1)
ygrid <-  seq(floor(min(dat_final$manus)), ceiling(max(dat_final$manus)), 1)
zgrid <-  phylo_order$species
data.fit       <- expand.grid(elbow = xgrid, manus = ygrid, species_order = zgrid)
data.fit$prop_q  <-  predict(test, newdata = data.fit)

# Visualize the wings as required - For each calibration verify that the axis is RH
m = 1:nrow(dat_wing_curr)
max = 0.4
plot(dat_wing_curr$pt2_Y[m],dat_wing_curr$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_wing_curr$pt3_Y[m],dat_wing_curr$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_wing_curr$pt4_Y[m], dat_wing_curr$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_wing_curr$pt1_Y[m],dat_wing_curr$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_wing_curr$pt8_Y[m],dat_wing_curr$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_wing_curr$pt10_Y[m],dat_wing_curr$pt10_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

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

plot(dat_complete$pt2_Y[m],dat_complete$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_complete$pt3_Y[m],dat_complete$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_complete$pt4_Y[m], dat_complete$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_complete$pt1_Y[m],dat_complete$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_complete$pt8_Y[m],dat_complete$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_complete$pt9_Y[m],dat_complete$pt9_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

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




chulls_xz <- ddply(dat_final[,c("species","BirdID","full_CGx_specific","full_CGz_specific")], .(species, BirdID),
                   function(df) df[chull(df$full_CGx_specific, df$full_CGz_specific), ])

CG_side_view <- ggplot()+
  geom_polygon(data=chulls_xz, aes(x=-full_CGx_specific, y=-full_CGz_specific, group=interaction(species,BirdID), fill = species), alpha = 0.5) +
  geom_point(data=dat_comp, aes(x=-mean_CGx_specific, y=-mean_CGz_specific, col = species)) +

  scale_y_continuous(name = "z (% of full body length)",limits = c(-0.11,0.11), breaks = c(-0.1,0,0.1), labels = c(-10,0,10)) +
  scale_x_continuous(name = "x (% of full body length)",limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  coord_fixed() + th +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.1, yend = 0.1) +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

chulls_xy <- ddply(dat_final[,c("species","BirdID","wing_CGx_specific","wing_CGy_specific")], .(species, BirdID),
                   function(df) df[chull(df$wing_CGx_specific, df$wing_CGy_specific), ])

CG_top_view <- ggplot()+
  geom_polygon(data=chulls_xy, aes(x=wing_CGy_specific, y=wing_CGx_specific, group=interaction(species,BirdID), fill = species), alpha = 0.15) +
  geom_point(data=dat_comp, aes(x=mean_wing_CGy_specific, y=mean_wing_CGx_specific, col = species)) +
  scale_y_continuous(name = "x (% of half wing span)",limits = c(-0.20,0.05), breaks = c(-0.2,-0.1,0), labels = c(20,10,0)) +
  scale_x_continuous(name = "y (% of half wing span)",limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  coord_fixed() + th +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.2, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))





q_dot_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = wing_CGy_specific, col = species), pch = 15) + th

del_M_plot <- ggplot()+
  geom_point(data = subset(dat_all_e, BirdID != "21_0203" &species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade)) +
  geom_point(data = subset(dat_all_t, BirdID != "21_0203" &species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade), pch = 15) +
  th +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

cg_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = full_CGz_specific, col = clade)) + th

CG_range_plot <- ggplot()+
  geom_point(data = dat_bird, aes(x = total_bird_mass, y = range_CGx_specific, col = clade)) + th +
  scale_x_continuous(trans='log10')


