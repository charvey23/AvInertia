### --- Script containing the data manipulations needed to create all plots -----


## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 1 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------
chulls_elbman <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                       function(df) df[chull(df$elbow, df$manus), ])


## --------------------------------------------------------------------------------------------------------------------
## ------------------------------------------- Figure 2 ---------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------

## ---------------- Panel B data -----------

# Compute the alphahull
for (i in 1:nrow(dat_bird)){
  curr_species = dat_bird$species[i]
  curr_BirdID = dat_bird$BirdID[i]
  tmp = dat_final[which(dat_final$species == curr_species & dat_final$BirdID == curr_BirdID), c("full_CGx_specific_orgShoulder", "full_CGz_specific_orgShoulder")]
  # fit the convex hull with an alpha factor
  alphashape <- ahull(tmp, alpha = 0.1)
  # save all the given vertices
  vertices <- as.data.frame(tmp[alphashape$ashape.obj$alpha.extremes,])
  # Need to order the points appropriately
  # calculate the mean value
  centerpt <- c(mean(vertices[,1]),mean(vertices[,2]))
  # calculate the angle from the mean of each point
  vertices$angle <- atan2(vertices[,2]-centerpt[2],vertices[,1]-centerpt[1])
  # sort by the angle
  vertices <- vertices[order(vertices$angle),]
  vertices$species = curr_species
  vertices$BirdID = curr_BirdID
  if(i == 1){
    vertices_cgxz = vertices
  }else{
    vertices_cgxz = rbind(vertices_cgxz,vertices)
  }
}
vertices_cgxz$species_order = factor(vertices_cgxz$species, levels = phylo_order$species)

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

## ---------------- Panel C data -----------

# Compute the alphahull
for (i in 1:nrow(dat_bird)){
  curr_species = dat_bird$species[i]
  curr_BirdID = dat_bird$BirdID[i]
  tmp = dat_final[which(dat_final$species == curr_species & dat_final$BirdID == curr_BirdID), c("wing_CGx_specific_orgShoulder", "wing_CGy_specific_orgShoulder")]
  # fit the convex hull with an alpha factor
  alphashape <- ahull(tmp, alpha = 0.1)
  # save all the given vertices
  vertices <- as.data.frame(tmp[alphashape$ashape.obj$alpha.extremes,])
  # Need to order the points appropriately
  # calculate the mean value
  centerpt <- c(mean(vertices[,1]),mean(vertices[,2]))
  # calculate the angle from the mean of each point
  vertices$angle <- atan2(vertices[,2]-centerpt[2],vertices[,1]-centerpt[1])
  # sort by the angle
  vertices <- vertices[order(vertices$angle),]
  vertices$species = curr_species
  vertices$BirdID = curr_BirdID
  if(i == 1){
    vertices_cgxy = vertices
  }else{
    vertices_cgxy = rbind(vertices_cgxy,vertices)
  }
}
vertices_cgxy$species_order = factor(vertices_cgxy$species, levels = phylo_order$species)


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

## ------------------------------------------------------------------------------------
## --------------------------------- Figure 3 -----------------------------------------
## ------------------------------------------------------------------------------------

dat_I_sp_plot <- merge(aggregate(list(min_Ixx_specific = dat_comp$min_Ixx_specific,
                                      min_Iyy_specific = dat_comp$min_Iyy_specific,
                                      min_Izz_specific = dat_comp$min_Izz_specific,
                                      min_Ixz_specific = dat_comp$min_Ixz_specific),  by=list(species = dat_comp$species), min),
                       aggregate(list(max_Ixx_specific = dat_comp$max_Ixx_specific,
                                      max_Iyy_specific = dat_comp$max_Iyy_specific,
                                      max_Izz_specific = dat_comp$max_Izz_specific,
                                      max_Ixz_specific = dat_comp$max_Ixz_specific),  by=list(species = dat_comp$species), max), by = "species")

## -------------- I - isometry -------------
Ixx_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_Ixx) ~ log(full_m),
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
Ixx_model_mcmc_output  = summary(Ixx_model_mcmc)

Iyy_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_Iyy) ~ log(full_m),
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
Iyy_model_mcmc_output  = summary(Iyy_model_mcmc)

Izz_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_Izz) ~ log(full_m),
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
Izz_model_mcmc_output  = summary(Izz_model_mcmc)


## ------------------------------------------------------------------------------------
## --------------------------------- Figure 4 -----------------------------------------
## ------------------------------------------------------------------------------------

tmp <- melt(dat_feat, id = c("species","feather","BirdID"))
tmp$comb <- paste(tmp$species,tmp$BirdID,sep = "_")
tmp <- subset(tmp, species %in% c("acc_str","acc_coo","chr_amh","lop_nyc","lop_imp","cyp_nig","bra_can","fal_per","lar_gla","pel_ery","cor_cor") |
                comb == "tyt_alb_19_265" | comb == "meg_alc_20_3495" | comb == "fal_col_20_1016" | comb == "ana_pla_20_0291" | comb == "cho_min_20_1027" | comb == "aec_occ_20_0922" |
                comb == "cya_ste_20_0925" | comb == "oce_leu_20_1015" | comb == "ard_her_20_0284" | comb == "col_liv_20_0303" | comb == "col_aur_20_0907")
dat_feat_means <- dcast(tmp[,c("species", "feather","variable","value")], species ~ feather + variable, value.var = "value", fun.aggregate = mean)
dat_feat_means <- dat_feat_means[ , apply(dat_feat_means, 2, function(x) !any(is.na(x)))]

output_data_means <- aggregate(select_if(dat_comp, is.numeric), by = list(species = dat_comp$species), mean)
morpho_data_means <- aggregate(select_if(dat_bird, is.numeric), by = list(species = dat_bird$species), mean)

all_data_means              <- merge(dat_feat_means, output_data_means, id = "species")
all_data_means              <- merge(all_data_means, morpho_data_means, id = "species")
all_data_means              <- merge(all_data_means, unique(dat_comp[,c("species", "binomial")]), id = "species")
all_data_means$torso_length <- all_data_means$torsotail_length - all_data_means$tail_length
all_data_means$torso_mass   <- all_data_means$torsotail_mass - all_data_means$tail_mass_g
all_data_means$leg_mass     <- 0.5*(all_data_means$left_leg_mass_g+all_data_means$right_leg_mass)

# create a matrix from the final data as it is required for fitContinuous()
all_data_means_mat           <- as.matrix(select_if(all_data_means, is.numeric))
rownames(all_data_means_mat) <- all_data_means$binomial
colnames(all_data_means_mat) <- colnames(select_if(all_data_means, is.numeric))

# list all colnames that we want to fit model to
model_names <- c("head_height","head_length","head_mass",
                 "body_height_max","body_width_max","torso_length","torso_mass",
                 "tail_width","tail_length","tail_mass_g",
                 "leg_mass",
                 # Wing bones
                 "humerus_diameter_mm","humerus_length_mm","humerus_mass_g",
                 "radius_diameter_mm","radius_length_mm","radius_mass_g",
                 "ulna_diameter_mm","ulna_length_mm","ulna_mass_g",
                 "cmc_diameter_mm","cmc_length_mm","cmc_mass_g",
                 # Wing muscle groups
                 "brachial_muscle_mass","antebrachial_muscle_mass","manus_muscle_mass",
                 #skin/coverts/terts
                 "all_skin_coverts_mass","tertiary_mass",
                 #full CG data
                 "mean_CGx_orgBeak","mean_CGz_orgDorsal","mean_CGx_specific_orgBeak","mean_CGz_specific_orgDorsal",
                 #full I data
                 "mean_Ixx","mean_Iyy","mean_Izz",
                 "mean_Ixx_specific","mean_Iyy_specific","mean_Izz_specific")

dat_var <- fitContinuous(phy = pruned_mcc, dat = log(abs(all_data_means_mat[,model_names])), model ="BM")
dat_var_feat <- fitContinuous(phy = pruned_mcc, dat = log(abs(all_data_means_mat[,4:123])), model ="BM")

dat_var_tot  <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(dat_var_tot) <- c("component","mean_sig_sq","min_val")

dat_var_tot = rbind(dat_var_tot,list(component = "fakeforcolour",
                                     mean_sig_sq = 0,
                                     min_val = 0.01))

dat_var_tot = rbind(dat_var_tot,list(component = "head",
                                     mean_sig_sq = mean(dat_var$head_height$opt$sigsq,
                                                        dat_var$head_length$opt$sigsq,
                                                        dat_var$head_mass$opt$sigsq),
                                     min_val = min(c(all_data_means$head_height,all_data_means$head_length,all_data_means$head_mass))))

tmp = mean(dat_var$body_width_max$opt$sigsq,dat_var$body_height_max$opt$sigsq)
dat_var_tot = rbind(dat_var_tot,list(component = "torso",
                                     mean_sig_sq = mean(tmp,
                                                        dat_var$torso_length$opt$sigsq,
                                                        dat_var$torso_mass$opt$sigsq),
                                     min_val = min(c(all_data_means$body_width_max,all_data_means$torso_length,all_data_means$torso_mass))))
dat_var_tot = rbind(dat_var_tot,list(component = "tail",
                                     mean_sig_sq = mean(dat_var$tail_width$opt$sigsq,
                                                        dat_var$tail_length$opt$sigsq,
                                                        dat_var$tail_mass_g$opt$sigsq),
                                     min_val = min(c(all_data_means$tail_width,all_data_means$tail_length,all_data_means$tail_mass_g)))
)
dat_var_tot = rbind(dat_var_tot,list(component = "humerus",
                                     mean_sig_sq = mean(dat_var$humerus_diameter_mm$opt$sigsq,
                                                        dat_var$humerus_length_mm$opt$sigsq,
                                                        dat_var$humerus_mass_g$opt$sigsq),
                                     min_val = min(c(all_data_means$humerus_diameter_mm,all_data_means$humerus_length_mm,all_data_means$humerus_mass_g))))
dat_var_tot = rbind(dat_var_tot,list(component = "radius",
                                     mean_sig_sq = mean(dat_var$radius_diameter_mm$opt$sigsq,
                                                        dat_var$radius_length_mm$opt$sigsq,
                                                        dat_var$radius_mass_g$opt$sigsq),
                                     min_val = min(c(all_data_means$radius_diameter_mm,all_data_means$radius_length_mm,all_data_means$radius_mass_g))))
dat_var_tot = rbind(dat_var_tot,list(component = "ulna",
                                     mean_sig_sq = mean(dat_var$ulna_diameter_mm$opt$sigsq,
                                                        dat_var$ulna_length_mm$opt$sigsq,
                                                        dat_var$ulna_mass_g$opt$sigsq),
                                     min_val = min(c(all_data_means$ulna_diameter_mm,all_data_means$ulna_length_mm,all_data_means$ulna_mass_g))))
dat_var_tot = rbind(dat_var_tot,list(component = "cmc",
                                     mean_sig_sq = mean(dat_var$cmc_diameter_mm$opt$sigsq,
                                                        dat_var$cmc_length_mm$opt$sigsq,
                                                        dat_var$cmc_mass_g$opt$sigsq),
                                     min_val = min(c(all_data_means$cmc_diameter_mm,all_data_means$cmc_length_mm,all_data_means$cmc_mass_g))))
dat_var_tot = rbind(dat_var_tot,list(component = "muscles_brach",
                                     mean_sig_sq = dat_var$brachial_muscle_mass$opt$sigsq,
                                     min_val = min(c(all_data_means$brachial_muscle_mass))))
dat_var_tot = rbind(dat_var_tot,list(component = "muscles_abrach",
                                     mean_sig_sq = dat_var$antebrachial_muscle_mass$opt$sigsq,
                                     min_val = min(c(all_data_means$antebrachial_muscle_mass))))
dat_var_tot = rbind(dat_var_tot,list(component = "muscles_manus",
                                     mean_sig_sq = dat_var$manus_muscle_mass$opt$sigsq,
                                     min_val = min(c(all_data_means$manus_muscle_mass))))
dat_var_tot = rbind(dat_var_tot,list(component = "tertiaries",
                                     mean_sig_sq = dat_var$tertiary_mass$opt$sigsq,
                                     min_val = min(c(all_data_means$tertiary_mass))))
dat_var_tot = rbind(dat_var_tot,list(component = "skin/coverts",
                                     mean_sig_sq = dat_var$all_skin_coverts_mass$opt$sigsq,
                                     min_val = min(c(all_data_means$all_skin_coverts_mass))))
dat_var_tot = rbind(dat_var_tot,list(component = "full_CG",
                                     mean_sig_sq = mean(dat_var$mean_CGx_orgBeak$opt$sigsq,
                                                        dat_var$mean_CGz_orgDorsal$opt$sigsq),
                                     min_val = min(abs(c(all_data_means$mean_CGx_orgBeak,all_data_means$mean_CGz_orgDorsal)))))
dat_var_tot = rbind(dat_var_tot,list(component = "full_CG_specific",
                                     mean_sig_sq = mean(dat_var$mean_CGx_specific_orgBeak$opt$sigsq,
                                                        dat_var$mean_CGz_specific_orgDorsal$opt$sigsq),
                                     min_val = min(abs(c(all_data_means$mean_CGx_specific_orgBeak,all_data_means$mean_CGz_specific_orgDorsal)))))
dat_var_tot = rbind(dat_var_tot,list(component = "full_I",
                                     mean_sig_sq = mean(dat_var$mean_Ixx$opt$sigsq,
                                                        dat_var$mean_Iyy$opt$sigsq,
                                                        dat_var$mean_Izz$opt$sigsq),
                                     min_val = min(c(all_data_means$mean_Ixx,all_data_means$mean_Iyy,all_data_means$mean_Izz))))
dat_var_tot = rbind(dat_var_tot,list(component = "full_I_specific",
                                     mean_sig_sq = mean(dat_var$mean_Ixx_specific$opt$sigsq,
                                                        dat_var$mean_Iyy_specific$opt$sigsq,
                                                        dat_var$mean_Izz_specific$opt$sigsq),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P10",
                                     mean_sig_sq = mean(dat_var_feat$P10_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P10_l_vane$opt$sigsq,dat_var_feat$P10_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P10_w_cal$opt$sigsq,dat_var_feat$P10_w_vd$opt$sigsq,dat_var_feat$P10_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P9",
                                     mean_sig_sq = mean(dat_var_feat$P9_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P9_l_vane$opt$sigsq,dat_var_feat$P9_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P9_w_cal$opt$sigsq,dat_var_feat$P9_w_vd$opt$sigsq,dat_var_feat$P9_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P8",
                                     mean_sig_sq = mean(dat_var_feat$P8_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P8_l_vane$opt$sigsq,dat_var_feat$P8_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P8_w_cal$opt$sigsq,dat_var_feat$P8_w_vd$opt$sigsq,dat_var_feat$P8_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P7",
                                     mean_sig_sq = mean(dat_var_feat$P7_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P7_l_vane$opt$sigsq,dat_var_feat$P7_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P7_w_cal$opt$sigsq,dat_var_feat$P7_w_vd$opt$sigsq,dat_var_feat$P7_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P6",
                                     mean_sig_sq = mean(dat_var_feat$P6_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P6_l_vane$opt$sigsq,dat_var_feat$P6_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P6_w_cal$opt$sigsq,dat_var_feat$P6_w_vd$opt$sigsq,dat_var_feat$P6_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P5",
                                     mean_sig_sq = mean(dat_var_feat$P5_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P5_l_vane$opt$sigsq,dat_var_feat$P5_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P5_w_cal$opt$sigsq,dat_var_feat$P5_w_vd$opt$sigsq,dat_var_feat$P5_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))
dat_var_tot = rbind(dat_var_tot,list(component = "P4",
                                     mean_sig_sq = mean(dat_var_feat$P4_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P4_l_vane$opt$sigsq,dat_var_feat$P4_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P4_w_cal$opt$sigsq,dat_var_feat$P4_w_vd$opt$sigsq,dat_var_feat$P4_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P3",
                                     mean_sig_sq = mean(dat_var_feat$P3_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P3_l_vane$opt$sigsq,dat_var_feat$P3_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P3_w_cal$opt$sigsq,dat_var_feat$P3_w_vd$opt$sigsq,dat_var_feat$P3_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P2",
                                     mean_sig_sq = mean(dat_var_feat$P2_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P2_l_vane$opt$sigsq,dat_var_feat$P2_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P2_w_cal$opt$sigsq,dat_var_feat$P2_w_vd$opt$sigsq,dat_var_feat$P2_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "P1",
                                     mean_sig_sq = mean(dat_var_feat$P1_m_f$opt$sigsq,
                                                        mean(dat_var_feat$P1_l_vane$opt$sigsq,dat_var_feat$P1_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$P1_w_cal$opt$sigsq,dat_var_feat$P1_w_vd$opt$sigsq,dat_var_feat$P1_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "S1",
                                     mean_sig_sq = mean(dat_var_feat$S1_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S1_l_vane$opt$sigsq,dat_var_feat$S1_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S1_w_cal$opt$sigsq,dat_var_feat$S1_w_vd$opt$sigsq,dat_var_feat$S1_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "S2",
                                     mean_sig_sq = mean(dat_var_feat$S2_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S2_l_vane$opt$sigsq,dat_var_feat$S2_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S2_w_cal$opt$sigsq,dat_var_feat$S2_w_vd$opt$sigsq,dat_var_feat$S2_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "S3",
                                     mean_sig_sq = mean(dat_var_feat$S3_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S3_l_vane$opt$sigsq,dat_var_feat$S3_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S3_w_cal$opt$sigsq,dat_var_feat$S3_w_vd$opt$sigsq,dat_var_feat$S3_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))
dat_var_tot = rbind(dat_var_tot,list(component = "S4",
                                     mean_sig_sq = mean(dat_var_feat$S4_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S4_l_vane$opt$sigsq,dat_var_feat$S4_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S4_w_cal$opt$sigsq,dat_var_feat$S4_w_vd$opt$sigsq,dat_var_feat$S4_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "S5",
                                     mean_sig_sq = mean(dat_var_feat$S5_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S5_l_vane$opt$sigsq,dat_var_feat$S5_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S5_w_cal$opt$sigsq,dat_var_feat$S5_w_vd$opt$sigsq,dat_var_feat$S5_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "S6",
                                     mean_sig_sq = mean(dat_var_feat$S6_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S6_l_vane$opt$sigsq,dat_var_feat$S6_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S6_w_cal$opt$sigsq,dat_var_feat$S6_w_vd$opt$sigsq,dat_var_feat$S6_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

dat_var_tot = rbind(dat_var_tot,list(component = "S7",
                                     mean_sig_sq = mean(dat_var_feat$S7_m_f$opt$sigsq,
                                                        mean(dat_var_feat$S7_l_vane$opt$sigsq,dat_var_feat$S7_l_cal$opt$sigsq),
                                                        mean(dat_var_feat$S7_w_cal$opt$sigsq,dat_var_feat$S7_w_vd$opt$sigsq,dat_var_feat$S7_w_vp$opt$sigsq)),
                                     min_val = min(c(all_data_means$mean_Ixx_specific,all_data_means$mean_Iyy_specific,all_data_means$mean_Izz_specific))))

