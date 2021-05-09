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
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, pr = TRUE, pl = TRUE)
CGx_model_mcmc_output  = summary(CGx_model_mcmc)

CGx_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(-mean_CGx_specific_orgBeak) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGx_sp_model_mcmc_output  = summary(CGx_sp_model_mcmc)

CGy_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_wing_CGy) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGy_model_mcmc_output  = summary(CGy_model_mcmc)

CGy_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_wing_CGy_specific) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGy_sp_model_mcmc_output  = summary(CGy_sp_model_mcmc)

CGz_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_CGz_orgDorsal) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, pr = TRUE, pl = TRUE)
CGz_model_mcmc_output  = summary(CGz_model_mcmc)

CGz_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(mean_CGz_specific_orgDorsal) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
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
    data = dat_comp,
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
    data = dat_comp,
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
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
Izz_model_mcmc_output  = summary(Izz_model_mcmc)


## ------------------------------------------------------------------------------------
## --------------------------------- Figure 4 -----------------------------------------
## ------------------------------------------------------------------------------------
tmp        <- subset(dat_feat, feather != "alula")
tmp$length <- tmp$l_vane + tmp$l_cal
tmp$width  <- tmp$w_cal + tmp$w_vd + tmp$w_vp
tmp$comb <- paste(tmp$species,tmp$BirdID,sep = "_")
tmp <- subset(tmp, species %in% c("acc_str","acc_coo","chr_amh","lop_nyc","lop_imp","cyp_nig","bra_can","fal_per","lar_gla","pel_ery","cor_cor") |
                comb == "tyt_alb_19_265" | comb == "meg_alc_20_3495" | comb == "fal_col_20_1016" | comb == "ana_pla_20_0291" | comb == "cho_min_20_1027" | comb == "aec_occ_20_0922" |
                comb == "cya_ste_20_0925" | comb == "oce_leu_20_1015" | comb == "ard_her_20_0284" | comb == "col_liv_20_0303" | comb == "col_aur_20_0907")
tmp$grouping = NA
tmp$grouping[which(grepl("P", tmp$feather, fixed = TRUE))] = "P"
tmp$grouping[which(grepl("S", tmp$feather, fixed = TRUE))] = "S"

dat_feat_means <- aggregate(list(feat_length = tmp$length,
                                 feat_width = tmp$width,
                                 feat_mass = tmp$m_f), by = list(species = tmp$species, grouping = tmp$grouping), mean)

dat_feat_means <- reshape(dat_feat_means,idvar = "species",timevar = "grouping",direction = "wide")

output_data_means <- aggregate(select_if(dat_comp, is.numeric), by = list(species = dat_comp$species), mean)
morpho_data_means <- aggregate(select_if(dat_bird, is.numeric), by = list(species = dat_bird$species), mean)

all_data_means              <- merge(dat_feat_means, output_data_means, id = "species")
all_data_means              <- merge(all_data_means, morpho_data_means, id = "species")
all_data_means              <- merge(all_data_means, unique(dat_comp[,c("species", "binomial")]), id = "species")
all_data_means$torso_length <- all_data_means$torsotail_length - all_data_means$tail_length
all_data_means$torso_mass   <- all_data_means$torsotail_mass - all_data_means$tail_mass_g
all_data_means$leg_mass     <- 0.5*(all_data_means$left_leg_mass_g+all_data_means$right_leg_mass)

test  <- aggregate(list(mean_chord = sqrt((dat_wing$pt12_X - dat_wing$pt11_X)^2+(dat_wing$pt12_Y - dat_wing$pt11_Y)^2+(dat_wing$pt12_Z - dat_wing$pt11_Z)^2)), by = list(species = dat_wing$species, BirdID = dat_wing$BirdID), max)
test2 <- merge(test,dat_comp, id = c("species","BirdID"))
test2$c_l      = test2$mean_chord/test2$full_length
test2$c_l_true = abs(0.25*test2$mean_chord)/abs(test2$mean_CGx_specific_orgShoulder*test2$full_length)
test2[,c("species","mean_CGx_specific_orgShoulder","c_l","c_l_true")]
test2$species_order = factor(test2$species, levels = phylo_order$species)

# create a matrix from the final data as it is required for fitContinuous()
all_data_means_mat           <- as.matrix(select_if(all_data_means, is.numeric))
rownames(all_data_means_mat) <- all_data_means$binomial
colnames(all_data_means_mat) <- colnames(select_if(all_data_means, is.numeric))
# log transform all data - absolute is just for CGx which is negative
all_data_means_mat <- apply(all_data_means_mat,2,function(x) log(abs(x)))

source("influ_continuous_BMmult.R")
# sets all birds to be within the same group
gp.end        = rep(1,22)
names(gp.end) = pruned_mcc$tip.label
#compare.evol.rates(A=all_data_means_mat[,c("head_height","head_length","head_mass")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# both computes the final rate and the influence per species
dat_rate_head  <- influ_continuous_BMmult(data=all_data_means_mat[,c("head_height","head_length","head_mass")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_torso <- influ_continuous_BMmult(data=all_data_means_mat[,c("body_height_max","body_width_max","torso_length","torso_mass")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_tail  <- influ_continuous_BMmult(data=all_data_means_mat[,c("tail_width","tail_length","tail_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_hum   <- influ_continuous_BMmult(data=all_data_means_mat[,c("humerus_diameter_mm","humerus_length_mm","humerus_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_rad   <- influ_continuous_BMmult(data=all_data_means_mat[,c("radius_diameter_mm","radius_length_mm","radius_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_uln   <- influ_continuous_BMmult(data=all_data_means_mat[,c("ulna_diameter_mm","ulna_length_mm","ulna_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_cmc   <- influ_continuous_BMmult(data=all_data_means_mat[,c("cmc_diameter_mm","cmc_length_mm","cmc_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_P     <- influ_continuous_BMmult(data=all_data_means_mat[,c("feat_length.P","feat_width.P","feat_mass.P")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_S     <- influ_continuous_BMmult(data=all_data_means_mat[,c("feat_length.S","feat_width.S","feat_mass.S")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_CG    <- influ_continuous_BMmult(data=all_data_means_mat[,c("mean_CGx_orgShoulder","mean_CGz_orgDorsal")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_I     <- influ_continuous_BMmult(data=all_data_means_mat[,c("mean_Ixx","mean_Iyy","mean_Izz")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_CG_sp <- influ_continuous_BMmult(data=all_data_means_mat[,c("mean_CGx_specific_orgShoulder","mean_CGz_specific_orgDorsal")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_I_sp  <- influ_continuous_BMmult(data=all_data_means_mat[,c("mean_Ixx_specific","mean_Iyy_specific","mean_Izz_specific")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_CG_wing    <- influ_continuous_BMmult(data=all_data_means_mat[,c("max_wing_CGx","max_wing_CGy")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
dat_rate_CG_wing_sp <- influ_continuous_BMmult(data=all_data_means_mat[,c("max_wing_CGx_specific","max_wing_CGy_specific")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)

dat_var_tot           <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dat_var_tot) <- c("component","mean_sig_sq")

dat_var_tot = rbind(dat_var_tot,list(component = "head", mean_sig_sq = dat_rate_head$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "torso",mean_sig_sq = dat_rate_torso$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "tail",mean_sig_sq = dat_rate_tail$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "humerus",mean_sig_sq = dat_rate_hum$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "radius",mean_sig_sq = dat_rate_rad$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "ulna",mean_sig_sq = dat_rate_uln$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "cmc",mean_sig_sq = dat_rate_cmc$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "full_CG",mean_sig_sq = dat_rate_CG$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "full_I",mean_sig_sq = dat_rate_I$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "full_CG_sp",mean_sig_sq = dat_rate_CG_sp$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "full_I_sp",mean_sig_sq = dat_rate_I_sp$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "wing_CG",mean_sig_sq = dat_rate_CG_wing$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "wing_CG_sp",mean_sig_sq = dat_rate_CG_wing_sp$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "P",mean_sig_sq = dat_rate_P$full.model.estimates$sigsq))
dat_var_tot = rbind(dat_var_tot,list(component = "S",mean_sig_sq = dat_rate_S$full.model.estimates$sigsq))

dat_var_range  <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(dat_var_range) <- c("component","species","mean_sig_sq")

dat_var_range = rbind(dat_var_range,list(component = rep("head",22), species = dat_rate_head$sensi.estimates$species, mean_sig_sq = dat_rate_head$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("torso",22), species = dat_rate_torso$sensi.estimates$species, mean_sig_sq = dat_rate_torso$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("tail",22), species = dat_rate_tail$sensi.estimates$species, mean_sig_sq = dat_rate_tail$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("humerus",22), species = dat_rate_hum$sensi.estimates$species, mean_sig_sq = dat_rate_hum$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("radius",22), species = dat_rate_rad$sensi.estimates$species, mean_sig_sq = dat_rate_rad$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("ulna",22), species = dat_rate_uln$sensi.estimates$species, mean_sig_sq = dat_rate_uln$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("cmc",22), species = dat_rate_cmc$sensi.estimates$species, mean_sig_sq = dat_rate_cmc$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("full_CG",22), species = dat_rate_CG$sensi.estimates$species, mean_sig_sq = dat_rate_CG$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("full_I",22), species = dat_rate_I$sensi.estimates$species, mean_sig_sq = dat_rate_I$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("full_CG_sp",22), species = dat_rate_CG_sp$sensi.estimates$species, mean_sig_sq = dat_rate_CG_sp$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("full_I_sp",22), species = dat_rate_I_sp$sensi.estimates$species, mean_sig_sq = dat_rate_I_sp$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("P",22), species = dat_rate_P$sensi.estimates$species, mean_sig_sq = dat_rate_P$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("S",22), species = dat_rate_S$sensi.estimates$species, mean_sig_sq = dat_rate_S$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("wing_CG",22), species = dat_rate_CG_wing$sensi.estimates$species, mean_sig_sq = dat_rate_CG_wing$sensi.estimates$sigsq))
dat_var_range = rbind(dat_var_range,list(component = rep("wing_CG_sp",22), species = dat_rate_CG_wing_sp$sensi.estimates$species, mean_sig_sq = dat_rate_CG_wing_sp$sensi.estimates$sigsq))

dat_var_range$component = factor(dat_var_range$component, levels = c("head","torso","tail","humerus","ulna","radius","cmc","P","S","full_CG","full_CG_sp","wing_CG","wing_CG_sp","full_I","full_I_sp"))

fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("max_wing_CGy_specific")], model = "OU")
library(mvMORPH)
fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("mean_CGx_specific_orgShoulder")], model = "OU")
mvOU(tree = pruned_mcc, data = all_data_means_mat[,c("mean_CGx_specific_orgShoulder")])
mvOU(tree = pruned_mcc, data = all_data_means_mat[,c("max_wing_CGy_specific")])


library(pmc)

# to verify that the OU model is an appropriate fit given the size of our data
#out <- pmc(tree = pruned_mcc, data = all_data_means_mat[,c("mean_CGx_specific_orgShoulder")], "BM", "OU", nboot = 1000, mc.cores = 1)
filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_pmc_output.RData",sep="")
save(out,file = filename)

library("ggplot2")
library("tidyr")
library("dplyr")
dists <- data.frame(null = out$null, test = out$test)
dists %>%
  gather(dist, value) %>%
  ggplot(aes(value, fill = dist)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out$lr)
