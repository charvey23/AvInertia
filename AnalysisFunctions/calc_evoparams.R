## Calculate the evolutionary models for the data
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
# Take the mean of each metric by species
output_data_means <- aggregate(select_if(dat_comp, is.numeric), by = list(species = dat_comp$species), mean)
morpho_data_means <- aggregate(select_if(dat_bird, is.numeric), by = list(species = dat_bird$species), mean)

all_data_means              <- merge(dat_feat_means, output_data_means, id = "species")
all_data_means              <- merge(all_data_means, morpho_data_means, id = "species")
all_data_means              <- merge(all_data_means, unique(dat_comp[,c("species", "binomial")]), id = "species")
all_data_means$torso_length <- all_data_means$torsotail_length - all_data_means$tail_length
all_data_means$torso_mass   <- all_data_means$torsotail_mass - all_data_means$tail_mass_g
all_data_means$leg_mass     <- 0.5*(all_data_means$left_leg_mass_g+all_data_means$right_leg_mass)
all_data_means$body_inertia <- (all_data_means$full_length*all_data_means$max_wingspan*all_data_means$full_m)

# create a matrix from the final data as it is required for fitContinuous()
all_data_means_mat           <- as.matrix(select_if(all_data_means, is.numeric))
rownames(all_data_means_mat) <- all_data_means$binomial
colnames(all_data_means_mat) <- colnames(select_if(all_data_means, is.numeric))

# calculate the phylogenetic signal for the stability metrics
# physignal(all_data_means_mat[,c("max_q_nd","min_q_nd")], pruned_mcc)

BM_xcg = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("mean_CGx_specific_orgShoulder")], model = "BM")
OU_xcg = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("mean_CGx_specific_orgShoulder")], model = "OU")

BM_maxsm = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("max_sm_nd")], model = "BM")
OU_maxsm = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("max_sm_nd")], model = "OU")
BM_minsm = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("min_sm_nd")], model = "BM")
OU_minsm = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("min_sm_nd")], model = "OU")

BM_maxag = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("max_q_nd")], model = "BM")
OU_maxag = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("max_q_nd")], model = "OU")

BM_minag = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("min_q_nd")], model = "BM")
OU_minag = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,c("min_q_nd")], model = "OU")
#
# # Calculate the individual models of OU and BM respectively
#
# morpho_traits = c("head_height","head_length","head_mass",
#                   "body_height_max","body_width_max","torso_length","torso_mass",
#                   "tail_width","tail_length","tail_mass_g",
#                   "humerus_diameter_mm","humerus_length_mm","humerus_mass_g",
#                   "radius_diameter_mm","radius_length_mm","radius_mass_g",
#                   "ulna_diameter_mm","ulna_length_mm","ulna_mass_g",
#                   "cmc_diameter_mm","cmc_length_mm","cmc_mass_g",
#                   "feat_length.P","feat_width.P","feat_mass.P",
#                   "feat_length.S","feat_width.S","feat_mass.S")
# check_OU <- as.data.frame(matrix(nrow = 28,ncol = 6))
# colnames(check_OU) <- c("trait","aicc","del_aicc","sigsq","alpha","w")
# for (i in 1:length(morpho_traits)){
#   BM_model = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,morpho_traits[i]], model = "BM")
#   OU_model = fitContinuous(phy = pruned_mcc, dat = all_data_means_mat[,morpho_traits[i]], model = "OU")
#   delta <- c(BM_model$opt$aicc,OU_model$opt$aicc) - min(c(BM_model$opt$aicc,OU_model$opt$aicc))
#   L <- exp(-0.5 * delta)            # relative likelihoods of models
#   w <- L/sum(L)                     # Akaike weights
#   # Save all OU data
#   check_OU$trait[i] = morpho_traits[i]
#   check_OU$aicc[i]  = OU_model$opt$aicc
#   check_OU$del_aicc[i]  = BM_model$opt$aicc-OU_model$opt$aicc
#   check_OU$sigsq[i] = OU_model$opt$sigsq
#   check_OU$alpha[i] = OU_model$opt$alpha
#   check_OU$w[i]     = w[2]
# }
#

# ----- Below is code to calculate the multivariate BM results for each body component ----
#
# source("influ_continuous_BMmult.R")
# # sets all birds to be within the same group
# gp.end        = rep(1,22)
# names(gp.end) = pruned_mcc$tip.label
# #compare.evol.rates(A=all_data_means_mat[,c("head_height","head_length","head_mass")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# # both computes the final rate and the influence per species
# dat_rate_head  <- influ_continuous_BMmult(data=all_data_means_mat[,c("head_height","head_length","head_mass")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_torso <- influ_continuous_BMmult(data=all_data_means_mat[,c("body_height_max","body_width_max","torso_length","torso_mass")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_tail  <- influ_continuous_BMmult(data=all_data_means_mat[,c("tail_width","tail_length","tail_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_hum   <- influ_continuous_BMmult(data=all_data_means_mat[,c("humerus_diameter_mm","humerus_length_mm","humerus_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_rad   <- influ_continuous_BMmult(data=all_data_means_mat[,c("radius_diameter_mm","radius_length_mm","radius_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_uln   <- influ_continuous_BMmult(data=all_data_means_mat[,c("ulna_diameter_mm","ulna_length_mm","ulna_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_cmc   <- influ_continuous_BMmult(data=all_data_means_mat[,c("cmc_diameter_mm","cmc_length_mm","cmc_mass_g")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_P     <- influ_continuous_BMmult(data=all_data_means_mat[,c("feat_length.P","feat_width.P","feat_mass.P")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_S     <- influ_continuous_BMmult(data=all_data_means_mat[,c("feat_length.S","feat_width.S","feat_mass.S")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_CG    <- influ_continuous_BMmult(data=all_data_means_mat[,c("mean_CGx_orgBeak","mean_CGz_orgDorsal")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_CG_sh <- influ_continuous_BMmult(data=all_data_means_mat[,c("mean_CGx_orgShoulder","mean_CGz_orgShoulder")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_CG_wing    <- influ_continuous_BMmult(data=all_data_means_mat[,c("max_wing_CGx","max_wing_CGy")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
# dat_rate_fullbird <- influ_continuous_BMmult(data=all_data_means_mat[,c("full_length","max_wingspan","full_m")], phy=pruned_mcc, method="simulation",gp=gp.end,iter=999)
#
# dat_var_tot           <- data.frame(matrix(nrow = 0, ncol = 2))
# colnames(dat_var_tot) <- c("component","mean_sig_sq")
#
# dat_var_tot = rbind(dat_var_tot,list(component = "head", mean_sig_sq = dat_rate_head$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "torso",mean_sig_sq = dat_rate_torso$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "tail",mean_sig_sq = dat_rate_tail$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "humerus",mean_sig_sq = dat_rate_hum$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "radius",mean_sig_sq = dat_rate_rad$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "ulna",mean_sig_sq = dat_rate_uln$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "cmc",mean_sig_sq = dat_rate_cmc$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "full_CG",mean_sig_sq = dat_rate_CG$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "full_CG_sh",mean_sig_sq = dat_rate_CG_sh$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "full_I",mean_sig_sq = dat_rate_I$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "wing_CG",mean_sig_sq = dat_rate_CG_wing$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "P",mean_sig_sq = dat_rate_P$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "S",mean_sig_sq = dat_rate_S$full.model.estimates$sigsq))
# dat_var_tot = rbind(dat_var_tot,list(component = "full_bird",mean_sig_sq = dat_rate_fullbird$full.model.estimates$sigsq))
#
# dat_var_range  <- data.frame(matrix(nrow = 0, ncol = 3))
# colnames(dat_var_range) <- c("component","species","mean_sig_sq")
#
# dat_var_range = rbind(dat_var_range,list(component = rep("head",22), species = dat_rate_head$sensi.estimates$species, mean_sig_sq = dat_rate_head$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("torso",22), species = dat_rate_torso$sensi.estimates$species, mean_sig_sq = dat_rate_torso$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("tail",22), species = dat_rate_tail$sensi.estimates$species, mean_sig_sq = dat_rate_tail$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("humerus",22), species = dat_rate_hum$sensi.estimates$species, mean_sig_sq = dat_rate_hum$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("radius",22), species = dat_rate_rad$sensi.estimates$species, mean_sig_sq = dat_rate_rad$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("ulna",22), species = dat_rate_uln$sensi.estimates$species, mean_sig_sq = dat_rate_uln$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("cmc",22), species = dat_rate_cmc$sensi.estimates$species, mean_sig_sq = dat_rate_cmc$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("full_CG",22), species = dat_rate_CG$sensi.estimates$species, mean_sig_sq = dat_rate_CG$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("full_CG_sh",22), species = dat_rate_CG_sh$sensi.estimates$species, mean_sig_sq = dat_rate_CG_sh$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("full_I",22), species = dat_rate_I$sensi.estimates$species, mean_sig_sq = dat_rate_I$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("P",22), species = dat_rate_P$sensi.estimates$species, mean_sig_sq = dat_rate_P$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("S",22), species = dat_rate_S$sensi.estimates$species, mean_sig_sq = dat_rate_S$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("wing_CG",22), species = dat_rate_CG_wing$sensi.estimates$species, mean_sig_sq = dat_rate_CG_wing$sensi.estimates$sigsq))
# dat_var_range = rbind(dat_var_range,list(component = rep("full_bird",22), species = dat_rate_fullbird$sensi.estimates$species, mean_sig_sq = dat_rate_fullbird$sensi.estimates$sigsq))
#
# dat_var_range$component = factor(dat_var_range$component, levels = rev(c("full_bird","head","torso","tail","humerus","ulna","radius","cmc","P","S","full_CG","full_CG_sh","full_I","wing_CG")))

