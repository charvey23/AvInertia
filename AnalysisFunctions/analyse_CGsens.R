path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

library(AvInertia)

### -------------------------------------------------------------------------
### -------------------- CG Error Sensitivity Analysis ----------------------
### -------------------------------------------------------------------------

# the more negative the lower stability this is the maximum agility with the most stability
dat_sens_max_ag   = dat_final[which(dat_final$prop_q_dot_nd %in% dat_comp$min_q_nd),]
dat_sens_max_stab = dat_final[which(dat_final$sm_nd %in% dat_comp$max_sm_nd),]
dat_sens_min_stab = dat_final[which(dat_final$sm_nd %in% dat_comp$min_sm_nd),]

shift <- seq(-0.15,0.15, by = 0.01)
count = 1
dat_sens_out = as.data.frame(matrix(0, nrow = 10000, ncol = 0))

# Note: some results with different shifts will give the same final CG placement
# due to the optimization routine with the body density calculation. Need to show
# maximum shift for an error in the input measurement separately from the shift
# due to a tail metric

for(m in 1:36){
  species_curr = dat_bird$species[m]
  birdid_curr  = dat_bird$BirdID[m]

  dat_bird_curr = subset(dat_bird, species == species_curr & BirdID == birdid_curr)

  true_CGx      = dat_bird_curr$x_loc_TorsotailCoG
  for (l in 1:3){
    if(l == 1){
      dat_wing_curr = subset(dat_sens_max_ag, species == species_curr & BirdID == birdid_curr)
    }else{
      if(l == 2){
        dat_wing_curr = subset(dat_sens_max_stab, species == species_curr & BirdID == birdid_curr)
      }else{
        dat_wing_curr = subset(dat_sens_min_stab, species == species_curr & BirdID == birdid_curr)
      }
    }

    for (i in 1:length(shift)){

      # --- Compute the new torso measurements for the given data ---
      dat_bird_curr$x_loc_TorsotailCoG = true_CGx + shift[i]*dat_wing_curr$torso_length
      if(shift[i]*dat_wing_curr$torso_length > 0.04){next} # overestimation of the error even for the large birds, the absolute error is highly unlikely to exceed +-4cm
      if(dat_bird_curr$x_loc_TorsotailCoG > 0){next} # skip to next if this adjustment would move the torsotail CG in front of the beginning of the torso - non-physical
      if(-dat_bird_curr$x_loc_TorsotailCoG > dat_bird_curr$torsotail_length){next} # skip to next if this adjustment would move the torsotail CG behind the end of the tail - non-physical

      curr_torsotail_data     = tryCatch(massprop_restbody(dat_wing_curr, dat_bird_curr),
                                         error = function(e){err_catch = 0})
      if(length(curr_torsotail_data) == 1){next}

      dat_sens_out$species[count] = species_curr
      dat_sens_out$BirdID[count]  = birdid_curr
      if(l == 1){
        dat_sens_out$group[count]  = "max_agility"
      }else{
        if(l == 2){
          dat_sens_out$group[count]  = "max_stability"
        }else{
          dat_sens_out$group[count]  = "min_stability"
        }
      }
      dat_sens_out$iter[count]   = i

      # Compute the full bird results
      fullbird    = list()
      fullbird$I  = matrix(0, nrow = 3, ncol = 3)
      fullbird$CG = matrix(0, nrow = 3, ncol = 1)
      fullbird$I[1,1]   = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Ixx")]) + 2*dat_wing_curr$wing_Ixx
      fullbird$I[2,2]   = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Iyy")]) + 2*dat_wing_curr$wing_Iyy
      fullbird$I[3,3]   = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Izz")]) + 2*dat_wing_curr$wing_Izz
      fullbird$I[1,3]   = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Ixz")]) + 2*dat_wing_curr$wing_Ixz
      fullbird$I[3,1]   = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Ixz")]) + 2*dat_wing_curr$wing_Ixz

      if(length(subset(curr_torsotail_data, object == "m" & component == "neck")$value) == 0){
        dat_sens_out$CGx[count]   = (subset(curr_torsotail_data, object == "CGx" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                                       subset(curr_torsotail_data, object == "CGx" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                                       subset(curr_torsotail_data, object == "CGx" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                                       2*dat_wing_curr$wing_CGx*dat_wing_curr$wing_m)/dat_wing_curr$full_m}else{
                                         dat_sens_out$CGx[count]   = (subset(curr_torsotail_data, object == "CGx" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                                                                        subset(curr_torsotail_data, object == "CGx" & component == "neck")$value*subset(curr_torsotail_data, object == "m" & component == "neck")$value +
                                                                        subset(curr_torsotail_data, object == "CGx" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                                                                        subset(curr_torsotail_data, object == "CGx" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                                                                        2*dat_wing_curr$wing_CGx*dat_wing_curr$wing_m)/dat_wing_curr$full_m}

      #Shift the I to the CG
      fullbird$I         = parallelaxis(fullbird$I,-c(dat_sens_out$CGx[count],0,dat_wing_curr$full_CGz),dat_wing_curr$full_m,"A")

      # ------ Save all data -------
      dat_sens_out$Iyy[count]            = fullbird$I[2,2]
      dat_sens_out$mean_CGx_specific_orgShoulder[count]  = (dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)/dat_wing_curr$full_length
      dat_sens_out$mean_CGx_specific_orgShoulder_root[count]  = (dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)/dat_wing_curr$c_root

      # ------ Calculate all agility and stability metrics --------
      dat_sens_out$sm[count]             = (dat_sens_out$CGx[count]-dat_wing_curr$pt1_X) - dat_wing_curr$x_np_est_orgShoulder
      dat_sens_out$prop_q_dot[count]     = (-dat_sens_out$sm[count]*dat_wing_curr$S_max*dat_wing_curr$full_m^0.24)/dat_sens_out$Iyy[count]
      dat_sens_out$prop_q_dot_nd[count]  = (-dat_sens_out$sm[count]*dat_wing_curr$S_max*dat_wing_curr$full_length^2)/dat_sens_out$Iyy[count]
      dat_sens_out$sm_nd[count]          = dat_sens_out$sm[count]/dat_wing_curr$c_root_max

      ## ------------ Check that this range includes the possible shift that would be caused by a tail neutral point ---------
      # to estimate downwash in subsonic flow with z_t = 0 - Raymer Fig. 16.12
      dat_sens_out$r_Raymer[count]  = abs(-(dat_wing_curr$torso_length+0.25*dat_wing_curr$tail_length)-dat_wing_curr$full_CGx)/(0.5*dat_wing_curr$span)
      dat_sens_out$AR[count]        = dat_wing_curr$AR
      ## assume that a_t/a*(1-d(epsilon)/d(alpha))*eta_t = 0.73
      dat_sens_out$sm_TailNP[count]       = (dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)-(dat_wing_curr$x_np_est_orgShoulder-(dat_wing_curr$Vh*0.73*dat_wing_curr$c_mean_max))
      dat_sens_out$sm_TailNP_nd[count]    = dat_sens_out$sm_TailNP[count]/dat_wing_curr$c_root_max

      dat_sens_out$binomial[count]        = dat_wing_curr$binomial
      dat_sens_out$elbow[count]           = dat_wing_curr$elbow
      dat_sens_out$manus[count]           = dat_wing_curr$manus
      count = count+1
      remove(curr_torsotail_data)
    }
  }
}
# this line requires the plotting_info.R to be run
dat_sens_out$species_order = factor(dat_sens_out$species, levels = phylo_order$species)
dat_sens_out <- dat_sens_out[c(1:count),]

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_CGsensana.RData",sep="")
save(dat_sens_out,file = filename)
# iter 16 is shift = 0
tmp   <- aggregate(list(min_stab = subset(dat_sens_out, iter == 16)$sm_TailNP_nd),
                   by = list(binomial = subset(dat_sens_out, iter == 16)$binomial,
                             group = subset(dat_sens_out, iter == 16)$group), mean)

data_in        <- as.vector(subset(tmp, group == "max_stability")[,3])
names(data_in) <- subset(tmp, group == "max_stability")$binomial
maxstab_tail_OU_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "OU")
maxstab_tail_BM_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "BM")
maxstab_tail_OU_outputs$opt$aicc-maxstab_tail_BM_outputs$opt$aicc

data_in        <- as.vector(subset(tmp, group == "min_stability")[,3])
names(data_in) <- subset(tmp, group == "min_stability")$binomial
minstab_tail_OU_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "OU")
minstab_tail_BM_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "BM")
minstab_tail_OU_outputs$opt$aicc-minstab_tail_BM_outputs$opt$aicc

### -------------------------------------------------------------------------
### ------------------------ CG Error Bootstrapping -------------------------
### -------------------------------------------------------------------------

### ------------------- Bootstrap results for maximum stability -------------

no_boot = 5000
dat_sens_opt_max = as.data.frame(matrix(nrow = no_boot,ncol = 4))
colnames(dat_sens_opt_max) <- c("z0","alpha","sigsq","err")
for (j in 1:no_boot){
  #must be re-established each time
  dat_sens_opt = as.data.frame(matrix(nrow = 36,ncol = 3))
  colnames(dat_sens_opt) <- c("max_stab","binomial","err")
  for(i in 1:36){
    dat_curr = subset(dat_sens_out, group == "max_stability" & species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])
    iter_curr  = sample(dat_curr$iter, 1, replace = FALSE, prob = NULL)
    dat_sens_opt$max_stab[i] = dat_curr$sm_nd[which(dat_curr$iter == iter_curr)]
    dat_sens_opt$binomial[i] = dat_bird$binomial[i]
    dat_sens_opt$err[i]      = abs(shift[iter_curr])
  }
  # fit the OU model for the data
  tmp   <- aggregate(list(max_stab = dat_sens_opt$max_stab,
                          err = dat_sens_opt$err), by = list(binomial = dat_sens_opt$binomial), mean)
  data_in        <- as.vector(tmp[,2])
  names(data_in) <- tmp$binomial
  stab_OU_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "OU")

  # save the data
  dat_sens_opt_max$z0[j]     = stab_OU_outputs$opt$z0
  dat_sens_opt_max$alpha[j]  = stab_OU_outputs$opt$alpha
  dat_sens_opt_max$sigsq[j]  = stab_OU_outputs$opt$sigsq
  dat_sens_opt_max$err[j]    = mean(tmp$err) # error is the mean shift of the torso CG across all species
}

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_maxstab_z0.RData",sep="")
save(dat_sens_opt_max,file = filename)

### ------------------- Bootstrap results for minimum stability -------------
no_boot = 5000
dat_sens_opt_min = as.data.frame(matrix(nrow = no_boot,ncol = 4))
colnames(dat_sens_opt_min) <- c("z0","alpha","sigsq","err")
for (j in 1:no_boot){
  #must be re-established each time
  dat_sens_opt = as.data.frame(matrix(nrow = 36,ncol = 3))
  colnames(dat_sens_opt) <- c("min_stab","binomial","err")
  for(i in 1:36){
    dat_curr = subset(dat_sens_out, group == "min_stability" & species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])
    iter_curr  = sample(dat_curr$iter, 1, replace = FALSE, prob = NULL)
    dat_sens_opt$min_stab[i] = dat_curr$sm_nd[which(dat_curr$iter == iter_curr)]
    dat_sens_opt$binomial[i] = dat_bird$binomial[i]
    dat_sens_opt$err[i]      = abs(shift[iter_curr])
  }

  # fit the OU model for the data
  tmp   <- aggregate(list(min_stab = dat_sens_opt$min_stab,
                          err = dat_sens_opt$err), by = list(binomial = dat_sens_opt$binomial), mean)
  data_in        <- as.vector(tmp[,2])
  names(data_in) <- tmp$binomial
  stab_OU_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "OU")

  # save the data
  dat_sens_opt_min$z0[j]     = stab_OU_outputs$opt$z0
  dat_sens_opt_min$alpha[j]  = stab_OU_outputs$opt$alpha
  dat_sens_opt_min$sigsq[j]  = stab_OU_outputs$opt$sigsq
  dat_sens_opt_min$err[j]    = mean(tmp$err) # error is the mean shift of the torso CG across all species
}

plot(dat_sens_opt_max$err,dat_sens_opt_max$z0)

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_minstab_z0.RData",sep="")
save(dat_sens_opt_min,file = filename)

##### --------------------------------------------------------------------------------------------
##### ------------------------------------- PMC analysis ----------------------------------------
##### --------------------------------------------------------------------------------------------

# to verify that the OU model is an appropriate fit given the size of our data
out_xcg <- pmc(tree = pruned_mcc, data = all_data_means_mat[,c("mean_CGx_specific_orgShoulder")], "BM", "OU", nboot = 5000, mc.cores = 1)
filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_pmc_xcg_output.RData",sep="")
save(out_xcg,file = filename)
out_maxstab <- pmc(tree = pruned_mcc, data = all_data_means_mat[,c("max_sm_nd")], "BM", "OU", nboot = 5000, mc.cores = 1)
filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_pmc_maxstab_output.RData",sep="")
save(out_maxstab,file = filename)
out_minstab <- pmc(tree = pruned_mcc, data = all_data_means_mat[,c("min_sm_nd")], "BM", "OU", nboot = 5000, mc.cores = 1)
filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_pmc_minstab_output.RData",sep="")
save(out_minstab,file = filename)

# ##### --------------------------------------------------------------------------------------------
# ##### ------------------------------ Supplementary Figures ---------------------------------------
# ##### --------------------------------------------------------------------------------------------

## ----------------------------------- Figure S1 -------------------------------------------------
library("ggplot2")
library("tidyr")
library("dplyr")

load("/Users/christinaharvey/Documents/AvInertia/AnalysisData/2021_08_20_pmc_xcg_output.RData")
dists_x_cg <- data.frame(null = out_xcg$null, test = out_xcg$test)
plot_x_cg <- dists_x_cg %>%
  gather(dists_x_cg, value) %>%
  ggplot(aes(value, fill = dists_x_cg)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out_xcg$lr, linetype = "dashed") +
  scale_fill_manual(values = c("null" = "black","test" = "#FFCC00"), name = "Model fit", labels = c("BM (null)","OU")) +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,1.5), name = "Density", breaks = c(0,0.5,1,1.5)) +
  scale_x_continuous(limits = c(-2,33), name = "Maximum likelihood estimates", breaks = c(0,10,20,30)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 30, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 1.5, x = log(0), xend = log(0))

load("/Users/christinaharvey/Documents/AvInertia/AnalysisData/2021_08_21_pmc_maxstab_output.RData")
dists_maxstab <- data.frame(null = out_maxstab$null, test = out_maxstab$test)
plot_maxstab <- dists_maxstab %>%
  gather(dists_maxstab, value) %>%
  ggplot(aes(value, fill = dists_maxstab)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out_maxstab$lr, linetype = "dashed") +
  scale_fill_manual(values = c("null" = "black","test" = "#FFCC00"), name = "Model fit", labels = c("BM (null)","OU")) +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,1.5), name = "Density", breaks = c(0,0.5,1,1.5)) +
  scale_x_continuous(limits = c(-2,33), name = "Maximum likelihood estimates", breaks = c(0,10,20,30)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 30, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 1.5, x = log(0), xend = log(0))

load("/Users/christinaharvey/Documents/AvInertia/AnalysisData/2021_08_21_pmc_minstab_output.RData")
dists_minstab <- data.frame(null = out_minstab$null, test = out_minstab$test)
plot_minstab <- dists_minstab %>%
  gather(dists_minstab, value) %>%
  ggplot(aes(value, fill = dists_minstab)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out_minstab$lr, linetype = "dashed") +
  scale_fill_manual(values = c("null" = "black","test" = "#FFCC00"), name = "Model fit", labels = c("BM (null)","OU")) +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,1.5), name = "Density", breaks = c(0,0.5,1,1.5)) +
  scale_x_continuous(limits = c(-2,33), name = "Maximum likelihood estimates", breaks = c(0,10,20,30)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 30, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 1.5, x = log(0), xend = log(0))

#exported as 7x5
fig_pmc <- plot_grid(plot_x_cg,plot_maxstab,plot_minstab,
                      #arrangement data
                      ncol = 1,
                      rel_heights = c(1,1,1),
                      #labels
                      labels = c("a","b","c"),
                      label_size = 10,
                      label_fontfamily = "sans")

# calculate power of tests
length(which(out_xcg$test > quantile(out_xcg$null, probs = 0.95)))/5000
length(which(out_maxstab$test > quantile(out_maxstab$null, probs = 0.95)))/5000
length(which(out_minstab$test > quantile(out_minstab$null, probs = 0.95)))/5000

# extract the 95% confidence intervals
dist_alpha_xcg = out_xcg$par_dists$value[which(out_xcg$par_dists$comparison == "BB" & out_xcg$par_dists$parameter =="alpha")]
quantile(dist_alpha, probs = c(0.05,0.95))


#dat_sens_opt_max <- dat_sens_opt_out - run this after loading saved data
1-length(which(dat_sens_opt_max$z0 < 1))/1000
1-length(which(dat_sens_opt_min$z0 > 1))/1000

## ----------------------------------- Figure S2 -------------------------------------------------

load("/Users/christinaharvey/Documents/AvInertia/AnalysisData/2021_08_20_maxstab_z0.RData")
load("/Users/christinaharvey/Documents/AvInertia/AnalysisData/2021_08_20_minstab_z0.RData")

#exported as 3.5x5
boot_stab <- ggplot()+
  # add background info
  geom_density(data = dat_sens_opt_max,aes(x = z0), fill = "#00510A", alpha = 0.8)  +
  geom_density(data = dat_sens_opt_min,aes(x = z0), fill = "#1D0747", alpha = 0.8)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,50), name = "Density", breaks = c(0,25,50)) +
  scale_x_continuous(limits = c(-0.3,0.3), breaks = seq(-0.3,0.3,0.15), name = expression(paste("Phenotypic optimum, ",theta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.3, xend = 0.3, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 50, x = log(0), xend = log(0))

## ----------------------------------- Figure S3 -------------------------------------------------
load("/Users/christinaharvey/Documents/AvInertia/AnalysisData/2021_08_23_CGsensana.RData")

sens_agility_nd <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1.5, xmax = 1.5), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = subset(dat_sens_out, group == "max_agility"), aes(y = species_order, x = prop_q_dot_nd, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16 & group == "max_agility"), aes(y = species_order, x = prop_q_dot_nd, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(limits = c(-1.5,1.5), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5),
                     name = "Normalized maximum stable pitch agility") +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = -1.5, xend = 1.5, y = log(0), yend = log(0))

sens_stability_max <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = subset(dat_sens_out, group == "max_stability"), aes(y = species_order, x = sm_nd, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16 & group == "max_stability"), aes(y = species_order, x = sm_nd, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0, col = "black") +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(limits = c(-0.65,0.65), breaks = seq(-0.6,0.6,0.3), labels = c("-60","-30","0","30","60"),name = "maximum static margin (% of c   )") +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = -0.6, xend = 0.6, y = log(0), yend = log(0))

sens_stability_min <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = subset(dat_sens_out, group == "min_stability"), aes(y = species_order, x = sm_nd, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16 & group == "min_stability"), aes(y = species_order, x = sm_nd, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0, col = "black") +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(limits = c(-0.65,0.65), breaks = seq(-0.6,0.6,0.3), labels = c("-60","-30","0","30","60"),name = "minimum static margin (% of c   )") +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = -0.6, xend = 0.6, y = log(0), yend = log(0))

#exported as 4x12
fig_sens <- plot_grid(sens_agility_nd,sens_stability_min,sens_stability_max,
                     #arrangement data
                     nrow = 1,
                     rel_heights = c(1,1,1),
                     #labels
                     labels = c("a","b","c"),
                     label_size = 10,
                     label_fontfamily = "sans")

