
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"
# write in a save function into process_data that pulls out the maximum agility wing configs for each individual
# read in here
devtools::load_all()

dat_sens_max_ag   = dat_final[which(dat_final$prop_q_dot_nd %in% dat_comp$max_q_nd),]
dat_sens_max_stab = dat_final[which((0.25*dat_final$chord/-dat_final$full_CGx_orgShoulder) %in% dat_comp$max_stab),]
dat_sens_min_stab = dat_final[which((0.25*dat_final$chord/-dat_final$full_CGx_orgShoulder) %in% dat_comp$min_stab),]

shift <- seq(-0.15,0.15, by = 0.01)
count = 1
dat_sens_out = as.data.frame(matrix(0, nrow = 10000, ncol = 0))

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
      dat_sens_out$Iyy[count]            <- fullbird$I[2,2]
      dat_sens_out$prop_q_dot[count]     <- (abs((dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)-0.25*dat_wing_curr$chord)*dat_wing_curr$S_max*dat_wing_curr$full_m^0.24)/dat_sens_out$Iyy[count]
      dat_sens_out$prop_q_dot_nd[count]  <- (abs((dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)-0.25*dat_wing_curr$chord)*dat_wing_curr$S_max*dat_wing_curr$full_length^2)/dat_sens_out$Iyy[count]
      dat_sens_out$ac_range[count]       <- 0.25*dat_wing_curr$chord/-(dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)
      dat_sens_out$mean_CGx_specific_orgShoulder[count]  <- dat_sens_out$CGx[count]-dat_wing_curr$pt1_X
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
### ------------------- Bootstrap results for maximum stability -------------

nboot = 1000
dat_sens_opt_out = as.data.frame(matrix(nrow = nboot,ncol = 2))
colnames(dat_sens_opt_out) <- c("z0","err")
for (j in 1:nboot){
  #must be re-established each time
  dat_sens_opt = as.data.frame(matrix(nrow = 36,ncol = 3))
  colnames(dat_sens_opt) <- c("max_stab","binomial","err")
  for(i in 1:36){
    dat_curr = subset(dat_sens_out, group == "max_stability" & species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])
    iter_curr  = sample(dat_curr$iter, 1, replace = FALSE, prob = NULL)
    dat_sens_opt$max_stab[i] = dat_curr$ac_range[which(dat_curr$iter == iter_curr)]
    dat_sens_opt$binomial[i] = dat_bird$binomial[i]
    dat_sens_opt$err[i]      = abs(shift[iter_curr])
  }
  tmp   <- aggregate(list(max_stab = dat_sens_opt$max_stab,
                          err = dat_sens_opt$err), by = list(binomial = dat_sens_opt$binomial), mean)
  data_in        <- as.vector(log(tmp[,2]))
  names(data_in) <- tmp$binomial

  stab_OU_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "OU")
  dat_sens_opt_out$z0[j]  = exp(stab_OU_outputs$opt$z0)
  dat_sens_opt_out$err[j] = sum(tmp$err)
}

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_maxstab_z0.RData",sep="")
save(dat_sens_opt_out,file = filename)

### ------------------- Bootstrap results for minimum stability -------------
nboot = 1000
dat_sens_opt_out = as.data.frame(matrix(nrow = nboot,ncol = 2))
colnames(dat_sens_opt_out) <- c("z0","err")
for (j in 1:nboot){
  #must be re-established each time
  dat_sens_opt = as.data.frame(matrix(nrow = 36,ncol = 3))
  colnames(dat_sens_opt) <- c("min_stab","binomial","err")
  for(i in 1:36){
    dat_curr = subset(dat_sens_out, group == "min_stability" & species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])
    iter_curr  = sample(dat_curr$iter, 1, replace = FALSE, prob = NULL)
    dat_sens_opt$min_stab[i] = dat_curr$ac_range[which(dat_curr$iter == iter_curr)]
    dat_sens_opt$binomial[i] = dat_bird$binomial[i]
    dat_sens_opt$err[i]      = abs(shift[iter_curr])
  }
  tmp   <- aggregate(list(min_stab = dat_sens_opt$min_stab,
                          err = dat_sens_opt$err), by = list(binomial = dat_sens_opt$binomial), mean)
  data_in        <- as.vector(log(tmp[,2]))
  names(data_in) <- tmp$binomial

  stab_OU_outputs = fitContinuous(phy = pruned_mcc, dat = data_in, model = "OU")
  dat_sens_opt_out$z0[j]  = exp(stab_OU_outputs$opt$z0)
  dat_sens_opt_out$err[j] = sum(tmp$err)
}

plot(dat_sens_opt_max$err,dat_sens_opt_max$z0)

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_minstab_z0.RData",sep="")
save(dat_sens_opt_out,file = filename)

##### --------------------------------------------------------------------------------------------
##### ------------------------------ Supplementary Figures ---------------------------------------
##### --------------------------------------------------------------------------------------------

## ----------------------------------- Figure S1 -------------------------------------------------
library("ggplot2")
library("tidyr")
library("dplyr")
dists_x_cg <- data.frame(null = out_xcg$null, test = out_xcg$test)
plot_x_cg <- dists_x_cg %>%
  gather(dists_x_cg, value) %>%
  ggplot(aes(value, fill = dists_x_cg)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out_xcg$lr) +
  scale_fill_manual(values = c("null" = "black","test" = "gray70"), name = "Model fit", labels = c("BM (null)","OU")) +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,1.5), name = "Density", breaks = c(0,0.5,1,1.5)) +
  scale_x_continuous(limits = c(-1,30), name = "Maximum likelihood estimates", breaks = c(0,10,20,30)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 30, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 1.5, x = log(0), xend = log(0))

dists_maxstab <- data.frame(null = out_maxstab$null, test = out_maxstab$test)
plot_maxstab <- dists_maxstab %>%
  gather(dists_maxstab, value) %>%
  ggplot(aes(value, fill = dists_maxstab)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out_xcg$lr) +
  scale_fill_manual(values = c("null" = "black","test" = "gray70"), name = "Model fit", labels = c("BM (null)","OU")) +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,1), name = "Density", breaks = c(0,0.5,1)) +
  scale_x_continuous(limits = c(-1,30), name = "Maximum likelihood estimates", breaks = c(0,10,20,30)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 30, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 1, x = log(0), xend = log(0))

dists_minstab <- data.frame(null = out_minstab$null, test = out_minstab$test)
plot_minstab <- dists_minstab %>%
  gather(dists_minstab, value) %>%
  ggplot(aes(value, fill = dists_minstab)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = out_xcg$lr) +
  scale_fill_manual(values = c("null" = "black","test" = "gray70"), name = "Model fit", labels = c("BM (null)","OU")) +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,1), name = "Density", breaks = c(0,0.5,1)) +
  scale_x_continuous(limits = c(-1,30), name = "Maximum likelihood estimates", breaks = c(0,10,20,30)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 30, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 1, x = log(0), xend = log(0))

#dat_sens_opt_max <- dat_sens_opt_out - run this after loading saved data
1-length(which(dat_sens_opt_max$z0 < 1))/1000
1-length(which(dat_sens_opt_min$z0 > 1))/1000

#exported as 7x5
fig_pmc <- plot_grid(plot_x_cg,plot_maxstab,plot_minstab,
                      #arrangement data
                      ncol = 1,
                      rel_heights = c(1,1,1),
                      #labels
                      labels = c("a","b","c"),
                      label_size = 10,
                      label_fontfamily = "sans")
#exported as 3.5x5
boot_stab <- ggplot()+
  # add background info
  geom_density(data = dat_sens_opt_max,aes(x = z0), fill = "#00510A", alpha = 0.5)  +
  geom_density(data = dat_sens_opt_min,aes(x = z0), fill = "#1D0747", alpha = 0.5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # theme control
  th +
  # axis control
  scale_y_continuous(limits = c(0,12.5), name = "Density", breaks = c(0,4,8,12)) +
  scale_x_continuous(limits = c(0.5,1.5), name = expression(paste("Phenotypic optimum, ",theta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0.5, xend = 1.5, y = log(0), yend = log(0))+
  annotate(geom = "segment", y = 0, yend = 12, x = log(0), xend = log(0))



## ----------------------------------- Figure S3 -------------------------------------------------

sens_agility <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0.1, xmax = 20), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = subset(dat_sens_out, group == "max_agility"), aes(y = species_order, x = prop_q_dot, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16 & group == "max_agility"), aes(y = species_order, x = prop_q_dot, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(trans = "log10", limits = c(0.1,20), breaks = c(0.1,1,10),
                     labels = c(expression(10^-1),expression(10^0),expression(10^1)), name = expression(paste("Maximum pitch agility (s"^{-2},")"))) +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = 0.1, xend = 10, y = log(0), yend = log(0))

sens_agility_nd <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0.1, xmax = 1.3), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = subset(dat_sens_out, group == "max_agility"), aes(y = species_order, x = prop_q_dot_nd, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16 & group == "max_agility"), aes(y = species_order, x = prop_q_dot_nd, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(trans = "log10", limits = c(0.1,1.3), breaks = c(0.1,1,10),
                     labels = c(expression(10^-1),expression(10^0),expression(10^1)),  name = expression(paste("Normalized maximum pitch agility (s"^{-2},")"))) +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = 0.1, xend = 1, y = log(0), yend = log(0))

sens_stability <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = subset(dat_sens_out, group == "max_stability"), aes(y = species_order, x = log(ac_range), col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16 & group == "max_stability"), aes(y = species_order, x = log(ac_range), group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0, col = "black") +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(limits = c(-1,2), breaks = c(-1,0,1,2),
                     labels = c(-1,0,expression(10^1),2), name = expression(paste("log(pitch stability ratio) (s"^{-2},")"))) +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))

#exported as 4x12
fig_sens <- plot_grid(sens_agility,sens_agility_nd,sens_stability,
                     #arrangement data
                     nrow = 1,
                     rel_heights = c(1,1,1),
                     #labels
                     labels = c("a","b","c"),
                     label_size = 10,
                     label_fontfamily = "sans")
