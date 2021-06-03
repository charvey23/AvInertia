
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"
# write in a save function into process_data that pulls out the maximum agility wing configs for each individual
# read in here
dat_sens = read.csv(file = paste(path_data_folder,"wings_sensitivityanalysis.csv",sep= ""))

# also read in the output dat_bird_all from the run_script
filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))

shift <- seq(-0.15,0.15, by = 0.01)
count = 1
dat_sens_out = as.data.frame(matrix(nrow = 1476, ncol = 0))

for(m in 1:nrow(dat_sens)){
  species_curr = dat_bird$species[m]
  birdid_curr  = dat_bird$BirdID[m]

  dat_bird_curr = subset(dat_bird, species == species_curr & BirdID == birdid_curr)
  dat_wing_curr = subset(dat_sens, species == species_curr & BirdID == birdid_curr)
  true_CGx      = dat_bird_curr$x_loc_TorsotailCoG

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
    dat_sens_out$iter[count]   = i
    dat_sens_out$Iyy[count]   = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Iyy")]) + 2*dat_wing_curr$wing_Iyy
    if(length(subset(curr_torsotail_data, object == "m" & component == "neck")$value) == 0){
      dat_sens_out$CGx[count]   = (subset(curr_torsotail_data, object == "CGx" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                                     subset(curr_torsotail_data, object == "CGx" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                                     subset(curr_torsotail_data, object == "CGx" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                                     2*dat_wing_curr$wing_CGx*dat_wing_curr$wing_m)/dat_wing_curr$full_m}
    else{
      dat_sens_out$CGx[count]   = (subset(curr_torsotail_data, object == "CGx" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                                     subset(curr_torsotail_data, object == "CGx" & component == "neck")$value*subset(curr_torsotail_data, object == "m" & component == "neck")$value +
                                     subset(curr_torsotail_data, object == "CGx" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                                     subset(curr_torsotail_data, object == "CGx" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                                     2*dat_wing_curr$wing_CGx*dat_wing_curr$wing_m)/dat_wing_curr$full_m}
    dat_sens_out$prop_q_dot[count]     <- (abs((dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)-0.25*dat_wing_curr$chord)*dat_wing_curr$S_max*dat_wing_curr$full_m^0.24)/dat_sens_out$Iyy[count]
    dat_sens_out$prop_q_dot_nd[count]  <- (abs((dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)-0.25*dat_wing_curr$chord)*dat_wing_curr$S_max*dat_wing_curr$full_length^2)/dat_sens_out$Iyy[count]
    dat_sens_out$ac_range[count]       <- 0.25*dat_wing_curr$chord/-(dat_sens_out$CGx[count]-dat_wing_curr$pt1_X)

    count = count+1
  }
}
# this line requires the plotting_info.R to be run
dat_sens_out$species_order = factor(dat_sens_out$species, levels = phylo_order$species)


sens_agility <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0.1, xmax = 10), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = dat_sens_out, aes(y = species_order, x = prop_q_dot, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16), aes(y = species_order, x = prop_q_dot, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(trans = "log10", limits = c(0.1,10), breaks = c(0.1,1,10),
                     labels = c(expression(10^-1),expression(10^0),expression(10^1)), name = expression(paste("Maximum pitch agility (s"^{-2},")"))) +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = 0.1, xend = 10, y = log(0), yend = log(0))

sens_agility_nd <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0.06, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = dat_sens_out, aes(y = species_order, x = prop_q_dot_nd, col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16), aes(y = species_order, x = prop_q_dot_nd, group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(trans = "log10", limits = c(0.06,1), breaks = c(0.1,1,10),
                     labels = c(expression(10^-1),expression(10^0),expression(10^1)),  name = expression(paste("Normalized maximum pitch agility (s"^{-2},")"))) +
  scale_y_discrete(expand = c(0.1,0), limits = rev(phylo_order$species), name = "") +
  annotate(geom = "segment", x = 0.1, xend = 1, y = log(0), yend = log(0))

sens_stability <- ggplot()+
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_point(data = dat_sens_out, aes(y = species_order, x = log(ac_range), col = shift[iter], group = as.factor(BirdID)),
             position=position_dodge(width=0.5)) +
  geom_point(data = subset(dat_sens_out,iter==16), aes(y = species_order, x = log(ac_range), group = as.factor(BirdID)),
             fill = "white", col = "black", size = 1.5, pch = 23, position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0, col = "black") +
  scale_colour_gradient2(low = "#D01B1B", mid = "white",high = "#95D2EC",midpoint = 0, name = "CG shift (%)", breaks = c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), labels = c(-15,-10,-5,0,5,10,15)) +
  th +
  scale_x_continuous(limits = c(-1,2), breaks = c(-1,0,1,2),
                     labels = c(-1,0,expression(10^1),2), name = expression(paste("log(static margin ratio) (s"^{-2},")"))) +
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

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_sensitivitydata.csv",sep="")
write.csv(dat_sens,filename)
