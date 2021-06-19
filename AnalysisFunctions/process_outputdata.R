# This script was written to analyse the data that is returned from the AvInertia package
## Load up necessary packages
library(phytools)
library(ape)
library(geiger)
library(MCMCglmm)
library(tidyverse)
library(effectsize) # needed for eta_squared calculation
library(pracma)
library(ggplot2)
library(dplyr)
library(reshape2)

### ----------------------------------------------------------
### --------------------- Read in data -----------------------
### ----------------------------------------------------------

# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

## ------------- Read in the wing shape data --------------
filename_wing  = list.files(path = path_data_folder, pattern = paste("allspecimen_winginfo"))
dat_wing       = read.csv(file = paste(path_data_folder,filename_wing,sep= ""))
dat_wing       = dat_wing[,-c(1)]
# ---- Calculate the base aerodynamic properties of the wings ----
# don't want the y direction included - just want the distance from leading to trailing edge in x-z plane
# using max chord because we know for a gull that the AC is between 19%-43% of the max chord so this seems like a fair approach
dat_wing$chord        = sqrt((dat_wing$pt12_X - dat_wing$pt11_X)^2+(dat_wing$pt12_Z - dat_wing$pt11_Z)^2)
# span is calculated as the furthest distance of the humerus to either P10 or P7.
# The x component is neglected. Measured from the body center line
dat_wing$span         = apply(cbind(2*sqrt(dat_wing$pt8_Y^2+dat_wing$pt8_Z^2),2*sqrt(dat_wing$pt9_Y^2+dat_wing$pt9_Z^2),2*sqrt(dat_wing$pt7_Y^2+dat_wing$pt7_Z^2)), 1, max)
dat_wing$AR           = (dat_wing$span^2)/(2*dat_wing$S)
dat_wing$AR_proj      = (dat_wing$span^2)/(2*dat_wing$S_proj)

## ------------- Read in the feather data --------------
filename_feat = list.files(path = path_data_folder, pattern = paste("allspecimen_featherinfo"))
dat_feat      = read.csv(file = paste(path_data_folder,filename_feat,sep= ""))
dat_feat      = dat_feat[,-c(1,3:7)]

## ------------- Read in the body morphology data --------------
filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))
names(dat_bird)[names(dat_bird) == "binomial.x"] = "binomial"
dat_bird     = dat_bird[,-c(1)]

# set neck length to zero for those modeled with a retracted neck
for (i in 1:nrow(dat_bird)){
  if(!dat_bird$extend_neck[i]){
    dat_bird$neck_length[i] = 0
  }
}

## ------------- Read in body results from AvInertia --------------
filename_body = list.files(path = path_data_folder, pattern = paste("bodyCGandMOI_correct"))
dat_body      = read.csv(file = paste(path_data_folder,filename_body,sep= ""))
dat_body      = reshape2::dcast(dat_body, species + BirdID + TestID + FrameID ~ component + object, value.var="value")

## ------------- Read in wing configuration results from AvInertia --------------
filename_results = list.files(path = path_data_folder, pattern = paste("results"))
dat_results = read.csv(paste(path_data_folder,filename_results[1],sep=""))
for (i in 2:length(filename_results)){
  dat_results = rbind(dat_results,read.csv(paste(path_data_folder,filename_results[i],sep="")))
}

# clean up environment
remove(filename_wing,filename_feat,filename_bird,filename_body,filename_results)

## ------------- Read in full tree phylogeny --------------
full_tree <- read.nexus("vikROM_passerines_403sp.tre")

### --------------------------------------------------------------
### --------------------- Combine all data -----------------------
### --------------------------------------------------------------

dat_final = merge(dat_results,dat_wing[,c("species","BirdID","TestID","FrameID","BirdID_FrameSpec",
                                          "elbow","manus","S","S_proj","chord","span","AR","AR_proj",
                                          "pt1_X","pt1_Y","pt1_Z")], by = c("species","BirdID","TestID","FrameID"))
dat_final = merge(dat_final,dat_bird, by = c("species","BirdID"))
dat_final = merge(dat_final,dat_body[,-c(3,4)], by = c("species","BirdID"))

### ------------------------------------------------------------------
### --------------------- Calculate key params -----------------------
### ------------------------------------------------------------------

dat_final$full_length  = (dat_final$torsotail_length+dat_final$head_length+dat_final$neck_length)
dat_final$torso_length = dat_final$torsotail_length - dat_final$tail_length

dat_final$full_CGx_orgBeak             = (dat_final$full_CGx-dat_final$head_length-dat_final$neck_length)
dat_final$full_CGx_specific_orgBeak    = (dat_final$full_CGx-dat_final$head_length-dat_final$neck_length)/dat_final$full_length
dat_final$full_CGz_orgDorsal           = (dat_final$full_CGz+dat_final$z_dist_to_veh_ref_point_cm)
dat_final$full_CGz_specific_orgDorsal  = (dat_final$full_CGz+dat_final$z_dist_to_veh_ref_point_cm)/dat_final$full_length
dat_final$shoulderx_specific_orgBeak   = (dat_final$pt1_X-dat_final$head_length-dat_final$neck_length)/dat_final$full_length
dat_final$shoulderz_specific_orgDorsal = (dat_final$pt1_Z+dat_final$z_dist_to_veh_ref_point_cm)/dat_final$full_length

dat_final$full_CGx_specific_orgShoulder = (dat_final$full_CGx-dat_final$pt1_X)/dat_final$full_length
dat_final$full_CGx_orgShoulder          = (dat_final$full_CGx-dat_final$pt1_X)
dat_final$full_CGz_specific_orgShoulder = (dat_final$full_CGz-dat_final$pt1_Z)/dat_final$full_length
dat_final$full_CGz_orgShoulder          = (dat_final$full_CGz-dat_final$pt1_Z)
dat_final$BeakTipx_orgShoulder          = (dat_final$head_length+dat_final$neck_length+dat_final$pt1_X)/dat_final$full_length
dat_final$Centrez_orgShoulder           = (dat_final$pt1_Z)/dat_final$full_length
dat_final$TailTipx_orgShoulder          = (-dat_final$torsotail_length-dat_final$pt1_X)/dat_final$full_length
dat_final$MaxWidthx_orgShoulder         = (-dat_final$x_loc_of_body_max-dat_final$pt1_X)/dat_final$full_length
dat_final$Dorsalz_orgShoulder           = -(dat_final$z_dist_to_veh_ref_point_cm+dat_final$pt1_Z)/dat_final$full_length
dat_final$Ventralz_orgShoulder          = (dat_final$body_height_max-(dat_final$z_dist_to_veh_ref_point_cm+dat_final$pt1_Z))/dat_final$full_length

test <- aggregate(list(humerus_per = dat_final$humerus_length_mm/dat_final$span),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), min)

dat_final$elbow_scaled = dat_final$elbow*0.001
dat_final$manus_scaled = dat_final$manus*0.001

### ---------------------------------------------------------------------------------------
### --------------------------- Compute summed quantities ---------------------------------
### ---------------------------------------------------------------------------------------
# Maximum projected wing area and maximum wing span
test       <- aggregate(list(S_proj_max = dat_final$wing_S_proj,
                             S_max = dat_final$S,
                             b_max = dat_final$span,
                             c_max = dat_final$chord),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird   <- merge(dat_bird,test, by = c("species","BirdID"))
dat_final  <- merge(dat_final,test, by = c("species","BirdID"))

dat_final$wing_CGy_orgShoulder           = (dat_final$wing_CGy-dat_final$pt1_Y)
dat_final$wing_CGy_specific              = (dat_final$wing_CGy)/(0.5*dat_final$b_max)
dat_final$wing_CGy_specific_orgShoulder  = (dat_final$wing_CGy-dat_final$pt1_Y)/(0.5*dat_final$b_max)
dat_final$wing_CGx_specific_orgBeak      = (dat_final$wing_CGx-dat_final$head_length-dat_final$neck_length)/(0.5*dat_final$b_max)
dat_final$wing_CGx_specific_orgShoulder  = (dat_final$wing_CGx-dat_final$pt1_X)/(0.5*dat_final$b_max)
dat_final$wing_CGz_specific              = dat_final$wing_CGz/(0.5*dat_final$b_max)

dat_final$full_Ixx_specific  = dat_final$full_Ixx/(dat_final$full_m*dat_final$b_max*dat_final$full_length)
dat_final$full_Iyy_specific  = dat_final$full_Iyy/(dat_final$full_m*dat_final$b_max*dat_final$full_length)
dat_final$full_Izz_specific  = dat_final$full_Izz/(dat_final$full_m*dat_final$b_max*dat_final$full_length)
dat_final$full_Ixz_specific  = dat_final$full_Ixz/(dat_final$full_m*dat_final$b_max*dat_final$full_length)

# Shoulder position specific and mass
test     <- aggregate(list(full_m = dat_final$full_m),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))

### ---------------------------------------------------------------------------------------
### ------------- Compute extremes of the CG position due to shoulder motion --------------
### ---------------------------------------------------------------------------------------

angle_shoulder = 90
# the most effect will be felt for the wing that has the most distal CG position (independent from it's x or z position)
dat_final$shoulderCG_dist <- sqrt((dat_final$wing_CGy-dat_final$pt1_Y)^2)
tmp = aggregate(list(max_wingCG = dat_final$shoulderCG_dist),
                       by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
shoulder_motion = dat_final[which(dat_final$shoulderCG_dist %in% tmp$max_wingCG),c("species","BirdID","elbow","manus","wing_CGx","wing_CGy","wing_CGz","pt1_X","pt1_Y","pt1_Z","b_max",
                                                                             "wing_m","full_m","full_CGx","full_CGy","full_CGz","full_length","head_length","neck_length","binomial")]
shoulder_motion$rest_m   = (shoulder_motion$full_m-2*shoulder_motion$wing_m)
shoulder_motion$rest_CGx = (shoulder_motion$full_m*shoulder_motion$full_CGx - 2*(shoulder_motion$wing_m*shoulder_motion$wing_CGx))/shoulder_motion$rest_m
shoulder_motion$rest_CGz = (shoulder_motion$full_m*shoulder_motion$full_CGz - 2*(shoulder_motion$wing_m*shoulder_motion$wing_CGz))/shoulder_motion$rest_m

# - Rotate the wing forward about the shoulder in the x-y plane - allows a rotation about Pt1
new_wing_CGx = (cosd(angle_shoulder)*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + sind(angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_X
new_wing_CGy = (-sind(angle_shoulder)*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + cosd(angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
shoulder_motion$forward_CGx = (shoulder_motion$rest_m*shoulder_motion$rest_CGx + 2*shoulder_motion$wing_m*new_wing_CGx)/shoulder_motion$full_m

# - Rotate the wing backwards about the shoulder in the x-y plane
new_wing_CGx = (cosd(-angle_shoulder)*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + sind(-angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_X
new_wing_CGy = (-sind(-angle_shoulder)*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + cosd(-angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
shoulder_motion$backwards_CGx = (shoulder_motion$rest_m*shoulder_motion$rest_CGx + 2*shoulder_motion$wing_m*new_wing_CGx)/shoulder_motion$full_m

shoulder_motion$forward_CGx_specific   = (shoulder_motion$forward_CGx-shoulder_motion$pt1_X)/shoulder_motion$full_length
shoulder_motion$backwards_CGx_specific = (shoulder_motion$backwards_CGx-shoulder_motion$pt1_X)/shoulder_motion$full_length
shoulder_motion$range_CGx              = (shoulder_motion$forward_CGx-shoulder_motion$backwards_CGx)
shoulder_motion$range_CGx_specific     = shoulder_motion$range_CGx/(shoulder_motion$full_length)

# - Rotate the wing up about the shoulder in the y-z plane
new_wing_CGz = (cosd(angle_shoulder)*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z) - sind(angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Z
new_wing_CGy = (sind(angle_shoulder)*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z) + cosd(angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
shoulder_motion$upwards_CGz   = (shoulder_motion$rest_m*shoulder_motion$rest_CGz + 2*shoulder_motion$wing_m*new_wing_CGz)/shoulder_motion$full_m

# - Rotate the wing down about the shoulder in the y-z plane
new_wing_CGz = (cosd(-angle_shoulder)*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z) - sind(-angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Z
new_wing_CGy = (sind(-angle_shoulder)*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z) + cosd(-angle_shoulder)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
shoulder_motion$downwards_CGz = (shoulder_motion$rest_m*shoulder_motion$rest_CGz + 2*shoulder_motion$wing_m*new_wing_CGz)/shoulder_motion$full_m

shoulder_motion$upwards_CGz_specific   = (shoulder_motion$upwards_CGz-shoulder_motion$pt1_Z)/shoulder_motion$full_length
shoulder_motion$downwards_CGz_specific = (shoulder_motion$downwards_CGz-shoulder_motion$pt1_Z)/shoulder_motion$full_length
shoulder_motion$range_CGz              = (shoulder_motion$downwards_CGz-shoulder_motion$upwards_CGz)
shoulder_motion$range_CGz_specific     = shoulder_motion$range_CGz/(shoulder_motion$full_length)
shoulder_motion$span_ratio             = shoulder_motion$b_max/shoulder_motion$full_length
shoulder_motion                        = mutate(shoulder_motion, phylo = binomial)
### ---------------------- Compute inertial metrics -------------------------

# uses scaling from:
# Alerstam, T., Rosén, M., Bäckman, J., Ericson, P. G., & Hellgren, O. (2007).
# Flight speeds among bird species: allometric and phylogenetic effects. PLoS Biol, 5(8), e197.
dat_final$prop_q_dot     <- (abs((dat_final$full_CGx-dat_final$pt1_X)-0.25*dat_final$chord)*dat_final$S_max*dat_final$full_m^0.24)/dat_final$full_Iyy
dat_final$prop_q_dot_nd  <- (abs((dat_final$full_CGx-dat_final$pt1_X)-0.25*dat_final$chord)*dat_final$S_max*dat_final$full_length^2)/dat_final$full_Iyy
dat_final$del_M_specific <- dat_final$prop_q_dot*dat_final$full_Iyy/(dat_final$full_m*dat_final$full_length)

dat_final$sachs_pred_Ixx <- dat_final$full_m*(sqrt((0.14*dat_final$wing_m/dat_final$full_m))*dat_final$span*0.5)^2

dat_final$pitch_div <- (dat_final$full_Izz - dat_final$full_Ixx)/dat_final$full_Iyy
dat_final$yaw_div   <- (dat_final$full_Iyy - dat_final$full_Ixx)/dat_final$full_Izz
## ---- Compute the regression coefficients for each species for each variable -------

no_specimens <- nrow(dat_bird)
dat_model_out        <- data.frame(matrix(nrow = no_specimens*9*2, ncol = 5))
names(dat_model_out) <- c("species","model_variable","elb","man","elbman")
varlist_sp       <- c("full_Ixx_specific","full_Iyy_specific","full_Izz_specific","full_Ixz_specific",
                   "full_CGx_specific_orgBeak","full_CGx_specific_orgShoulder","wing_CGy_specific_orgShoulder","full_CGz_specific_orgDorsal","wing_CGx_specific_orgBeak")
short_varlist_sp <- c("Ixx_sp","Iyy_sp","Izz_sp","Ixz_sp","CGx_sp","CGx_sp_sh","CGy_sp","CGz_sp","CGx_wing_sp")
varlist_abs      <- c("full_Ixx","full_Iyy","full_Izz","full_Ixz",
                      "full_CGx_orgBeak","full_CGx_orgShoulder","wing_CGy_orgShoulder","full_CGz","wing_CGx")
short_varlist_abs<- c("Ixx","Iyy","Izz","Ixz","CGx","CGx_sh","CGy","CGz", "CGx_wing")
dat_bird$species <- as.character(dat_bird$species)
success = TRUE
count = 1
for (i in 1:no_specimens){
  # subset data to the current species
  tmp = subset(dat_final, species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])
  tmp$elbow_centered = tmp$elbow-mean(tmp$elbow)
  tmp$manus_centered = tmp$manus-mean(tmp$manus)
  for (m in 1:2){

    if(m == 1){
      varlist = varlist_abs
      short_varlist = short_varlist_abs
    } else{
      varlist = varlist_sp
      short_varlist = short_varlist_sp
    }

    # --------------- Fit models ------------------
    # full models
    models <- lapply(varlist, function(x) {lm(substitute(k ~ elbow_scaled*manus_scaled, list(k = as.name(x))), data = tmp)})
    # compute the effect size from the centered values of elbow and wrist to avoid conflating the differences between species in absolute range
    models_adj <- lapply(varlist, function(x) {lm(substitute(k*100 ~ elbow_centered*manus_centered, list(k = as.name(x))), data = tmp)})

    CI_values      = lapply(models, confint)
    coef_values    = lapply(models, coef)
    # Compute the effect size
    etap_values    = lapply(models_adj, function(x){eta_squared(car::Anova(x, type = 3))}) # partial = TRUE by default
    count2 = 1
    for (j in 1:length(models)){
      dat_model_out$species[count]        <- dat_bird$species[i]
      dat_model_out$BirdID[count]         <- dat_bird$BirdID[i]
      dat_model_out$model_variable[count] <- short_varlist[j]
      dat_model_out$int[count]            <- coef_values[[j]]["(Intercept)"]
      dat_model_out$elb[count]            <- coef_values[[j]]["elbow_scaled"]
      dat_model_out$man[count]            <- coef_values[[j]]["manus_scaled"]
      dat_model_out$r2[count]             <- summary(models[[j]])$r.squared
      dat_model_out$elb_p[count]          <- summary(models[[j]])$coefficients["elbow_scaled",4]
      dat_model_out$man_p[count]          <- summary(models[[j]])$coefficients["manus_scaled",4]
      dat_model_out$elbman_p[count]       <- summary(models[[j]])$coefficients["elbow_scaled:manus_scaled",4]
      dat_model_out$elbman[count]         <- coef_values[[j]]["elbow_scaled:manus_scaled"]
      dat_model_out$int_lb[count]         <- CI_values[[j]]["(Intercept)",1]
      dat_model_out$elb_lb[count]         <- CI_values[[j]]["elbow_scaled",1]
      dat_model_out$man_lb[count]         <- CI_values[[j]]["manus_scaled",1]
      dat_model_out$elbman_lb[count]      <- CI_values[[j]]["elbow_scaled:manus_scaled",1]
      dat_model_out$int_ub[count]         <- CI_values[[j]]["(Intercept)",2]
      dat_model_out$elb_ub[count]         <- CI_values[[j]]["elbow_scaled",2]
      dat_model_out$man_ub[count]         <- CI_values[[j]]["manus_scaled",2]
      dat_model_out$elbman_ub[count]      <- CI_values[[j]]["elbow_scaled:manus_scaled",2]
      # Compute the effect sizes
      dat_model_out$elb_etap[count]        <- etap_values[[j]][1,2]
      dat_model_out$man_etap[count]        <- etap_values[[j]][2,2]
      dat_model_out$elbman_etap[count]     <- etap_values[[j]][3,2]
      count = count + 1
    }
  }
}

tmp           = reshape2::melt(dat_model_out, id = c("species","BirdID","model_variable"))
dat_model_out = reshape2::dcast(tmp, species + BirdID ~ model_variable + variable, value.var="value")

# Include basic geometry effects
tmp       = aggregate(list(full_m = dat_bird$full_m),  by=list(species = dat_bird$species, binomial = dat_bird$binomial, BirdID = dat_bird$BirdID), mean)
dat_comp  = merge(tmp, dat_model_out, by = c("species","BirdID"))
dat_model_out  = merge(tmp, dat_model_out, by = c("species","BirdID"))
tmp       = aggregate(list(torsotail_length = dat_bird$torsotail_length),  by=list(species = dat_bird$species, BirdID = dat_bird$BirdID), mean)
dat_comp  = merge(dat_comp, tmp, by = c("species","BirdID"))

# Include other important factors
# Range of each component - doesn't matter where the origin is for the range
test     <- aggregate(list(range_CGx               = dat_final$full_CGx,
                           range_CGx_specific      = dat_final$full_CGx_specific_orgBeak,
                           range_wing_CGy          = dat_final$wing_CGy,
                           range_wing_CGy_specific = dat_final$wing_CGy_specific,
                           range_CGz               = dat_final$full_CGz,
                           range_CGz_specific      = dat_final$full_CGz_specific_orgDorsal,
                           range_Ixx               = dat_final$full_Ixx,
                           range_Ixx_specific      = dat_final$full_Ixx_specific,
                           range_Iyy               = dat_final$full_Iyy,
                           range_Iyy_specific      = dat_final$full_Iyy_specific,
                           range_Izz               = dat_final$full_Izz,
                           range_Izz_specific      = dat_final$full_Izz_specific),
                      by=list(species = dat_final$species, BirdID = dat_final$BirdID), FUN=function(x){max(x)-min(x)})
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))

# Include other important factors
# Range of each component
test     <- aggregate(list(mean_CGx_orgBeak      = dat_final$full_CGx_orgBeak,
                           mean_CGz_orgDorsal    = dat_final$full_CGz_orgDorsal,
                           mean_CGx_orgShoulder  = dat_final$full_CGx_orgShoulder,
                           mean_CGz_orgShoulder  = dat_final$full_CGz_orgShoulder,
                           mean_wing_CGy                 = dat_final$wing_CGy_orgShoulder,
                           mean_CGx_specific_orgBeak     = dat_final$full_CGx_specific_orgBeak,
                           mean_CGz_specific_orgDorsal   = dat_final$full_CGz_specific_orgDorsal,
                           mean_CGx_specific_orgShoulder = dat_final$full_CGx_specific_orgShoulder,
                           mean_CGz_specific_orgShoulder = dat_final$full_CGz_specific_orgShoulder,
                           mean_wing_CGy_specific        = dat_final$wing_CGy_specific_orgShoulder,
                           mean_Ixx_specific = dat_final$full_Ixx_specific,
                           mean_Iyy_specific = dat_final$full_Iyy_specific,
                           mean_Izz_specific = dat_final$full_Izz_specific,
                           mean_Ixz_specific = dat_final$full_Ixz_specific,
                           mean_Ixx          = dat_final$full_Ixx,
                           mean_Iyy          = dat_final$full_Iyy,
                           mean_Izz          = dat_final$full_Izz,
                           mean_Ixz          = dat_final$full_Ixz,
                           full_length = dat_final$full_length),
                      by=list(species = dat_final$species, BirdID = dat_final$BirdID), mean)
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))
# Maximum values
test     <- aggregate(list(max_CGx_orgBeak       = dat_final$full_CGx_orgBeak,
                           max_CGz_orgDorsal     = dat_final$full_CGz_orgDorsal,
                           max_CGx_specific      = dat_final$full_CGx_specific_orgBeak,
                           max_CGx_orgShoulder   = dat_final$full_CGx_orgShoulder,
                           max_wing_CGy          = dat_final$wing_CGy_orgShoulder,
                           max_wing_CGy_specific = dat_final$wing_CGy_specific_orgShoulder,
                           max_CGz_specific      = dat_final$full_CGz_specific_orgDorsal,
                           max_Ixx               = dat_final$full_Ixx,
                           max_wing_Ixx          = dat_final$wing_Ixx,
                           sachs_pred_Ixx        = dat_final$sachs_pred_Ixx,
                           max_Ixx_specific      = dat_final$full_Ixx_specific,
                           max_Iyy               = dat_final$full_Iyy,
                           max_Iyy_specific      = dat_final$full_Iyy_specific,
                           max_Izz               = dat_final$full_Izz,
                           max_Izz_specific      = dat_final$full_Izz_specific,
                           max_Ixz_specific      = dat_final$full_Ixz_specific,
                           max_q                 = dat_final$prop_q_dot,
                           max_q_nd              = dat_final$prop_q_dot_nd,
                           max_wingspan          = dat_final$span,
                           max_length            = dat_final$full_length,
                           max_S                 = dat_final$S_max,
                           max_S_proj            = dat_final$S_proj_max,
                           max_stab              = 0.25*dat_final$chord/-dat_final$full_CGx_orgShoulder),
                      by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))
# Minimum values
test     <- aggregate(list(min_CGx_orgBeak       = dat_final$full_CGx_orgBeak,
                           min_CGz_orgDorsal     = dat_final$full_CGz_orgDorsal,
                           min_CGx_specific      = dat_final$full_CGx_specific_orgBeak,
                           min_CGx_orgShoulder   = dat_final$full_CGx_orgShoulder,
                           min_wing_CGy          = dat_final$wing_CGy_orgShoulder,
                           min_wing_CGy_specific = dat_final$wing_CGy_specific_orgShoulder,
                           min_CGz_specific      = dat_final$full_CGz_specific_orgDorsal,
                           min_wing_Ixx          = dat_final$wing_Ixx,
                           min_Ixx               = dat_final$full_Ixx,
                           min_Ixx_specific      = dat_final$full_Ixx_specific,
                           min_Iyy               = dat_final$full_Iyy,
                           min_Iyy_specific      = dat_final$full_Iyy_specific,
                           min_Izz               = dat_final$full_Izz,
                           min_q                 = dat_final$prop_q_dot,
                           min_q_nd              = dat_final$prop_q_dot_nd,
                           min_wingspan          = dat_final$span,
                           min_Izz_specific      = dat_final$full_Izz_specific,
                           min_Ixz_specific      = dat_final$full_Ixz_specific,
                           hum_len               = (dat_final$humerus_length_mm+dat_final$ulna_length_mm+dat_final$radius_length_mm+dat_final$cmc_length_mm),
                           min_stab              = 0.25*dat_final$chord/-dat_final$full_CGx_orgShoulder),
                      by=list(species = dat_final$species, BirdID = dat_final$BirdID), min)
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))

dat_comp$max_wing_CGx <- dat_final$wing_CGx[which((dat_final$wing_CGy - dat_final$pt1_Y) %in% dat_comp$max_wing_CGy)]
dat_comp$max_wing_CGx_specific <- dat_final$wing_CGx_specific_orgShoulder[which((dat_final$wing_CGy - dat_final$pt1_Y) %in% dat_comp$max_wing_CGy)]

# --------- Trim the trees -----------
# critical for PGLS models
dat_comp <- mutate(dat_comp, phylo = binomial)
## Prune down the tree to the relevant species
sp_mean_matched <- keep.tip(phy = full_tree, tip = dat_comp$binomial)
## ladderization rotates nodes to make it easier to see basal vs derived
pruned_mcc      <- ape::ladderize(sp_mean_matched)
# plot plot(pruned_mcc)
## The phylogeny will need to be re-formatted for use within MCMCglmm
inv.phylo <- inverseA(pruned_mcc, nodes = "TIPS", scale = TRUE)
## This is the heirarcy of the univariate prior.
univ_prior <-
  list(G = list(G1 = list(V = 1, nu = 0.02)),
       R = list(V = 1, nu = 0.02))




# to compute Pagels lambda - (pgls_model_mcmc$VCV[, 1] / (pgls_model_mcmc$VCV[, 1] + pgls_model_mcmc$VCV[, 2])) %>% mean


filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_alldata.csv",sep="")
write.csv(dat_final,filename)
filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_compdata.csv",sep="")
write.csv(dat_comp,filename)

