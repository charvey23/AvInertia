# This script was written to analyse the data that is returned from the birdmoment package
## Load up necessary packages
library(phytools)
library(ape)
library(geiger)
library(MCMCglmm)
library(tidyverse)

library(pracma)
library(ggplot2)
library(dplyr)


# --------------------- Read in data -----------------------
# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

filename_wing = list.files(path = path_data_folder, pattern = paste("allspecimen_winginfo"))
dat_wing      = read.csv(file = paste(path_data_folder,filename_wing,sep= ""))

filename_feat = list.files(path = path_data_folder, pattern = paste("allspecimen_featherinfo"))
dat_feat      = read.csv(file = paste(path_data_folder,filename_feat,sep= ""))
dat_feat      = dat_feat[,-c(1,3:7)]

filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))
dat_bird$clade = c("Accipitriformes","Accipitriformes","Galloanserae","Galloanserae","Aequorlitornithes","Aequorlitornithes","Galloanserae","Strisores","Strisores",
                    "Galloanserae","Coraciimorphae","Coraciimorphae","Columbaves","Columbaves","Columbaves","Passeriformes","Passeriformes","Passeriformes","Passeriformes","Strisores",
                    "Australaves","Australaves","Australaves","Australaves","Galloanserae","Galloanserae","Coraciimorphae","Coraciimorphae",
                    "Coraciimorphae","Aequorlitornithes","Aequorlitornithes","Strigiformes","Strigiformes")
dat_bird$binomial = dat_bird$binomial.x

filename_body_e = list.files(path = path_data_folder, pattern = paste("bodyCGandMOI_extend"))
dat_body_e      = read.csv(file = paste(path_data_folder,filename_body_e,sep= ""))
dat_body_e      = reshape2::dcast(dat_body_e, species + BirdID + TestID + FrameID ~ component + object, value.var="value")

filename_body_t = list.files(path = path_data_folder, pattern = paste("bodyCGandMOI_tucked"))
dat_body_t      = read.csv(file = paste(path_data_folder,filename_body_t,sep= ""))
dat_body_t      = reshape2::dcast(dat_body_t, species + BirdID + TestID + FrameID ~ component + object, value.var="value")
remove(filename_wing,filename_feat,filename_bird,filename_body_e,filename_body_t)

filename_results = list.files(path = path_data_folder, pattern = paste("results"))
dat_results = read.csv(paste(path_data_folder,filename_results[1],sep=""))
for (i in 2:length(filename_results)){
  dat_results = rbind(dat_results,read.csv(paste(path_data_folder,filename_results[i],sep="")))
}
# adjust the sharp-shinned hawk 20_1016 to a coopers hawk
dat_wing$species[which(dat_wing$species == "acc_str" & dat_wing$BirdID == "20_1016")] = "acc_coo"
dat_body_e$species[which(dat_body_e$species == "acc_str" & dat_body_e$BirdID == "20_1016")] = "acc_coo"
dat_body_t$species[which(dat_body_t$species == "acc_str" & dat_body_t$BirdID == "20_1016")] = "acc_coo"
# remove wings where they accidentally got flipped around 44 from chr_amh and 7 from oce_leu CHECK THIS as you add more
# dat_results[which(dat_wing$pt2_X > 0 & dat_wing$species == "chr_amh")]

# Merge with all info we have about the individual
dat_all <- merge(dat_results,dat_wing[,c("species","BirdID","TestID","FrameID","elbow","manus","pt1_X","pt1_Y","pt1_Z",
                                         "pt6_X","pt6_Y","pt7_X","pt7_Y","pt8_X","pt8_Y","pt9_X","pt9_Y","pt10_X","pt10_Y",
                                         "pt11_X","pt11_Y","pt11_Z","pt12_X","pt12_Y","pt12_Z")], by = c("species","BirdID","TestID","FrameID"))
dat_all <- merge(dat_all,dat_bird[,-c(1,4)], by = c("species","BirdID"))
dat_all$torso_length <- dat_all$torsotail_length - dat_all$tail_length

dat_all_e <- merge(dat_all,dat_body_e[,-c(3,4)], by = c("species","BirdID"))
dat_all_t <- merge(dat_all,dat_body_t[,-c(3,4)], by = c("species","BirdID"))
dat_all_e$extend_neck = 1
dat_all_t$extend_neck = 0

### ----------- Combine all outputs ----------------

## ----------- Extended neck -----------
dat_all_e$full_VRP_Ixx <- dat_all_e$head_Ixx+dat_all_e$neck_Ixx+dat_all_e$tail_Ixx+dat_all_e$torso_Ixx+2*dat_all_e$wing_Ixx
dat_all_e$full_VRP_Iyy <- dat_all_e$head_Iyy+dat_all_e$neck_Iyy+dat_all_e$tail_Iyy+dat_all_e$torso_Iyy+2*dat_all_e$wing_Iyy
dat_all_e$full_VRP_Izz <- dat_all_e$head_Izz+dat_all_e$neck_Izz+dat_all_e$tail_Izz+dat_all_e$torso_Izz+2*dat_all_e$wing_Izz
dat_all_e$full_VRP_Ixz <- dat_all_e$head_Ixz+dat_all_e$neck_Ixz+dat_all_e$tail_Ixz+dat_all_e$torso_Ixz+2*dat_all_e$wing_Ixz

dat_all_e$full_CGx <- ((dat_all_e$head_CGx*dat_all_e$head_m)+(dat_all_e$neck_CGx*dat_all_e$neck_m)+(dat_all_e$tail_CGx*dat_all_e$tail_m)+
                     (dat_all_e$torso_CGx*dat_all_e$torso_m)+2*(dat_all_e$wing_CGx*dat_all_e$wing_m))/dat_all_e$full_m
dat_all_e$full_CGz <- ((dat_all_e$head_CGz*dat_all_e$head_m)+(dat_all_e$neck_CGz*dat_all_e$neck_m)+(dat_all_e$tail_CGz*dat_all_e$tail_m)+
                       (dat_all_e$torso_CGz*dat_all_e$torso_m)+2*(dat_all_e$wing_CGz*dat_all_e$wing_m))/dat_all_e$full_m

## ----------- Tucked neck -----------
dat_all_t$full_VRP_Ixx <- dat_all_t$head_Ixx+dat_all_t$tail_Ixx+dat_all_t$torso_Ixx+2*dat_all_t$wing_Ixx
dat_all_t$full_VRP_Iyy <- dat_all_t$head_Iyy+dat_all_t$tail_Iyy+dat_all_t$torso_Iyy+2*dat_all_t$wing_Iyy
dat_all_t$full_VRP_Izz <- dat_all_t$head_Izz+dat_all_t$tail_Izz+dat_all_t$torso_Izz+2*dat_all_t$wing_Izz
dat_all_t$full_VRP_Ixz <- dat_all_t$head_Ixz+dat_all_t$tail_Ixz+dat_all_t$torso_Ixz+2*dat_all_t$wing_Ixz

dat_all_t$full_CGx <- ((dat_all_t$head_CGx*dat_all_t$head_m)+(dat_all_t$tail_CGx*dat_all_t$tail_m)+
                         (dat_all_t$torso_CGx*dat_all_t$torso_m)+2*(dat_all_t$wing_CGx*dat_all_t$wing_m))/dat_all_t$full_m
dat_all_t$full_CGz <- ((dat_all_t$head_CGz*dat_all_t$head_m)+(dat_all_t$tail_CGz*dat_all_t$tail_m)+
                         (dat_all_t$torso_CGz*dat_all_t$torso_m)+2*(dat_all_t$wing_CGz*dat_all_t$wing_m))/dat_all_t$full_m

## -------------- Iterate through each wing shape to calculate
dat_all_e$S_proj <- 0
for(i in 1:nrow(dat_all_e)){
  # ---------- Extended neck results -------------
  I_vrp = matrix(0, nrow = 3, ncol = 3)
  I_vrp[1,1] = dat_all_e$full_VRP_Ixx[i]
  I_vrp[2,2] = dat_all_e$full_VRP_Iyy[i]
  I_vrp[3,3] = dat_all_e$full_VRP_Izz[i]
  I_vrp[3,1] = dat_all_e$full_VRP_Ixz[i]
  I_vrp[1,3] = dat_all_e$full_VRP_Ixz[i]

  CG = c(dat_all_e$full_CGx[i], 0, dat_all_e$full_CGz[i])
  I = parallelaxis(I_vrp,-CG,dat_all_e$full_m[i],"A")

  dat_all_e$full_Ixx[i] = I[1,1]
  dat_all_e$full_Iyy[i] = I[2,2]
  dat_all_e$full_Izz[i] = I[3,3]
  dat_all_e$full_Ixz[i] = I[3,1]
  # save the principal axes
  pri_axes = eigen(I)
  dat_all_e$intaxis_x = pri_axes$vectors[2,1]
  dat_all_e$intaxis_y = pri_axes$vectors[2,2]
  dat_all_e$intaxis_z = pri_axes$vectors[2,3]
  dat_all_e$majoraxis_x = pri_axes$vectors[1,1]
  dat_all_e$majoraxis_y = pri_axes$vectors[1,2]
  dat_all_e$majoraxis_z = pri_axes$vectors[1,3]
  dat_all_e$minoraxis_x = pri_axes$vectors[3,1]
  dat_all_e$minoraxis_y = pri_axes$vectors[3,2]
  dat_all_e$minoraxis_z = pri_axes$vectors[3,3]
  # ---------- Tucked neck results -------------
  I_vrp = matrix(0, nrow = 3, ncol = 3)
  I_vrp[1,1] = dat_all_t$full_VRP_Ixx[i]
  I_vrp[2,2] = dat_all_t$full_VRP_Iyy[i]
  I_vrp[3,3] = dat_all_t$full_VRP_Izz[i]
  I_vrp[3,1] = dat_all_t$full_VRP_Ixz[i]
  I_vrp[1,3] = dat_all_t$full_VRP_Ixz[i]

  CG = c(dat_all_t$full_CGx[i], 0, dat_all_t$full_CGz[i])
  I = parallelaxis(I_vrp,-CG,dat_all_t$full_m[i],"A")

  dat_all_t$full_Ixx[i] = I[1,1]
  dat_all_t$full_Iyy[i] = I[2,2]
  dat_all_t$full_Izz[i] = I[3,3]
  dat_all_t$full_Ixz[i] = I[3,1]
  # save the principal axes
  pri_axes = eigen(I)
  dat_all_t$intaxis_x = pri_axes$vectors[2,1]
  dat_all_t$intaxis_y = pri_axes$vectors[2,2]
  dat_all_t$intaxis_z = pri_axes$vectors[2,3]
  dat_all_t$majoraxis_x = pri_axes$vectors[1,1]
  dat_all_t$majoraxis_y = pri_axes$vectors[1,2]
  dat_all_t$majoraxis_z = pri_axes$vectors[1,3]
  dat_all_t$minoraxis_x = pri_axes$vectors[3,1]
  dat_all_t$minoraxis_y = pri_axes$vectors[3,2]
  dat_all_t$minoraxis_z = pri_axes$vectors[3,3]
  # Calculate the projected area for each wing - this is the correct order because X is negative and Y is positive
  x_vertices = c(dat_all_e$pt6_X[i],dat_all_e$pt7_X[i],dat_all_e$pt8_X[i],dat_all_e$pt9_X[i],dat_all_e$pt10_X[i],dat_all_e$pt11_X[i],dat_all_e$pt12_X[i])
  y_vertices = c(dat_all_e$pt6_Y[i],dat_all_e$pt7_Y[i],dat_all_e$pt8_Y[i],dat_all_e$pt9_Y[i],dat_all_e$pt10_Y[i],dat_all_e$pt11_Y[i],dat_all_e$pt12_Y[i])
  dat_all_e$S_proj[i] <- polyarea(x_vertices, y_vertices)
  dat_all_t$S_proj[i] <- dat_all_e$S_proj[i] # wing shapes have not changed between the neck positions
}

# Create the final combined data set
dat_final <- bind_rows(subset(dat_all_t, !(species %in% c("col_liv","lop_imp","lop_nyc","chr_amh","ana_pla","meg_alc","bra_can","aec_occ"))), # tucked neck species
                       subset(dat_all_e, species %in% c("col_liv","lop_imp","lop_nyc","chr_amh","ana_pla","meg_alc","bra_can","aec_occ"))) # extended neck species
dat_final$full_length = (dat_final$torsotail_length+dat_final$head_length+dat_final$neck_length)

dat_final <- dat_final[-c(which(dat_final$species == "cor_cor" & dat_final$BirdID == "20_2526")),]
dat_bird  <- dat_bird[-c(which(dat_bird$species == "cor_cor" & dat_bird$BirdID == "20_2526")),]

dat_final$shoulderx_specific  <- (dat_final$pt1_X-dat_final$head_length)/dat_final$full_length
dat_final$shouldery_specific  <- (dat_final$pt1_Y)/dat_final$full_length
dat_final$shoulderz_specific  <- (dat_final$pt1_Z)/dat_final$full_length

dat_final$full_Ixx_specific = dat_final$full_Ixx/(dat_final$full_m*dat_final$torso_length^2)
dat_final$full_Iyy_specific = dat_final$full_Iyy/(dat_final$full_m*dat_final$torso_length^2)
dat_final$full_Izz_specific = dat_final$full_Izz/(dat_final$full_m*dat_final$torso_length^2)
dat_final$full_Ixz_specific = dat_final$full_Ixz/(dat_final$full_m*dat_final$torso_length^2)

dat_final$full_CGx_specific <- (dat_final$full_CGx-dat_final$head_length)/dat_final$full_length
dat_final$full_CGz_specific <- (dat_final$full_CGz)/dat_final$body_height_max

dat_final$wing_span_cm[which(dat_final$species == "col_aur" & dat_final$BirdID == "21_0203")] = 0.489
dat_final$wing_CGy_specific = dat_final$wing_CGy/(0.5*dat_final$wing_span_cm)
dat_final$wing_CGx_specific = dat_final$wing_CGx/(0.5*dat_final$wing_span_cm)
dat_final$wing_CGz_specific = dat_final$wing_CGz/(0.5*dat_final$wing_span_cm)

dat_final$elbow_scaled = dat_final$elbow*0.001
dat_final$manus_scaled = dat_final$manus*0.001
### ------------- Compute summed quantities -----------------

# Maximum projected wing area
test       <- aggregate(list(S_proj_max = dat_final$S_proj),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird   <- merge(dat_bird,test, by = c("species","BirdID"))
dat_final  <- merge(dat_final,test, by = c("species","BirdID"))

# Shoulder position specific
test     <- aggregate(list(shoulderx_specific = dat_final$shoulderx_specific),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
test     <- aggregate(list(shoulderz_specific = dat_final$shoulderz_specific),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))

test     <- aggregate(list(full_m = dat_final$full_m),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))

### ---------------------- Compute inertial metrics -------------------------

# uses scaling from:
# Alerstam, T., Rosén, M., Bäckman, J., Ericson, P. G., & Hellgren, O. (2007).
# Flight speeds among bird species: allometric and phylogenetic effects. PLoS Biol, 5(8), e197.
dat_final$prop_q_dot     <- (-(dat_final$full_CGx-dat_final$pt1_X)*dat_final$S_proj_max*dat_final$full_m^0.24)/(dat_final$full_Iyy)
dat_final$prop_q_dot_nd  <- (-(dat_final$full_CGx-dat_final$pt1_X)*dat_final$S_proj_max*dat_final$torso_length^2)/(dat_final$full_Iyy)
dat_final$del_M_specific <- dat_final$prop_q_dot*dat_final$full_Iyy/(dat_final$full_m*dat_final$torso_length)

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_alldata.csv",sep="")
write.csv(dat_final,filename)
dat_final$elbow_scaled = dat_final$elbow*0.001
dat_final$manus_scaled = dat_final$manus*0.001

## ---- Compute the regression coefficients for each species for each variable -------
no_specimens <- nrow(dat_bird)
dat_model_out        <- data.frame(matrix(nrow = no_specimens*7*2, ncol = 5))
names(dat_model_out) <- c("species","model_variable","elb","man","elbman")
varlist_sp       <- c("full_Ixx_specific","full_Iyy_specific","full_Izz_specific","full_Ixz_specific",
                   "full_CGx_specific","wing_CGy_specific","full_CGz_specific")
short_varlist_sp <- c("Ixx_sp","Iyy_sp","Izz_sp","Ixz_sp","CGx_sp","CGy_sp","CGz_sp")
varlist_abs      <- c("full_Ixx","full_Iyy","full_Izz","full_Ixz",
                      "full_CGx","wing_CGy","full_CGz")
short_varlist_abs<- c("Ixx","Iyy","Izz","Ixz","CGx","CGy","CGz")
success = TRUE
count = 1
for (i in 1:no_specimens){
  # subset data to the current species
  tmp = subset(dat_final, species == dat_bird$species[i] & BirdID == dat_bird$BirdID[i])

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
    # null models
    model_null <- lapply(varlist, function(x) {lm(substitute(k ~ 1, list(k = as.name(x))), data = tmp)})
    # models with interactive effect removed
    model_red_elbman <- lapply(varlist, function(x) {lm(substitute(k ~ elbow_scaled+manus_scaled, list(k = as.name(x))), data = tmp)})
    # models with elbow only effect removed
    model_red_elb <- lapply(varlist, function(x) {lm(substitute(k ~ elbow_scaled:manus_scaled+manus_scaled, list(k = as.name(x))), data = tmp)})
    # models with wrist only effect removed
    model_red_man <- lapply(varlist, function(x) {lm(substitute(k ~ elbow_scaled:manus_scaled+elbow_scaled, list(k = as.name(x))), data = tmp)})

    CI_values      = lapply(models, confint)
    coef_values    = lapply(models, coef)

    for (j in 1:length(models)){
      dat_model_out$species[count]        <- dat_bird$species[i]
      dat_model_out$BirdID[count]         <- dat_bird$BirdID[i]
      dat_model_out$model_variable[count] <- short_varlist[j]
      dat_model_out$int[count]            <- coef_values[[j]]["(Intercept)"]
      dat_model_out$elb[count]            <- coef_values[[j]]["elbow_scaled"]
      dat_model_out$man[count]            <- coef_values[[j]]["manus_scaled"]
      dat_model_out$elbman[count]         <- coef_values[[j]]["elbow_scaled:manus_scaled"]
      dat_model_out$int_lb[count]         <- CI_values[[j]]["(Intercept)",1]
      dat_model_out$elb_lb[count]         <- CI_values[[j]]["elbow_scaled",1]
      dat_model_out$man_lb[count]         <- CI_values[[j]]["manus_scaled",1]
      dat_model_out$elbman_lb[count]      <- CI_values[[j]]["elbow_scaled:manus_scaled",1]
      dat_model_out$int_ub[count]         <- CI_values[[j]]["(Intercept)",2]
      dat_model_out$elb_ub[count]         <- CI_values[[j]]["elbow_scaled",2]
      dat_model_out$man_ub[count]         <- CI_values[[j]]["manus_scaled",2]
      dat_model_out$elbman_ub[count]      <- CI_values[[j]]["elbow_scaled:manus_scaled",2]

      dat_model_out$elb_eff[count]        <- cohen_f2(summary(models[[j]])$sigma,
                                                      summary(model_red_elb[[j]])$sigma,
                                                      summary(model_null[[j]])$sigma)
      dat_model_out$man_eff[count]        <- cohen_f2(summary(models[[j]])$sigma,
                                                      summary(model_red_man[[j]])$sigma,
                                                      summary(model_null[[j]])$sigma)
      dat_model_out$elbman_eff[count]     <- cohen_f2(summary(models[[j]])$sigma,
                                                      summary(model_red_elbman[[j]])$sigma,
                                                      summary(model_null[[j]])$sigma)
      count = count + 1
    }
  }
}

tmp           = reshape2::melt(dat_model_out, id = c("species","BirdID","model_variable"))
dat_model_out = reshape2::dcast(tmp, species + BirdID ~ model_variable + variable, value.var="value")
dat_model_out <- subset(dat_model_out, species != "oce_leu")

# Include basic geometry effects
tmp       = aggregate(list(mass = dat_bird$full_m),  by=list(species = dat_bird$species, binomial = dat_bird$binomial, BirdID = dat_bird$BirdID), mean)
dat_comp  = merge(tmp, dat_model_out, by = c("species","BirdID"))
dat_model_out  = merge(tmp, dat_model_out, by = c("species","BirdID"))
tmp       = aggregate(list(torsotail_length = dat_bird$torsotail_length),  by=list(species = dat_bird$species, BirdID = dat_bird$BirdID), mean)
dat_comp  = merge(dat_comp, tmp, by = c("species","BirdID"))

# Include other important factors
# Range of each component
test     <- aggregate(list(range_CGx = dat_final$full_CGx,
                           range_CGx_specific = dat_final$full_CGx_specific,
                           range_wing_CGy = dat_final$wing_CGy,
                           range_wing_CGy_specific = dat_final$wing_CGy_specific,
                           range_CGz = dat_final$full_CGz,
                           range_CGz_specific = dat_final$full_CGz_specific,
                           range_Ixx = dat_final$full_Ixx,
                           range_Ixx_specific = dat_final$full_Ixx_specific,
                           range_Iyy = dat_final$full_Iyy,
                           range_Iyy_specific = dat_final$full_Iyy_specific,
                           range_Izz = dat_final$full_Izz,
                           range_Izz_specific = dat_final$full_Izz_specific),
                      by=list(species = dat_final$species, BirdID = dat_final$BirdID), FUN=function(x){max(x)-min(x)})
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))

# Maximum values
test     <- aggregate(list(max_CGx = dat_final$full_CGx,
                           max_CGx_specific = dat_final$full_CGx_specific,
                           max_wing_CGy = dat_final$wing_CGy,
                           max_wing_CGy_specific = dat_final$wing_CGy_specific,
                           max_CGz = dat_final$full_CGz,
                           max_Ixx = dat_final$full_Ixx,
                           max_Ixx_specific = dat_final$full_Ixx_specific,
                           max_Iyy = dat_final$full_Ixx,
                           max_Iyy_specific = dat_final$full_Ixx_specific,
                           max_Izz = dat_final$full_Ixx,
                           max_Izz_specific = dat_final$full_Ixx_specific,
                           max_q = dat_final$prop_q_dot,
                           max_q_nd = dat_final$prop_q_dot_nd,
                           max_wingspan = dat_final$wing_span_cm,
                           max_length = dat_final$full_length),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))
# Minimum values
test     <- aggregate(list(min_CGx = dat_final$full_CGx,
                           min_CGx_specific = dat_final$full_CGx_specific,
                           min_CGz = dat_final$full_CGz,
                           min_wing_CGy_specific = dat_final$wing_CGy_specific,
                           min_CGz = dat_final$full_CGz,
                           min_CGz_specific = dat_final$full_CGz_specific),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_comp <- merge(dat_comp,test, by = c("species","BirdID"))

dat_comp <- subset(dat_comp, species != "oce_leu")

### --------------------- Phylogeny info ---------------------
## Read in tree
full_tree <-
  read.nexus("vikROM_passerines_403sp.tre")

## Species means
morpho_data_means <-
  mutate(dat_comp, phylo = binomial)# %>%
  ## turn binomial names into rownames
  #column_to_rownames(var = "binomial")

## Prune down the tree to the relevant species
sp_mean_matched <-
  geiger::treedata(
    phy = full_tree,
    data = morpho_data_means,
    sort = TRUE
  )

## pruned tree
## ladderization rotates nodes to make it easier to see basal vs derived
pruned_mcc <-
  ape::ladderize(
    sp_mean_matched$phy
  )
plot(pruned_mcc)

## matched species mean data
sp_mean_data <-
  subset(morpho_data_means,
         phylo %in% pruned_mcc$tip.label
  )

cohen_f2 <- function(V_ab, V_a, V_null){
  f2 = (((V_null-V_ab)/V_null)-((V_null-V_a)/V_null))/(1-((V_null-V_ab)/V_null))
  return(f2)
}
