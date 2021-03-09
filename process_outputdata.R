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
dat_wing     = dat_wing[,-c(1)]

filename_feat = list.files(path = path_data_folder, pattern = paste("allspecimen_featherinfo"))
dat_feat      = read.csv(file = paste(path_data_folder,filename_feat,sep= ""))
dat_feat      = dat_feat[,-c(1,3:7)]

filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))
dat_bird$clade = c("Accipitriformes","Accipitriformes","Aequorlitornithes","Aequorlitornithes","Galloanserae","Galloanserae","Aequorlitornithes","Aequorlitornithes","Galloanserae","Strisores","Strisores",
                    "Galloanserae","Coraciimorphae","Coraciimorphae","Columbaves","Columbaves","Columbaves","Passeriformes","Passeriformes","Passeriformes","Passeriformes","Strisores",
                    "Australaves","Australaves","Australaves","Australaves","Galloanserae","Galloanserae","Coraciimorphae","Coraciimorphae",
                    "Coraciimorphae","Aequorlitornithes","Aequorlitornithes","Strigiformes","Strigiformes")
names(dat_bird)[names(dat_bird) == "binomial.x"] = "binomial"
dat_bird     = dat_bird[,-c(1)]

filename_body = list.files(path = path_data_folder, pattern = paste("bodyCGandMOI_correct"))
dat_body      = read.csv(file = paste(path_data_folder,filename_body,sep= ""))
dat_body      = reshape2::dcast(dat_body, species + BirdID + TestID + FrameID ~ component + object, value.var="value")

# Read in each individual result
filename_results = list.files(path = path_data_folder, pattern = paste("results"))
dat_results = read.csv(paste(path_data_folder,filename_results[1],sep=""))
for (i in 2:length(filename_results)){
  dat_results = rbind(dat_results,read.csv(paste(path_data_folder,filename_results[i],sep="")))
}
# clean up environments
remove(filename_wing,filename_feat,filename_bird,filename_body,filename_results)

# -------------  Merge with all info -----------------
dat_final = merge(dat_results,dat_wing[,c("species","BirdID","TestID","FrameID","BirdID_FrameSpec","elbow","manus","pt1_X","pt1_Y","pt1_Z",
                                         "pt6_X","pt6_Y","pt7_X","pt7_Y","pt8_X","pt8_Y","pt8_Z","pt9_X","pt9_Y","pt10_X","pt10_Y",
                                         "pt11_X","pt11_Y","pt11_Z","pt12_X","pt12_Y","pt12_Z")], by = c("species","BirdID","TestID","FrameID"))
dat_final = merge(dat_final,dat_bird, by = c("species","BirdID"))

dat_final = merge(dat_final,dat_body[,-c(3,4)], by = c("species","BirdID"))

dat_final$full_length  = (dat_final$torsotail_length+dat_final$head_length+dat_final$neck_length)
dat_final$torso_length = dat_final$torsotail_length - dat_final$tail_length

dat_final$shoulderx_specific = (dat_final$pt1_X-dat_final$head_length)/dat_final$full_length
dat_final$shouldery_specific = (dat_final$pt1_Y)/dat_final$full_length
dat_final$shoulderz_specific = (dat_final$pt1_Z)/dat_final$full_length

dat_final$full_Ixx_specific = dat_final$full_Ixx/(dat_final$full_m*dat_final$wing_span_cm^2)
dat_final$full_Iyy_specific = dat_final$full_Iyy/(dat_final$full_m*dat_final$full_length^2)
dat_final$full_Izz_specific = dat_final$full_Izz/(dat_final$full_m*dat_final$wing_span_cm*dat_final$full_length)
dat_final$full_Ixz_specific = dat_final$full_Ixz/(dat_final$full_m*dat_final$wing_span_cm*dat_final$full_length)

dat_final$full_CGx_specific = (dat_final$full_CGx-dat_final$head_length)/dat_final$full_length
dat_final$full_CGz_specific = (dat_final$full_CGz)/dat_final$body_height_max

dat_final$wing_CGy_specific = dat_final$wing_CGy/(0.5*dat_final$wing_span_cm)
dat_final$wing_CGx_specific = dat_final$wing_CGx/(0.5*dat_final$wing_span_cm)
dat_final$wing_CGz_specific = dat_final$wing_CGz/(0.5*dat_final$wing_span_cm)

dat_final$elbow_scaled = dat_final$elbow*0.001
dat_final$manus_scaled = dat_final$manus*0.001

### ------------- Compute summed quantities -----------------

# Maximum projected wing area and maximum wing span
test       <- aggregate(list(S_proj_max = dat_final$wing_S_proj),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
dat_bird   <- merge(dat_bird,test, by = c("species","BirdID"))
dat_final  <- merge(dat_final,test, by = c("species","BirdID"))

# Shoulder position specific and mass
test     <- aggregate(list(shoulderx_specific = dat_final$shoulderx_specific,
                           shoulderz_specific = dat_final$shoulderz_specific,
                           full_m = dat_final$full_m),  by=list(species = dat_final$species, BirdID = dat_final$BirdID), max)
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

### --------------------- Phylogeny info ---------------------
## Read in tree
full_tree <-
  read.nexus("vikROM_passerines_403sp.tre")

## Species means
morpho_data_means <-
  mutate(dat_comp, phylo = binomial) #%>%
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
