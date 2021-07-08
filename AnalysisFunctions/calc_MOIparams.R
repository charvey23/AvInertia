## -------------------- Calculate the Ixx_wing about the humeral head ------------------
## necessary for a direct comparison
library(AvInertia)

shift_Iorigin <- function(input_I,input_origin,input_CG,input_cg_or_a,input_m,new_origin,name,dat_info){
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  colnames(mass_properties) = c("species","BirdID","TestID","FrameID",
                                "component","object","value")
  dat    = list()
  dat$I  = matrix(0, nrow = 3, ncol = 3)
  dat$CG = matrix(0, nrow = 3, ncol = 1)

  #I defined about the wing CG
  if(input_cg_or_a == "A"){
    I_CG  = parallelaxis(input_I,(input_origin-input_CG),input_m,"A")
  }else{
    I_CG = input_I
  }

  dat$CG = new_origin - input_CG

  dat$I  = parallelaxis(I_CG,-dat$CG,input_m,"CG")

  dat$m = input_m
  new_row = store_data(dat_info,dat,mass_properties,name)
}


for (i in which(dat_final$wing_Ixx %in% dat_comp$max_wing_Ixx)){

  wing    = list()
  wing$I  = matrix(0, nrow = 3, ncol = 3)
  wing$CG = matrix(0, nrow = 3, ncol = 1)

  #current I is about the VRP
  wing$I[1,1] = dat_final$wing_Ixx[i]
  wing$I[2,2] = dat_final$wing_Iyy[i]
  wing$I[3,3] = dat_final$wing_Izz[i]

  wing$I[1,2] = dat_final$wing_Ixy[i]
  wing$I[2,1] = dat_final$wing_Ixy[i]
  wing$I[2,3] = dat_final$wing_Iyz[i]
  wing$I[3,2] = dat_final$wing_Iyz[i]
  wing$I[1,3] = dat_final$wing_Ixz[i]
  wing$I[3,1] = dat_final$wing_Ixz[i]

  wing$CG[1] = dat_final$wing_CGx[i]
  wing$CG[2] = dat_final$wing_CGy[i]
  wing$CG[3] = dat_final$wing_CGz[i]

  wing$m = dat_final$wing_m[i]

  hum_head = c(dat_final$pt1_X[i],
               dat_final$pt1_Y[i],
               dat_final$pt1_Z[i])

  new_row     = shift_Iorigin(wing$I,c(0,0,0),wing$CG,"A",wing$m,hum_head,"wing_hum",dat_final[i,])

  if(i == 804){
    Ixx_max = new_row
  } else{
    Ixx_max = rbind(Ixx_max,new_row)
  }
}
Ixx_max = reshape2::dcast(Ixx_max, species + BirdID + TestID + FrameID ~ component + object, value.var="value")
Ixx_max = merge(Ixx_max,dat_comp, id = c("species","BirdID"))

## --------- Fit the PGLS model to extract relationship between body mass and elements ------
Ixx_val_mcmc <-
  MCMCglmm::MCMCglmm(
    log(wing_hum_Ixx) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = Ixx_max,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
Ixx_val_mcmc_output  = summary(Ixx_val_mcmc)


## --------------------- Calculate the inertia PGLS relationship ----------------------
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


### ----------------- Compute the contributions of each component about the moment of inertia for the most extended configuration -----------------

dat_final$dist_joints = sqrt(dat_final$elbow^2+dat_final$manus^2)
tmp1 <- aggregate(list(dist_joints = dat_final$dist_joints), by = list(species = dat_final$species), max)
tmp1 <- merge(tmp1,dat_final, id = c("species","dist_joints"))
tmp1 <- tmp1[-which(duplicated(tmp1$species)),]
tmp1$type_metric = "max"
tmp2 <- aggregate(list(dist_joints = dat_final$dist_joints), by = list(species = dat_final$species), min)
tmp2 <- merge(tmp2,dat_final, id = c("species","dist_joints"))
tmp2$type_metric = "min"
tmp2 <- tmp2[-which(duplicated(tmp2$species)),]

I_contr <- rbind(tmp1,tmp2)

I_contr$species_order = factor(I_contr$species, levels = phylo_order$species)

##---------------------------------------------------------------------------------------------------------------------------------
##-------------------------------- Compute the contributions of each component to FULL CG origin I --------------------------------
##---------------------------------------------------------------------------------------------------------------------------------
# --------------------- Initialize variables -----------------------
mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
column_names = c("species","BirdID","TestID","FrameID","prop_type","component","value")
for (i in 1:nrow(I_contr)){
  #current I of each component is about the VRP will first calculate the contribution to the full body flight

  full_CG = c(I_contr$full_CGx[i], I_contr$full_CGy[i], I_contr$full_CGz[i])
  dat    = list()
  dat$I  = matrix(0, nrow = 3, ncol = 3)
  dat$CG = matrix(0, nrow = 3, ncol = 1)

  ##-------------------------------------------
  ##---------------- Wing ---------------------
  ##-------------------------------------------

  ## ----------- Bones ------------
  bones = dat

  #current I is about the VRP
  bones$I[1,1] = I_contr$bones_Ixx[i]
  bones$I[2,2] = I_contr$bones_Iyy[i]
  bones$I[3,3] = I_contr$bones_Izz[i]

  bones$I[1,2] = I_contr$bones_Ixy[i]
  bones$I[2,1] = I_contr$bones_Ixy[i]
  bones$I[2,3] = I_contr$bones_Iyz[i]
  bones$I[3,2] = I_contr$bones_Iyz[i]
  bones$I[1,3] = I_contr$bones_Ixz[i]
  bones$I[3,1] = I_contr$bones_Ixz[i]

  bones$CG[1] = I_contr$bones_CGx[i]
  bones$CG[2] = I_contr$bones_CGy[i]
  bones$CG[3] = I_contr$bones_CGz[i]

  bones$m     = I_contr$bones_m[i]

  bones_row     = shift_Iorigin(bones$I,c(0,0,0),bones$CG,"A",bones$m,full_CG,"bone",I_contr[i,])

  ## ----------- Muscles ------------
  muscles    = dat

  #current I is about the VRP
  muscles$I[1,1] = I_contr$muscles_Ixx[i]
  muscles$I[2,2] = I_contr$muscles_Iyy[i]
  muscles$I[3,3] = I_contr$muscles_Izz[i]

  muscles$I[1,2] = I_contr$muscles_Ixy[i]
  muscles$I[2,1] = I_contr$muscles_Ixy[i]
  muscles$I[2,3] = I_contr$muscles_Iyz[i]
  muscles$I[3,2] = I_contr$muscles_Iyz[i]
  muscles$I[1,3] = I_contr$muscles_Ixz[i]
  muscles$I[3,1] = I_contr$muscles_Ixz[i]

  muscles$CG[1] = I_contr$muscles_CGx[i]
  muscles$CG[2] = I_contr$muscles_CGy[i]
  muscles$CG[3] = I_contr$muscles_CGz[i]

  muscles$m     = I_contr$muscles_m[i]

  muscles_row     = shift_Iorigin(muscles$I,c(0,0,0),muscles$CG,"A",muscles$m,full_CG,"muscle",I_contr[i,])


  ## ----------- Feathers ------------
  feathers    = dat

  #current I is about the VRP
  feathers$I[1,1] = I_contr$feathers_Ixx[i]
  feathers$I[2,2] = I_contr$feathers_Iyy[i]
  feathers$I[3,3] = I_contr$feathers_Izz[i]

  feathers$I[1,2] = I_contr$feathers_Ixy[i]
  feathers$I[2,1] = I_contr$feathers_Ixy[i]
  feathers$I[2,3] = I_contr$feathers_Iyz[i]
  feathers$I[3,2] = I_contr$feathers_Iyz[i]
  feathers$I[1,3] = I_contr$feathers_Ixz[i]
  feathers$I[3,1] = I_contr$feathers_Ixz[i]

  feathers$CG[1] = I_contr$feathers_CGx[i]
  feathers$CG[2] = I_contr$feathers_CGy[i]
  feathers$CG[3] = I_contr$feathers_CGz[i]

  feathers$m     = I_contr$feathers_m[i]

  feathers_row     = shift_Iorigin(feathers$I,c(0,0,0),feathers$CG,"A",feathers$m,full_CG,"feathers",I_contr[i,])

  ## ----------- Skin ------------
  skin    = dat

  #current I is about the VRP
  skin$I[1,1] = I_contr$skin_Ixx[i]
  skin$I[2,2] = I_contr$skin_Iyy[i]
  skin$I[3,3] = I_contr$skin_Izz[i]

  skin$I[1,2] = I_contr$skin_Ixy[i]
  skin$I[2,1] = I_contr$skin_Ixy[i]
  skin$I[2,3] = I_contr$skin_Iyz[i]
  skin$I[3,2] = I_contr$skin_Iyz[i]
  skin$I[1,3] = I_contr$skin_Ixz[i]
  skin$I[3,1] = I_contr$skin_Ixz[i]

  skin$CG[1] = I_contr$skin_CGx[i]
  skin$CG[2] = I_contr$skin_CGy[i]
  skin$CG[3] = I_contr$skin_CGz[i]

  skin$m     = I_contr$skin_m[i]

  skin_row     = shift_Iorigin(skin$I,c(0,0,0),skin$CG,"A",skin$m,full_CG,"skin",I_contr[i,])

  ##-------------------------------------------
  ##---------------- Body ---------------------
  ##-------------------------------------------

  ## ----------- Head ------------
  head    = dat

  #current I is about the VRP
  head$I[1,1] = I_contr$head_Ixx[i]
  head$I[2,2] = I_contr$head_Iyy[i]
  head$I[3,3] = I_contr$head_Izz[i]

  head$I[1,2] = I_contr$head_Ixy[i]
  head$I[2,1] = I_contr$head_Ixy[i]
  head$I[2,3] = I_contr$head_Iyz[i]
  head$I[3,2] = I_contr$head_Iyz[i]
  head$I[1,3] = I_contr$head_Ixz[i]
  head$I[3,1] = I_contr$head_Ixz[i]

  head$CG[1] = I_contr$head_CGx[i]
  head$CG[2] = I_contr$head_CGy[i]
  head$CG[3] = I_contr$head_CGz[i]

  head$m     = I_contr$head_m[i]

  head_row     = shift_Iorigin(head$I,c(0,0,0),head$CG,"A",head$m,full_CG,"head",I_contr[i,])

  ## ----------- Torso ------------
  torso    = dat

  #current I is about the VRP
  torso$I[1,1] = I_contr$torso_Ixx[i]
  torso$I[2,2] = I_contr$torso_Iyy[i]
  torso$I[3,3] = I_contr$torso_Izz[i]

  torso$I[1,2] = I_contr$torso_Ixy[i]
  torso$I[2,1] = I_contr$torso_Ixy[i]
  torso$I[2,3] = I_contr$torso_Iyz[i]
  torso$I[3,2] = I_contr$torso_Iyz[i]
  torso$I[1,3] = I_contr$torso_Ixz[i]
  torso$I[3,1] = I_contr$torso_Ixz[i]

  torso$CG[1] = I_contr$torso_CGx[i]
  torso$CG[2] = I_contr$torso_CGy[i]
  torso$CG[3] = I_contr$torso_CGz[i]

  torso$m     = I_contr$torso_m[i]

  torso_row     = shift_Iorigin(torso$I,c(0,0,0),torso$CG,"A",torso$m,full_CG,"torso",I_contr[i,])

  ## ----------- Tail ------------
  tail    = dat

  #current I is about the VRP
  tail$I[1,1] = I_contr$tail_Ixx[i]
  tail$I[2,2] = I_contr$tail_Iyy[i]
  tail$I[3,3] = I_contr$tail_Izz[i]

  tail$I[1,2] = I_contr$tail_Ixy[i]
  tail$I[2,1] = I_contr$tail_Ixy[i]
  tail$I[2,3] = I_contr$tail_Iyz[i]
  tail$I[3,2] = I_contr$tail_Iyz[i]
  tail$I[1,3] = I_contr$tail_Ixz[i]
  tail$I[3,1] = I_contr$tail_Ixz[i]

  tail$CG[1] = I_contr$tail_CGx[i]
  tail$CG[2] = I_contr$tail_CGy[i]
  tail$CG[3] = I_contr$tail_CGz[i]

  tail$m     = I_contr$tail_m[i]

  tail_row     = shift_Iorigin(tail$I,c(0,0,0),tail$CG,"A",tail$m,full_CG,"tail",I_contr[i,])

  ## ----------- Neck ------------
  neck    = dat

  #current I is about the VRP
  neck$I[1,1] = I_contr$neck_Ixx[i]
  neck$I[2,2] = I_contr$neck_Iyy[i]
  neck$I[3,3] = I_contr$neck_Izz[i]

  neck$I[1,2] = I_contr$neck_Ixy[i]
  neck$I[2,1] = I_contr$neck_Ixy[i]
  neck$I[2,3] = I_contr$neck_Iyz[i]
  neck$I[3,2] = I_contr$neck_Iyz[i]
  neck$I[1,3] = I_contr$neck_Ixz[i]
  neck$I[3,1] = I_contr$neck_Ixz[i]

  neck$CG[1] = I_contr$neck_CGx[i]
  neck$CG[2] = I_contr$neck_CGy[i]
  neck$CG[3] = I_contr$neck_CGz[i]

  neck$m     = I_contr$neck_m[i]

  neck_row     = shift_Iorigin(neck$I,c(0,0,0),neck$CG,"A",neck$m,full_CG,"neck",I_contr[i,])

  ### ---- Just to check that the shift worked correctly ----
  I_contr_fullCG = rbind(feathers_row,skin_row,muscles_row,bones_row,head_row,neck_row,torso_row,tail_row)
  # replace any non-existent neck as 0
  I_contr_fullCG$value[which(is.na(I_contr_fullCG$value)& I_contr_fullCG$component == "neck")] <- 0
  # there are two wings which is why add two

  Ixx_CG = sum(subset(I_contr_fullCG, object == "Ixx")$value,
      subset(I_contr_fullCG, component == "feathers" & object == "Ixx")$value,
      subset(I_contr_fullCG, component == "skin" & object == "Ixx")$value,
      subset(I_contr_fullCG, component == "bone" & object == "Ixx")$value,
      subset(I_contr_fullCG, component == "muscle" & object == "Ixx")$value)
  if(abs(Ixx_CG - I_contr$full_Ixx[i]) > 10^-12){
    warning("Error: not correct match")
  }

  ### ------ Calculate the contributions ------
  # Bones
  I_contr$bone_con_Ixx[i] = 2*subset(I_contr_fullCG, component == "bone" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$bone_con_Iyy[i] = 2*subset(I_contr_fullCG, component == "bone" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$bone_con_Izz[i] = 2*subset(I_contr_fullCG, component == "bone" & object == "Izz")$value/I_contr$full_Izz[i]
  # Feathers
  I_contr$feat_con_Ixx[i] = 2*subset(I_contr_fullCG, component == "feathers" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$feat_con_Iyy[i] = 2*subset(I_contr_fullCG, component == "feathers" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$feat_con_Izz[i] = 2*subset(I_contr_fullCG, component == "feathers" & object == "Izz")$value/I_contr$full_Izz[i]
  # Muscles
  I_contr$musc_con_Ixx[i] = 2*subset(I_contr_fullCG, component == "muscle" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$musc_con_Iyy[i] = 2*subset(I_contr_fullCG, component == "muscle" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$musc_con_Izz[i] = 2*subset(I_contr_fullCG, component == "muscle" & object == "Izz")$value/I_contr$full_Izz[i]
  # Skin
  I_contr$skin_con_Ixx[i] = 2*subset(I_contr_fullCG, component == "skin" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$skin_con_Iyy[i] = 2*subset(I_contr_fullCG, component == "skin" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$skin_con_Izz[i] = 2*subset(I_contr_fullCG, component == "skin" & object == "Izz")$value/I_contr$full_Izz[i]
  # head
  I_contr$head_con_Ixx[i] = subset(I_contr_fullCG, component == "head" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$head_con_Iyy[i] = subset(I_contr_fullCG, component == "head" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$head_con_Izz[i] = subset(I_contr_fullCG, component == "head" & object == "Izz")$value/I_contr$full_Izz[i]
  # torso
  I_contr$torso_con_Ixx[i] = subset(I_contr_fullCG, component == "torso" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$torso_con_Iyy[i] = subset(I_contr_fullCG, component == "torso" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$torso_con_Izz[i] = subset(I_contr_fullCG, component == "torso" & object == "Izz")$value/I_contr$full_Izz[i]
  # tail
  I_contr$tail_con_Ixx[i] = subset(I_contr_fullCG, component == "tail" & object == "Ixx")$value/I_contr$full_Ixx[i]
  I_contr$tail_con_Iyy[i] = subset(I_contr_fullCG, component == "tail" & object == "Iyy")$value/I_contr$full_Iyy[i]
  I_contr$tail_con_Izz[i] = subset(I_contr_fullCG, component == "tail" & object == "Izz")$value/I_contr$full_Izz[i]
  # neck
  if(is.na(neck$I[1,1])){
    I_contr$neck_con_Ixx[i] = 0
    I_contr$neck_con_Iyy[i] = 0
    I_contr$neck_con_Izz[i] = 0
  }else{
    I_contr$neck_con_Ixx[i] = subset(I_contr_fullCG, component == "neck" & object == "Ixx")$value/I_contr$full_Ixx[i]
    I_contr$neck_con_Iyy[i] = subset(I_contr_fullCG, component == "neck" & object == "Iyy")$value/I_contr$full_Iyy[i]
    I_contr$neck_con_Izz[i] = subset(I_contr_fullCG, component == "neck" & object == "Izz")$value/I_contr$full_Izz[i]
  }
}
remove(feathers,skin,muscles,bones,head,neck,torso,tail,feathers_row,skin_row,muscles_row,bones_row,head_row,neck_row,torso_row,tail_row)



