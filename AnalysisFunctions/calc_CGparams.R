### --- Script containing the data manipulations needed to create all plots -----


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

CGx_sp_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(-mean_CGx_specific_orgShoulder) ~ log(full_m),
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGx_sp_model_mcmc_output  = summary(CGx_sp_model_mcmc)

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

## ---------- Calculate the shoulder motion range relation to span ratio ---------
CG_range_model_mcmc <-
  MCMCglmm::MCMCglmm(
    range_CGx_specific ~ span_ratio,
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = shoulder_motion,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CG_range_model_mcmc_output  = summary(CG_range_model_mcmc)

## ---------- Calculate the elbow and wrist CG range relation to arm to hand ratio ---------

CGy_range_armhand <-
  MCMCglmm::MCMCglmm(
    range_wing_CGy_specific ~ max_armhand_ratio,
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGy_range_armhand_output  = summary(CGy_range_armhand)

CGx_range_armhand <-
  MCMCglmm::MCMCglmm(
    range_CGx_specific ~ max_armhand_ratio,
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGx_range_armhand_output  = summary(CGx_range_armhand)

CGz_range_armhand <-
  MCMCglmm::MCMCglmm(
    range_CGz_specific ~ max_armhand_ratio,
    random = ~ phylo,
    scale = FALSE,
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"),
    data = dat_comp,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE,pr = TRUE, pl = TRUE)
CGz_range_armhand_output  = summary(CGz_range_armhand)
