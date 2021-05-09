
### -------- PAPER STATS ----------

# allometric relationship for Iwing max
pgls_model_Iwingmax <-
  MCMCglmm::MCMCglmm(
    log(max_wing_Ixx) ~ log(full_m),
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
summary(pgls_model_Iwingmax)

# ----- Max range of CG positions -------
# Cywing comparative discussion
# maximum Cywing relative to the shoulder
max(dat_final$wing_CGy_specific_orgShoulder)
dat_final$species[which.max(dat_final$wing_CGy_specific_orgShoulder)]

#maximum CGy for a pigeon
max(subset(dat_final, species == "col_liv")$wing_CGy_specific_orgShoulder)
#absolute differece
(max(subset(dat_final, species == "col_liv")$wing_CGy_specific_orgShoulder)-0.22)*max(subset(dat_final, species == "col_liv")$span)*0.5

# to highlight the link between relative humerus length and CGy range
pgls_model_hum_range <-
  MCMCglmm::MCMCglmm(
    log(range_wing_CGy_specific) ~ log(hum_maxspan),
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
summary(pgls_model_hum_range)

# Effect size of CGy
max(dat_comp$CGy_elb_p)
min(dat_comp$CGy_elb_etap)

# ---------- Full Bird CG -----------
# Full bird CGx range due to elbow and wrist
max(dat_comp$range_CGx)
max(dat_comp$range_CGx_specific)

dat_comp$species[which.max(dat_comp$range_CGx)]
dat_comp$species[which.max(dat_comp$range_CGx_specific)]

# Full bird CGz range due to elbow and wrist
max(dat_comp$range_CGz)
max(dat_comp$range_CGz_specific)

dat_comp$species[which.max(dat_comp$range_CGz)]
dat_comp$species[which.max(dat_comp$range_CGz_specific)]
dat_comp$range_CGz[which.max(dat_comp$range_CGz_specific)]

#Shoulder CGx range
max(shoulder_motion$range_CGx_specific)
shoulder_motion$species[which.max(shoulder_motion$range_CGx_specific)]
max(shoulder_motion$range_CGx)
shoulder_motion$species[which.max(shoulder_motion$range_CGx)]

#Shoulder CGz range
max(shoulder_motion$range_CGz_specific)
shoulder_motion$species[which.max(shoulder_motion$range_CGz_specific)]
max(shoulder_motion$range_CGz)
shoulder_motion$species[which.max(shoulder_motion$range_CGz)]

# Wing only CGy range due to elbow and wrist
max(dat_comp$range_wing_CGy)
dat_comp$species[which.max(dat_comp$range_wing_CGy)]
max(dat_comp$range_wing_CGy_specific)
dat_comp$species[which.max(dat_comp$range_wing_CGy_specific)]
min(dat_comp$range_wing_CGy_specific)
dat_comp$species[which.min(dat_comp$range_wing_CGy_specific)]


## check that the wing position isn't scaling with wingspan
pgls_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(max_wing_CGy) ~ log(max_wingspan),
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

## Determine the significance values for MOI
max(dat_comp$Ixx_elb_p)
min(dat_comp$Ixx_elb_etap)
dat_comp$species[which.max(dat_comp$Ixx_elb_p)]

## Longitudinal stablity
pgls_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(max_q) ~ log(full_m),
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
summary(pgls_model_mcmc)
## Create table of varied contributions for major body components to MOI
# Range of each component

pgls_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(max_q) ~ log(c_l_true),
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = test2,
    prior = univ_prior,
    nitt = 130000, thin = 100, burnin = 30000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
summary(pgls_model_mcmc)
