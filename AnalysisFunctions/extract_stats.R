
### -------- PAPER STATS ----------

# allometric relationship for Iwing max
Ixx_val_mcmc_output

# ----- Comparison of pigeon data -------
# absolute values of Ixxwing
Ixx_max[which(Ixx_max$species == "col_liv"),c("species","wing_hum_Ixx")]
# maximum yCGwing difference between our values and Berg and Rayner
max(subset(dat_final, species == "col_liv")$wing_CGy_specific_orgShoulder)-(0.071/0.323)

# ---------- Full Bird CG -----------
# Full bird CGx range due to elbow and wrist
max(dat_comp$range_CGx)
max(dat_comp$range_CGx_specific)

dat_comp$species[which.max(dat_comp$range_CGx)]
dat_comp$species[which.max(dat_comp$range_CGx_specific)]
dat_comp$range_CGx[which.max(dat_comp$range_CGx_specific)]
# Full bird CGz range due to elbow and wrist
max(dat_comp$range_CGz)
max(dat_comp$range_CGz_specific)

dat_comp$species[which.max(dat_comp$range_CGz)]
dat_comp$species[which.max(dat_comp$range_CGz_specific)]
dat_comp$range_CGz[which.max(dat_comp$range_CGz_specific)]

#Significance of wrist on CGx - will be the same independent of the origin position
max(dat_comp$CGx_man_p)
min(dat_comp$CGx_man_etap)
#Significance of elbow on CGx - will be the same independent of the origin position
max(dat_comp$CGx_elb_p)
min(dat_comp$CGx_elb_etap)

#Check the sign - if different than can't say if it always moved forwards or aftwards
max(dat_comp$CGx_man)
min(dat_comp$CGx_man)
#Check the sign - if different than can't say if it always moved forwards or aftwards
max(dat_comp$CGx_elb)
min(dat_comp$CGx_elb)
# check to say if fair to say "tend to move forward"
length(which(dat_comp$CGx_elb > 0))

#Significance of elbow on CGz
max(dat_comp$CGz_elb_p)
min(dat_comp$CGz_elb_etap)
#Significance of wrist on CGz
max(dat_comp$CGz_man_p)
min(dat_comp$CGz_man_etap)

#Check the sign - if different than can't say if it always moved dorsally or ventrally
max(dat_comp$CGz_elb)
min(dat_comp$CGz_elb)
#Check the sign - if different than can't say if it always moved dorsally or ventrally
max(dat_comp$CGz_man)
min(dat_comp$CGz_man)
# check to say if fair to say "tend to move dorsally"
length(which(dat_comp$CGz_man < 0))
length(which(dat_comp$CGz_elb < 0))

# ---------------------- Shoulder CG effects ---------------

#Shoulder maximum CG range
max(shoulder_motion$range_CGx_specific)
shoulder_motion$species[which.max(shoulder_motion$range_CGx_specific)]
max(shoulder_motion$range_CGx)
shoulder_motion$species[which.max(shoulder_motion$range_CGx)]
#Shoulder minimum CG range
min(shoulder_motion$range_CGx_specific)
shoulder_motion$species[which.min(shoulder_motion$range_CGx_specific)]
min(shoulder_motion$range_CGx)
shoulder_motion$species[which.min(shoulder_motion$range_CGx)]

# obtain p-value for range plot with wing length
CG_range_model_mcmc_output

# ---------------------- Wing only CG effects ---------------

# Wing only CGy range due to elbow and wrist
max(dat_comp$range_wing_CGy_specific)
dat_comp$species[which.max(dat_comp$range_wing_CGy_specific)]
min(dat_comp$range_wing_CGy_specific)
dat_comp$species[which.min(dat_comp$range_wing_CGy_specific)]

max(dat_final$wing_CGy_specific_orgShoulder)
dat_final$species[which.max(dat_final$wing_CGy_specific_orgShoulder)]

max(dat_comp$CGy_elb_p)
min(dat_comp$CGy_elb_etap)

max(dat_comp$CGy_man_p)
min(dat_comp$CGy_man_etap)

## check how the wing position is scaling with wingspan - per Rayner
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
summary(pgls_model_mcmc) # if 95% overlap 1 this is indistinguishable from isometry

## --------------------- Moment of Inertia -------------------

# verify if overlapping with isometry 5/3, 1.666..
Ixx_model_mcmc_output
Iyy_model_mcmc_output
Izz_model_mcmc_output

tmp = aggregate(list(Ixx_range =  dat_comp$max_Ixx/dat_comp$min_Ixx),  by=list(species = dat_comp$species), max)
View(tmp)
tmp = aggregate(list(Izz_range =  dat_comp$max_Izz/dat_comp$min_Izz),  by=list(species = dat_comp$species), max)
View(tmp)
## Determine the significance values for MOI
max(dat_comp$Ixx_elb_p)
min(dat_comp$Ixx_elb_etap)
tmp = aggregate(list(Ixx_range =  dat_comp$Ixx_elb_p),  by=list(species = dat_comp$species), max)
View(tmp)
dat_comp$species[which.max(dat_comp$Ixx_elb_p)]

max(dat_comp$Izz_elb_p)
min(dat_comp$Izz_elb_etap)

##CHECK Ixz shifts
tmp = aggregate(list(Iyy_range =  dat_comp$max_Iyy/dat_comp$min_Iyy),  by=list(species = dat_comp$species), max)
View(tmp)
tmp = aggregate(list(Ixz_range =  abs(dat_comp$max_Ixz-dat_comp$min_Ixz)),  by=list(species = dat_comp$species), max)
View(tmp)
## ---------------------- Agility -----------------------
max_q_model_mcmc_output

dat_comp$species[which.max(dat_comp$max_q_nd)]
View(dat_comp[,c("species","max_q","max_q_nd")])

# check that this does not scale with mass
max_q_nd_model_mcmc_output
## ---------------------- Evolution -----------------------
OU_maxstab$opt$aicc-BM_maxstab$opt$aicc
OU_minstab$opt$aicc-BM_minstab$opt$aicc

OU_maxstab$opt
OU_minstab$opt

OU_xcg$opt$aicc-BM_xcg$opt$aicc
OU_xcg$opt
exp(OU_xcg$opt$z0)
