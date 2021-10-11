### -----------------------------------------
### -------- MAIN MANUSCRIPT STATS ----------
### -----------------------------------------

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

# NOTE: all of these elbow and wrist effect sizes will be the same independent
# of the origin position or normalization as these are constants for each species

#Significance of wrist on CGx
max(dat_comp$CGx_man_p)
min(dat_comp$CGx_man_etap)
#Significance of elbow on CGx
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

# If significantly different than isometry will scale with mass
CGx_sp_model_mcmc_output
CGz_sp_model_mcmc_output

# ---------------------- Shoulder CG effects ---------------

#Shoulder maximum CG range
max(shoulder_motion$range_CGx_specific)
shoulder_motion$species[which.max(shoulder_motion$range_CGx_specific)]
shoulder_motion$range_CGx[which.max(shoulder_motion$range_CGx_specific)]
#Shoulder minimum CG range
min(shoulder_motion$range_CGx_specific)
shoulder_motion$species[which.min(shoulder_motion$range_CGx_specific)]
shoulder_motion$range_CGx[which.min(shoulder_motion$range_CGx_specific)]

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

# arm to hand wing ratio
CGy_range_armhand_output
# ranges provided in the reviewer response
max(dat_comp$max_armhand_ratio)
min(dat_comp$max_armhand_ratio)
hist(dat_comp$max_armhand_ratio, breaks = seq(0,1.5,0.1))

max(dat_comp$CGy_elb_p)
min(dat_comp$CGy_elb_etap)

max(dat_comp$CGy_man_p)
min(dat_comp$CGy_man_etap)

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
tmp = aggregate(list(Ixx_p =  dat_comp$Ixx_elb_p),  by=list(species = dat_comp$species), max)
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
range_q_model_mcmc_output

# check this is true if we remove the storm petrels
range_q_model_mcmc_output_adj
# lowest negative is the most agile stable flyer
dat_comp$species[which.min(dat_comp$min_q_nd)]
# check out the most agile in the stable range
View(dat_comp[,c("species","min_q","min_q_nd")])

range_sm_nd_model_mcmc_output

#mean(dat_final$S_max[which(dat_final$species == "cyp_nig")]/dat_final$full_Iyy[which(dat_final$species == "cyp_nig")]) # this value is lower
#mean(dat_final$S_max[which(dat_final$species == "oce_leu")]/dat_final$full_Iyy[which(dat_final$species == "oce_leu")])

dat_comp$species[which(dat_comp$min_sm>0 & dat_comp$max_sm>0)] # number of individuals that remain stable
dat_comp$species[which(dat_comp$min_sm<0 & dat_comp$max_sm<0)] # number of individuals that remain unstable
## ---------------------- Evolution -----------------------
OU_maxsm$opt$aicc-BM_maxsm$opt$aicc
OU_minsm$opt$aicc-BM_minsm$opt$aicc
(length(which(out_maxstab$null > out_maxstab$lr))/5000)# p_value
(length(which(out_minstab$null > out_minstab$lr))/5000)# p_value
# max static margin
OU_maxsm$opt$z0
OU_maxsm$opt$alpha
OU_maxsm$opt$sigsq
# min static margin
OU_minsm$opt$z0
OU_minsm$opt$alpha
OU_minsm$opt$sigsq
# xcg
OU_xcg$opt$aicc-BM_xcg$opt$aicc
OU_xcg$opt$z0
OU_xcg$opt$alpha
OU_xcg$opt$sigsq

(length(which(out_xcg$null > out_xcg$lr))/5000) # p-value


### -----------------------------------------
### ------------ METHODS STATS --------------
### -----------------------------------------

# ----- Validation of AvInertia -------
Ixx_val_mcmc_output # allometric relationship for Iwing max
# compare upper 95% CI to lower 95% CI from Kirkpatrick  and Berg and Rayner
Ixx_val_mcmc_output$solutions[2,3] > 1.89  # Kirkpatick
Ixx_val_mcmc_output$solutions[2,3] > 1.787 # Berg and Rayner

# ----- Comparison of pigeon data -------
# absolute values of Ixxwing
Ixx_max[which(Ixx_max$species == "col_liv"),c("species","wing_hum_Ixx")]
# maximum yCGwing difference between our values and Berg and Rayner
max(subset(dat_final, species == "col_liv")$wing_CGy_specific_orgShoulder)-(0.071/0.323)

# ----- Extended Data Table 1 -----
CGx_sp_model_mcmc_output
CGz_sp_model_mcmc_output
CG_range_model_mcmc_output

CGy_sp_model_mcmc_output
CGy_range_armhand_output
CGx_range_armhand_output
CGz_range_armhand_output

Ixx_model_mcmc_output
Iyy_model_mcmc_output
Izz_model_mcmc_output
Ixx_val_mcmc_output

range_q_model_mcmc_output
range_sm_nd_model_mcmc_output

#check that the agility range still decreases with mass even once the petrels are removed
range_q_model_mcmc_output_adj

# -------- Descriptive stats
(length(which(out_xcg$null > out_xcg$lr))/5000)
(length(which(out_maxstab$null > out_maxstab$lr))/5000)
(length(which(out_minstab$null > out_minstab$lr))/5000)

# calculate power of tests
length(which(out_xcg$test > quantile(out_xcg$null, probs = 0.95)))/5000
length(which(out_maxstab$test > quantile(out_maxstab$null, probs = 0.95)))/5000
length(which(out_minstab$test > quantile(out_minstab$null, probs = 0.95)))/5000

# extract the 95% confidence intervals

# alpha
dist_alpha_xcg = out_xcg$par_dists$value[which(out_xcg$par_dists$comparison == "BB" & out_xcg$par_dists$parameter =="alpha")]
quantile(dist_alpha_xcg, probs = c(0.05,0.95))
dist_alpha_maxstab = out_maxstab$par_dists$value[which(out_maxstab$par_dists$comparison == "BB" & out_maxstab$par_dists$parameter =="alpha")]
quantile(dist_alpha_maxstab, probs = c(0.05,0.95))
dist_alpha_minstab = out_minstab$par_dists$value[which(out_minstab$par_dists$comparison == "BB" & out_minstab$par_dists$parameter =="alpha")]
quantile(dist_alpha_minstab, probs = c(0.05,0.95))

# z0
dist_z0_xcg = out_xcg$par_dists$value[which(out_xcg$par_dists$comparison == "BB" & out_xcg$par_dists$parameter =="z0")]
quantile(dist_z0_xcg, probs = c(0.05,0.95))
dist_z0_maxstab = out_maxstab$par_dists$value[which(out_maxstab$par_dists$comparison == "BB" & out_maxstab$par_dists$parameter =="z0")]
quantile(dist_z0_maxstab, probs = c(0.05,0.95))
dist_z0_minstab = out_minstab$par_dists$value[which(out_minstab$par_dists$comparison == "BB" & out_minstab$par_dists$parameter =="z0")]
quantile(dist_z0_minstab, probs = c(0.05,0.95))

# sigsq - CAUTION- THIS IS DIVIDED BY 1E-3 TO SIMPLIFY TABLE
dist_sigsq_xcg = out_xcg$par_dists$value[which(out_xcg$par_dists$comparison == "BB" & out_xcg$par_dists$parameter =="sigsq")]
quantile(dist_sigsq_xcg, probs = c(0.05,0.95))/1e-3
dist_sigsq_maxstab = out_maxstab$par_dists$value[which(out_maxstab$par_dists$comparison == "BB" & out_maxstab$par_dists$parameter =="sigsq")]
quantile(dist_sigsq_maxstab, probs = c(0.05,0.95))/1e-3
dist_sigsq_minstab = out_minstab$par_dists$value[which(out_minstab$par_dists$comparison == "BB" & out_minstab$par_dists$parameter =="sigsq")]
quantile(dist_sigsq_minstab, probs = c(0.05,0.95))/1e-3



