
max_q_model_mcmc <-
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
max_q_model_mcmc_output  = summary(max_q_model_mcmc)

max_q_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    log(max_q_nd) ~ log(full_m),
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
max_q_nd_model_mcmc_output  = summary(max_q_nd_model_mcmc)


tmp1               = aggregate(list(c_l_max = 0.25*dat_final$chord/-(dat_final$full_CGx_orgShoulder)), by = list(species = dat_final$species, BirdID = dat_final$BirdID), max)
tmp2               = aggregate(list(c_l_min = 0.25*dat_final$chord/-(dat_final$full_CGx_orgShoulder)), by = list(species = dat_final$species, BirdID = dat_final$BirdID), min)

tmp3               = merge(tmp1,tmp2, id = c("species","BirdID"))
tmp3               = merge(tmp3,dat_comp[,c("species","BirdID",
                                            "max_q","max_q_nd",
                                            "min_q","min_q_nd")], id = c("species","BirdID"))

tmp1               = aggregate(list(c_l_min = tmp3$c_l_min,
                                    min_q_nd = tmp3$min_q_nd), by = list(species = tmp3$species), min)
tmp2               = aggregate(list(c_l_max = tmp3$c_l_max,
                                    max_q_nd = tmp3$max_q_nd), by = list(species = tmp3$species), max)
stab_check               = merge(tmp1,tmp2, id = c("species","BirdID"))
stab_check$species_order = factor(stab_check$species, levels = phylo_order$species)


max_sm_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    max_sm_nd ~ full_m,
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
max_sm_nd_model_mcmc_output  = summary(max_sm_nd_model_mcmc)

min_sm_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    min_sm_nd ~ full_m,
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
min_sm_nd_model_mcmc_output  = summary(min_sm_nd_model_mcmc)
