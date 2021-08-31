## ------ Unstable maximum agility --------
max_q_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    max_q_nd ~ full_m,
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = dat_comp,
    prior = univ_prior,
    nitt = 26000000, thin = 20000, burnin = 6000000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
max_q_nd_model_mcmc_output  = summary(max_q_nd_model_mcmc)

## ------ Stable maximum agility --------
min_q_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    min_q_nd ~ full_m,
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = dat_comp,
    prior = univ_prior,
    nitt = 26000000, thin = 20000, burnin = 6000000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
min_q_nd_model_mcmc_output  = summary(min_q_nd_model_mcmc)


# ---- Check if robust to the removal of the petrel ----
sp_mean_matched_adj <- keep.tip(phy = full_tree, tip = dat_comp$binomial[which(dat_comp$binomial != "oce_leu")])
## ladderization rotates nodes to make it easier to see basal vs derived
pruned_mcc_adj      <- ape::ladderize(sp_mean_matched)
inv.phylo_adj <- inverseA(pruned_mcc_adj, nodes = "TIPS", scale = TRUE)
max_q_nd_model_mcmc_adj <-
  MCMCglmm::MCMCglmm(
    max_q_nd ~ full_m,
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo_adj$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = dat_comp[which(dat_comp$binomial != "oce_leu"),],
    prior = univ_prior,
    nitt = 26000000, thin = 20000, burnin = 6000000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
max_q_nd_model_mcmc_output_adj  = summary(max_q_nd_model_mcmc_adj)

min_q_nd_model_mcmc_adj <-
  MCMCglmm::MCMCglmm(
    min_q_nd ~ full_m,
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo_adj$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = dat_comp[which(dat_comp$binomial != "oce_leu"),],
    prior = univ_prior,
    nitt = 26000000, thin = 20000, burnin = 6000000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
min_q_nd_model_mcmc_output_adj  = summary(min_q_nd_model_mcmc_adj)

range_q_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    q_range ~ full_m,
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = dat_comp,
    prior = univ_prior,
    nitt = 26000000, thin = 20000, burnin = 6000000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
range_q_nd_model_mcmc_output  = summary(range_q_nd_model_mcmc)

# create a data frame that pulls out the maximum and minimum agility per species not just per individual
tmp1               = aggregate(list(min_sm_nd = dat_comp$min_sm_nd,
                                    min_q_nd = dat_comp$min_q_nd), by = list(species = dat_comp$species), min)
tmp2               = aggregate(list(max_sm_nd = dat_comp$max_sm_nd,
                                    max_q_nd = dat_comp$max_q_nd), by = list(species = dat_comp$species), max)
stab_check               = merge(tmp1,tmp2, id = c("species"))
stab_check$species_order = factor(stab_check$species, levels = phylo_order$species)


range_sm_nd_model_mcmc <-
  MCMCglmm::MCMCglmm(
    sm_range ~ full_m,
    random = ~ phylo,
    scale = FALSE, ## whether you use this is up to you -- whatever is fair
    ginverse = list(phylo = inv.phylo$Ainv),
    family = c("gaussian"), ## errors are modeled as drawn from a Gaussian
    data = dat_comp,
    prior = univ_prior,
    nitt = 26000000, thin = 20000, burnin = 6000000,
    verbose = FALSE, ## switch this to TRUE if you feel like it
    pr = TRUE, pl = TRUE ## this saves some model output stuff
  )
range_sm_nd_model_mcmc_output  = summary(range_sm_nd_model_mcmc)
