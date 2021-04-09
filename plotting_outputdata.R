library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe
library(gridExtra) # for using grid arrange
library(ggtree)    # for plotting phylogeny

th <- ggplot2::theme_classic() +
  ggplot2::theme(
    # Text
    axis.title  = ggplot2::element_text(size = 10),
    axis.text   = ggplot2::element_text(size = 10, colour = "black"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, unit = "pt")),
    axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
    # Axis line
    axis.line   = ggplot2::element_blank(),
    axis.ticks.length = ggplot2::unit(-5,"pt"),
    # Background transparency
    # Background of panel
    panel.background = ggplot2::element_rect(fill = "transparent"),
    # Background behind actual data points
    plot.background  = ggplot2::element_rect(fill = "transparent", color = NA))

## ----------------- Plot of linear model outputs ------------------
# order species to match the phylogeny

man_fixed = c(80,100,120)*0.001
elb_fixed = c(80,100,120)*0.001
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
col_xz = "#af0069"

col_elb = "#55b3b1"
col_man = "#af0069"
col_elbman = "#5F0995"

phylo_label = c("Cooper's hawk","Sharp-shinned hawk", "Western grebe", "Mallard","Great blue heron","Canada goose","Common nighthawk",
                "Lady Amherst's pheasant","Northern flicker","Pigeon","Common raven","Steller's jay","Black swift",
                "Merlin","Peregrine falcon","Leach's storm petrel","Glaucous-winged gull","Himalayan monal","Silver pheasant",
                "Belted kingfisher","American white pelican","Barn owl")
dat_phylo_plot <- data.frame(label = pruned_mcc$tip.label, genus = phylo_label)
phylo_plot <- ggtree(pruned_mcc) %<+% dat_phylo_plot +
  geom_tiplab(aes(label=genus)) + xlim(0,150) + theme(plot.margin = unit(c(14,8,14,8), "mm"))

phylo_data <- ggtree::ggtree(pruned_mcc)[["data"]][,c("label","y")]
phylo_order <- merge(phylo_data[complete.cases(phylo_data),], unique(dat_bird[,c("binomial","species")]), by.x = "label", by.y = "binomial")
phylo_order <- arrange(phylo_order,y)

dat_comp$species_order = factor(dat_comp$species, levels = phylo_order$species)
dat_final$species_order = factor(dat_final$species, levels = phylo_order$species)

dat_model_out$species <- factor(dat_model_out$species, levels = phylo_order$species)
shading <- data.frame(col1 = levels(dat_model_out$species)[seq(from = 1, to = max(as.numeric(as.factor(dat_model_out$species)))-1, by = 2)],
                      col2 = levels(dat_model_out$species)[seq(from = 2, to = max(as.numeric(as.factor(dat_model_out$species))), by = 2)])


min_lb     <- aggregate(list(min_CGx_sp_int_lb = dat_model_out$CGx_sp_int_lb,
                           min_CGy_sp_int_lb = dat_model_out$CGy_sp_int_lb,
                           min_CGz_sp_int_lb = dat_model_out$CGz_sp_int_lb,
                           min_CGx_sp_elb_lb1 = dat_model_out$CGx_sp_elb_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[1],
                           min_CGx_sp_elb_lb2 = dat_model_out$CGx_sp_elb_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[3],
                           min_CGy_sp_elb_lb1 = dat_model_out$CGy_sp_elb_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[1],
                           min_CGy_sp_elb_lb2 = dat_model_out$CGy_sp_elb_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[3],
                           min_CGz_sp_elb_lb1 = dat_model_out$CGz_sp_elb_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[1],
                           min_CGz_sp_elb_lb2 = dat_model_out$CGz_sp_elb_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[3],
                           min_CGx_sp_man_lb1 = dat_model_out$CGx_sp_man_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[1],
                           min_CGx_sp_man_lb2 = dat_model_out$CGx_sp_man_lb+dat_model_out$CGx_sp_elbman_lb*man_fixed[3],
                           min_CGy_sp_man_lb1 = dat_model_out$CGy_sp_man_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[1],
                           min_CGy_sp_man_lb2 = dat_model_out$CGy_sp_man_lb+dat_model_out$CGy_sp_elbman_lb*man_fixed[3],
                           min_CGz_sp_man_lb1 = dat_model_out$CGz_sp_man_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[1],
                           min_CGz_sp_man_lb2 = dat_model_out$CGz_sp_man_lb+dat_model_out$CGz_sp_elbman_lb*man_fixed[3]),
                        by=list(species = dat_model_out$species), min)
max_ub     <- aggregate(list(max_CGx_sp_int_ub = dat_model_out$CGx_sp_int_ub,
                             max_CGy_sp_int_ub = dat_model_out$CGy_sp_int_ub,
                             max_CGz_sp_int_ub = dat_model_out$CGz_sp_int_ub,
                             max_CGx_sp_elb_ub1 = dat_model_out$CGx_sp_elb_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[1],
                             max_CGx_sp_elb_ub2 = dat_model_out$CGx_sp_elb_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[3],
                             max_CGy_sp_elb_ub1 = dat_model_out$CGy_sp_elb_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[1],
                             max_CGy_sp_elb_ub2 = dat_model_out$CGy_sp_elb_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[3],
                             max_CGz_sp_elb_ub1 = dat_model_out$CGz_sp_elb_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[1],
                             max_CGz_sp_elb_ub2 = dat_model_out$CGz_sp_elb_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[3],
                             max_CGx_sp_man_ub1 = dat_model_out$CGx_sp_man_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[1],
                             max_CGx_sp_man_ub2 = dat_model_out$CGx_sp_man_ub+dat_model_out$CGx_sp_elbman_ub*man_fixed[3],
                             max_CGy_sp_man_ub1 = dat_model_out$CGy_sp_man_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[1],
                             max_CGy_sp_man_ub2 = dat_model_out$CGy_sp_man_ub+dat_model_out$CGy_sp_elbman_ub*man_fixed[3],
                             max_CGz_sp_man_ub1 = dat_model_out$CGz_sp_man_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[1],
                             max_CGz_sp_man_ub2 = dat_model_out$CGz_sp_man_ub+dat_model_out$CGz_sp_elbman_ub*man_fixed[3]),
                        by=list(species = dat_model_out$species), max)
ci_bounds <- merge(min_lb,max_ub, by = c("species"))

dat_CG_plot     <- aggregate(list(mean_CGx_sp_int = dat_model_out$CGx_sp_int,
                                 mean_CGy_sp_int = dat_model_out$CGy_sp_int,
                                 mean_CGz_sp_int = dat_model_out$CGz_sp_int,
                                 mean_CGx_sp_elb = dat_model_out$CGx_sp_elb,
                                 mean_CGy_sp_elb = dat_model_out$CGy_sp_elb,
                                 mean_CGz_sp_elb = dat_model_out$CGz_sp_elb,
                                 mean_CGx_sp_man = dat_model_out$CGx_sp_man,
                                 mean_CGy_sp_man = dat_model_out$CGy_sp_man,
                                 mean_CGz_sp_man = dat_model_out$CGz_sp_man,
                                 mean_CGx_sp_elbman = dat_model_out$CGx_sp_elbman,
                                 mean_CGy_sp_elbman = dat_model_out$CGy_sp_elbman,
                                 mean_CGz_sp_elbman = dat_model_out$CGz_sp_elbman),  by=list(species = dat_model_out$species), mean)

test <- merge(aggregate(list(min_CGx_specific = dat_comp$min_CGx_specific,
                             min_CGy_specific = dat_comp$min_wing_CGy_specific,
                             min_CGz_specific = dat_comp$min_CGz_specific),  by=list(species = dat_comp$species), min),
              aggregate(list(max_CGx_specific = dat_comp$max_CGx_specific,
                             max_CGy_specific = dat_comp$max_wing_CGy_specific,
                             max_CGz_specific = dat_comp$max_CGz_specific),  by=list(species = dat_comp$species), max), by = "species")

CGx_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = test, aes(xmin = -min_CGx_specific, xmax = -max_CGx_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_CGx_specific = dat_comp$max_CGx_specific),  by=list(species = dat_comp$species), max),
             aes(x = -max_CGx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_CGx_specific = dat_comp$min_CGx_specific),  by=list(species = dat_comp$species), min),
             aes(x = -min_CGx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGx_specific = dat_comp$mean_CGx_specific),  by=list(species = dat_comp$species), mean),
             aes(x = -mean_CGx_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  geom_point(data = aggregate(list(shoulderx_specific = dat_bird$shoulderx_specific),  by=list(species = dat_bird$species), mean),
             aes(x = -shoulderx_specific, y = species), col = col_x, pch = 2, alpha = 1) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "x (% of full body length)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

CGy_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = test, aes(xmin = min_CGy_specific, xmax = max_CGy_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_CGy_specific = dat_comp$max_wing_CGy_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_CGy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_CGy_specific = dat_comp$min_wing_CGy_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_CGy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGy_specific = dat_comp$mean_wing_CGy_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_CGy_specific, y = species), fill = "white", col = col_x,pch = 21, alpha = 1) +
  geom_point(data = aggregate(list(shouldery = dat_final$shouldery_specific),  by=list(species = dat_final$species), mean),
             aes(x = shouldery, y = species), fill = col_x, col = "black", pch = 2, alpha = 1) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "y (% of half span)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

CGz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = 0, xmax = 1), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_point(data = aggregate(list(shoulderz_specific = dat_bird$shoulderz_specific),  by=list(species = dat_bird$species), mean),
             aes(x = shoulderz_specific, y = species), fill = col_x, col = "black", pch = 2, alpha = 1) +
  geom_errorbarh(data = test, aes(xmin = min_CGz_specific, xmax = max_CGz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_CGz_specific = dat_comp$max_CGz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_CGz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_CGz_specific = dat_comp$min_CGz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_CGz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_CGz_specific = dat_comp$mean_CGz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_CGz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = "z (% of max. body height)", limits= c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))


CGx_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1, xmax = 2.5), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.5, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = -pmin(min_CGx_sp_man_lb1,min_CGx_sp_man_lb2), xmax = -pmax(max_CGx_sp_man_ub1,max_CGx_sp_man_ub2), y = species),
                 col = col_z, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = -pmin(min_CGx_sp_elb_lb1,min_CGx_sp_elb_lb2), xmax = -pmax(max_CGx_sp_elb_ub1,max_CGx_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_man+mean_CGx_sp_elbman*elb_fixed[1]), y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_man+mean_CGx_sp_elbman*elb_fixed[2]), y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_man+mean_CGx_sp_elbman*elb_fixed[3]), y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_elb+mean_CGx_sp_elbman*man_fixed[1]), y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_elb+mean_CGx_sp_elbman*man_fixed[2]), y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = -(mean_CGx_sp_elb+mean_CGx_sp_elbman*man_fixed[3]), y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
    #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  # to get %/deg we need to divide the output slopes by 10 instead I will just adjust the labels accordingly
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-1,2.5) ,breaks = c(-1,0,1,2), labels = c(-0.1,0,0.1,0.2))+
  geom_rangeframe() +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))


CGy_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1, xmax = 2.5), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.50, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGy_sp_man_lb1,min_CGy_sp_man_lb2), xmax = pmax(max_CGy_sp_man_ub1,max_CGy_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGy_sp_elb_lb1,min_CGy_sp_elb_lb2), xmax = pmax(max_CGy_sp_elb_ub1,max_CGy_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  # add data
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_man+mean_CGy_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_man+mean_CGy_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_man+mean_CGy_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_elb+mean_CGy_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_elb+mean_CGy_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGy_sp_elb+mean_CGy_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  # to get %/deg we need to divide the output slopes by 10 instead I will just adjust the labels accordingly
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-1,2.5) ,breaks = c(-1,0,1,2), labels = c(-0.1,0,0.1,0.2))+
  geom_rangeframe() +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))


CGz_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -1, xmax = 2.4), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.50, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGz_sp_man_lb1,min_CGz_sp_man_lb2), xmax = pmax(max_CGz_sp_man_ub1,max_CGz_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_bounds, aes(xmin = pmin(min_CGz_sp_elb_lb1,min_CGz_sp_elb_lb2), xmax = pmax(max_CGz_sp_elb_ub1,max_CGz_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  # add data
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_man+mean_CGz_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_man+mean_CGz_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_man+mean_CGz_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_elb+mean_CGz_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_elb+mean_CGz_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_CG_plot, aes(x = mean_CGz_sp_elb+mean_CGz_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  # to get %/deg we need to divide the output slopes by 10 instead I will just adjust the labels accordingly
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-1,2.4) ,breaks = c(-1,0,1,2), labels = c(-0.1,0,0.1,0.2))+
  geom_rangeframe() +
  annotate(geom = "segment", x = -1, xend = 2, y = log(0), yend = log(0))

effect_CGx <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(CGx_elb_etap = dat_model_out$CGx_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGx_man_etap = dat_model_out$CGx_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGx_elbman_etap = dat_model_out$CGx_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGx_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_CGy <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(CGy_elb_etap = dat_model_out$CGy_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGy_man_etap = dat_model_out$CGy_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGy_elbman_etap = dat_model_out$CGy_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGy_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_CGz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(CGz_elb_etap = dat_model_out$CGz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGz_man_etap = dat_model_out$CGz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(CGz_elbman_etap = dat_model_out$CGz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = CGz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_Ixx <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Ixx_elb_etap = dat_model_out$Ixx_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixx_man_etap = dat_model_out$Ixx_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixx_elbman_etap = dat_model_out$Ixx_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixx_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_Iyy <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Iyy_elb_etap = dat_model_out$Iyy_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Iyy_man_etap = dat_model_out$Iyy_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Iyy_elbman_etap = dat_model_out$Iyy_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Iyy_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_Izz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Izz_elb_etap = dat_model_out$Izz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Izz_man_etap = dat_model_out$Izz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Izz_elbman_etap = dat_model_out$Izz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Izz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

effect_Ixz <- ggplot()+
  # add background info
  geom_density(data = aggregate(list(Ixz_elb_etap = dat_model_out$Ixz_elb_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_elb_etap), fill = col_elb, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixz_man_etap = dat_model_out$Ixz_man_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_man_etap), fill = col_man, alpha = 0.5)  +
  geom_density(data = aggregate(list(Ixz_elbman_etap = dat_model_out$Ixz_elbman_etap),  by=list(species = dat_model_out$species), mean),
               aes(x = Ixz_elbman_etap), fill = col_elbman, alpha = 0.5)  +
  # theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  # axis control
  scale_y_continuous(limits = c(0,12), name = "") +
  scale_x_continuous(limits = c(0,1), name = "Cohen's effect size") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

# -------- Combine panels into figure ------------
bottomrow <- plot_grid(phylo_plot,CGx_specific,CGy_specific, CGz_specific,
                       #arrangement data
                       ncol = 4,
                       rel_widths = c(1.5,1,1,1),
                       #labels
                       labels = c("A","B","C","D"),
                       label_size = 10,
                       label_fontfamily = "sans")

#exported as 6x12
figure_final <- plot_grid(toprow,bottomrow,
                           #arrangement data
                           ncol = 1, nrow = 2, rel_heights = c(1,1))
#expot as 6x10
bottomrow <- plot_grid(phylo_plot,CGx_angleeffect,CGy_angleeffect,
                       #arrangement data
                       ncol = 3,
                       rel_widths = c(1.5,1,1),
                       #labels
                       labels = c("A","B","C"),
                       label_size = 10,
                       label_fontfamily = "sans")

bottomrow <- plot_grid(cohens_effect_x,cohens_effect_y,cohens_effect_z,
                       #arrangement data
                       nrow = 3,
                       rel_heights = c(1,1,1),
                       #labels
                       labels = c("A","B","C"),
                       label_size = 10,
                       label_fontfamily = "sans")

effectsizes <- plot_grid(effect_CGx,effect_CGy,effect_CGz, effect_Ixx,effect_Iyy,effect_Izz,effect_Ixz,
                       #arrangement data
                       nrow = 7,
                       #labels
                       labels = c("CGx","CGy","CGz","Ixx","Iyy","Izz","Ixz"),
                       label_size = 10,
                       label_fontfamily = "sans")


## ----------------- Plot of linear model outputs ------------------
# order species to match the phylogeny
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
col_xz = "#af0069"

min_I_lb   <- aggregate(list(min_Ixx_sp_int_lb = dat_model_out$Ixx_sp_int_lb,
                             min_Iyy_sp_int_lb = dat_model_out$Iyy_sp_int_lb,
                             min_Izz_sp_int_lb = dat_model_out$Izz_sp_int_lb,
                             min_Ixz_sp_int_lb = dat_model_out$Ixz_sp_int_lb,
                             min_Ixx_sp_elb_lb1 = dat_model_out$Ixx_sp_elb_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[1],
                             min_Ixx_sp_elb_lb2 = dat_model_out$Ixx_sp_elb_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[3],
                             min_Iyy_sp_elb_lb1 = dat_model_out$Iyy_sp_elb_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[1],
                             min_Iyy_sp_elb_lb2 = dat_model_out$Iyy_sp_elb_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[3],
                             min_Izz_sp_elb_lb1 = dat_model_out$Izz_sp_elb_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[1],
                             min_Izz_sp_elb_lb2 = dat_model_out$Izz_sp_elb_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[3],
                             min_Ixz_sp_elb_lb1 = dat_model_out$Ixz_sp_elb_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[1],
                             min_Ixz_sp_elb_lb2 = dat_model_out$Ixz_sp_elb_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[3],
                             min_Ixx_sp_man_lb1 = dat_model_out$Ixx_sp_man_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[1],
                             min_Ixx_sp_man_lb2 = dat_model_out$Ixx_sp_man_lb+dat_model_out$Ixx_sp_elbman_lb*man_fixed[3],
                             min_Iyy_sp_man_lb1 = dat_model_out$Iyy_sp_man_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[1],
                             min_Iyy_sp_man_lb2 = dat_model_out$Iyy_sp_man_lb+dat_model_out$Iyy_sp_elbman_lb*man_fixed[3],
                             min_Izz_sp_man_lb1 = dat_model_out$Izz_sp_man_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[1],
                             min_Izz_sp_man_lb2 = dat_model_out$Izz_sp_man_lb+dat_model_out$Izz_sp_elbman_lb*man_fixed[3],
                             min_Ixz_sp_man_lb1 = dat_model_out$Ixz_sp_man_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[1],
                             min_Ixz_sp_man_lb2 = dat_model_out$Ixz_sp_man_lb+dat_model_out$Ixz_sp_elbman_lb*man_fixed[3]),
                        by=list(species = dat_model_out$species), min)
max_I_ub   <- aggregate(list(max_Ixx_sp_int_ub = dat_model_out$Ixx_sp_int_ub,
                             max_Iyy_sp_int_ub = dat_model_out$Iyy_sp_int_ub,
                             max_Izz_sp_int_ub = dat_model_out$Izz_sp_int_ub,
                             max_Ixz_sp_int_ub = dat_model_out$Ixz_sp_int_ub,
                             max_Ixx_sp_elb_ub1 = dat_model_out$Ixx_sp_elb_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[1],
                             max_Ixx_sp_elb_ub2 = dat_model_out$Ixx_sp_elb_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[3],
                             max_Iyy_sp_elb_ub1 = dat_model_out$Iyy_sp_elb_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[1],
                             max_Iyy_sp_elb_ub2 = dat_model_out$Iyy_sp_elb_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[3],
                             max_Izz_sp_elb_ub1 = dat_model_out$Izz_sp_elb_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[1],
                             max_Izz_sp_elb_ub2 = dat_model_out$Izz_sp_elb_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[3],
                             max_Ixz_sp_elb_ub1 = dat_model_out$Ixz_sp_elb_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[1],
                             max_Ixz_sp_elb_ub2 = dat_model_out$Ixz_sp_elb_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[3],
                             max_Ixx_sp_man_ub1 = dat_model_out$Ixx_sp_man_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[1],
                             max_Ixx_sp_man_ub2 = dat_model_out$Ixx_sp_man_ub+dat_model_out$Ixx_sp_elbman_ub*man_fixed[3],
                             max_Iyy_sp_man_ub1 = dat_model_out$Iyy_sp_man_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[1],
                             max_Iyy_sp_man_ub2 = dat_model_out$Iyy_sp_man_ub+dat_model_out$Iyy_sp_elbman_ub*man_fixed[3],
                             max_Izz_sp_man_ub1 = dat_model_out$Izz_sp_man_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[1],
                             max_Izz_sp_man_ub2 = dat_model_out$Izz_sp_man_ub+dat_model_out$Izz_sp_elbman_ub*man_fixed[3],
                             max_Ixz_sp_man_ub1 = dat_model_out$Ixz_sp_man_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[1],
                             max_Ixz_sp_man_ub2 = dat_model_out$Ixz_sp_man_ub+dat_model_out$Ixz_sp_elbman_ub*man_fixed[3]),
                        by=list(species = dat_model_out$species), max)
ci_I_bounds <- merge(min_I_lb,max_I_ub, by = c("species"))

dat_I_plot     <- aggregate(list(mean_Ixx_sp_int = dat_model_out$Ixx_sp_int,
                                 mean_Iyy_sp_int = dat_model_out$Iyy_sp_int,
                                 mean_Izz_sp_int = dat_model_out$Izz_sp_int,
                                 mean_Ixx_sp_elb = dat_model_out$Ixx_sp_elb,
                                 mean_Iyy_sp_elb = dat_model_out$Iyy_sp_elb,
                                 mean_Izz_sp_elb = dat_model_out$Izz_sp_elb,
                                 mean_Ixz_sp_elb = dat_model_out$Ixz_sp_elb,
                                 mean_Ixx_sp_man = dat_model_out$Ixx_sp_man,
                                 mean_Iyy_sp_man = dat_model_out$Iyy_sp_man,
                                 mean_Izz_sp_man = dat_model_out$Izz_sp_man,
                                 mean_Ixz_sp_man = dat_model_out$Ixz_sp_man,
                                 mean_Ixx_sp_elbman = dat_model_out$Ixx_sp_elbman,
                                 mean_Iyy_sp_elbman = dat_model_out$Iyy_sp_elbman,
                                 mean_Izz_sp_elbman = dat_model_out$Izz_sp_elbman,
                                 mean_Ixz_sp_elbman = dat_model_out$Ixz_sp_elbman),  by=list(species = dat_model_out$species), mean)


dat_I_sp_plot <- merge(aggregate(list(min_Ixx_specific = dat_comp$min_Ixx_specific,
                             min_Iyy_specific = dat_comp$min_Iyy_specific,
                             min_Izz_specific = dat_comp$min_Izz_specific,
                             min_Ixz_specific = dat_comp$min_Ixz_specific),  by=list(species = dat_comp$species), min),
              aggregate(list(max_Ixx_specific = dat_comp$max_Ixx_specific,
                             max_Iyy_specific = dat_comp$max_Iyy_specific,
                             max_Izz_specific = dat_comp$max_Izz_specific,
                             max_Ixz_specific = dat_comp$max_Ixz_specific),  by=list(species = dat_comp$species), max), by = "species")

Ixx_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Ixx_specific, xmax = max_Ixx_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Ixx_specific = dat_comp$max_Ixx_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Ixx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Ixx_specific = dat_comp$min_Ixx_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Ixx_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Ixx_specific = dat_comp$mean_Ixx_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Ixx_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["xx"]," (% of maximum inertia)", sep = "")), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))


Iyy_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Iyy_specific, xmax = max_Iyy_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Iyy_specific = dat_comp$max_Iyy_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Iyy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Iyy_specific = dat_comp$min_Iyy_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Iyy_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Iyy_specific = dat_comp$mean_Iyy_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Iyy_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["yy"]," (% of maximum inertia)", sep = "")), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Izz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Izz_specific, xmax = max_Izz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Izz_specific = dat_comp$max_Izz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Izz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Izz_specific = dat_comp$min_Izz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Izz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Izz_specific = dat_comp$mean_Izz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Izz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["zz"]," (% of maximum inertia)", sep = "")), limits= c(0,0.03), breaks = c(0,0.01,0.02,0.03), labels = c(0,1,2,3)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.03, y = log(0), yend = log(0))

Ixz_specific <- ggplot() +
  # add background info
  geom_rect(data = shading, aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  # add data
  geom_errorbarh(data = dat_I_sp_plot, aes(xmin = min_Ixz_specific, xmax = max_Ixz_specific, y = species), height = 0, alpha = 0.5)+
  geom_point(data = aggregate(list(max_Ixz_specific = dat_comp$max_Ixz_specific),  by=list(species = dat_comp$species), max),
             aes(x = max_Ixz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(min_Ixz_specific = dat_comp$min_Ixz_specific),  by=list(species = dat_comp$species), min),
             aes(x = min_Ixz_specific, y = species), fill = col_x, col = "black",pch = 22, alpha = 0.85) +
  geom_point(data = aggregate(list(mean_Ixz_specific = dat_comp$mean_Ixz_specific),  by=list(species = dat_comp$species), mean),
             aes(x = mean_Ixz_specific, y = species), fill = "white", col = col_x, pch = 21, alpha = 1, size = 1.5) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "") +
  scale_x_continuous(name = expression(paste("I"["xz"]," (% of maximum inertia)", sep = "")), limits= c(-0.015,0.015), breaks = c(-0.015,0,0.015), labels = c(-1.5,0,1.5)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.015, xend = 0.015, y = log(0), yend = log(0))

Ixx_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars

  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Ixx_sp_man_lb1,min_Ixx_sp_man_lb2), xmax = pmax(max_Ixx_sp_man_ub1,max_Ixx_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Ixx_sp_elb_lb1,min_Ixx_sp_elb_lb2), xmax = pmax(max_Ixx_sp_elb_ub1,max_Ixx_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_man+mean_Ixx_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_man+mean_Ixx_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_man+mean_Ixx_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_elb+mean_Ixx_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_elb+mean_Ixx_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Ixx_sp_elb+mean_Ixx_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
    #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-0.05,0.15), breaks = c(-0.05,0,0.05,0.1,0.15), labels = c(-0.005,0,0.005,0.01,0.015)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.05, xend = 0.15, y = log(0), yend = log(0))


Iyy_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars

  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Iyy_sp_man_lb1,min_Iyy_sp_man_lb2), xmax = pmax(max_Iyy_sp_man_ub1,max_Iyy_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Iyy_sp_elb_lb1,min_Iyy_sp_elb_lb2), xmax = pmax(max_Iyy_sp_elb_ub1,max_Iyy_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_man+mean_Iyy_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_man+mean_Iyy_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_man+mean_Iyy_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_elb+mean_Iyy_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_elb+mean_Iyy_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Iyy_sp_elb+mean_Iyy_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-0.05,0.15), breaks = c(-0.05,0,0.05,0.1,0.15), labels = c(-0.005,0,0.005,0.01,0.015)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.05, xend = 0.15, y = log(0), yend = log(0))



Izz_angleeffect <- ggplot()+
  # add background info
  geom_rect(data = shading,aes(ymin = col1, ymax = col2, xmin = -Inf, xmax = Inf), alpha = 0.1, position = position_nudge(y = -0.5)) +
  geom_segment(aes(y = 0, yend = 23, x = 0, xend = 0), alpha = 0.25, linetype = 2) +
  # add error bars
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Izz_sp_man_lb1,min_Izz_sp_man_lb2), xmax = pmax(max_Izz_sp_man_ub1,max_Izz_sp_man_ub2), y = species),
                 col = col_xz, position = position_nudge(y = 0.2), alpha = 0.2, height = 0.3, size = 1) +
  geom_errorbarh(data = ci_I_bounds, aes(xmin = pmin(min_Izz_sp_elb_lb1,min_Izz_sp_elb_lb2), xmax = pmax(max_Izz_sp_elb_ub1,max_Izz_sp_elb_ub2), y = species),
                 col = col_y, position = position_nudge(y = -0.2), alpha = 0.2, height = 0.3, size = 1) +
  #add data
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_man+mean_Izz_sp_elbman*elb_fixed[1], y = species), col = col_xz, alpha = 0.2, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_man+mean_Izz_sp_elbman*elb_fixed[2], y = species), col = col_xz, alpha = 0.5, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_man+mean_Izz_sp_elbman*elb_fixed[3], y = species), col = col_xz, alpha = 1, position = position_nudge(y = 0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_elb+mean_Izz_sp_elbman*man_fixed[1], y = species), col = col_y, alpha = 0.4, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_elb+mean_Izz_sp_elbman*man_fixed[2], y = species), col = col_y, alpha = 0.7, position = position_nudge(y = -0.2)) +
  geom_point(data = dat_I_plot, aes(x = mean_Izz_sp_elb+mean_Izz_sp_elbman*man_fixed[3], y = species), col = col_y, alpha = 1, position = position_nudge(y = -0.2)) +
  #theme control
  th +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) +
  #axis control
  scale_y_discrete(expand = c(0.1,0), limits = phylo_order$species, name = "")+
  # need to adjust
  scale_x_continuous(name = "Coefficient of Regression (%/°)", limits= c(-0.05,0.15), breaks = c(-0.05,0,0.05,0.1,0.15), labels = c(-0.005,0,0.005,0.01,0.015)) +
  geom_rangeframe() +
  annotate(geom = "segment", x = -0.05, xend = 0.15, y = log(0), yend = log(0))


# -------- Combine panels into figure ------------
#exported as 6x14
bottomrow <- plot_grid(phylo_plot,Ixx_specific,Iyy_specific, Izz_specific,Ixz_specific,
                       #arrangement data
                       ncol = 5,
                       rel_widths = c(1.5,1,1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")

bottomrow <- plot_grid(phylo_plot,Ixx_angleeffect,Iyy_angleeffect, Izz_angleeffect,
                       #arrangement data
                       ncol = 4,
                       rel_widths = c(1.5,1,1,1),
                       #labels
                       labels = c("","","",""),
                       label_size = 10,
                       label_fontfamily = "sans")

#exported as 7.5x9
figure_final <- plot_grid(toprow,bottomrow,
                          #arrangement data
                          ncol = 1, nrow = 2, rel_heights = c(1,1))


Ixx_iso_model <- lm(log(max_Ixx) ~ 1 + offset((5/3)*log(full_m)), data = dat_comp)
Iyy_iso_model <- lm(log(max_Iyy) ~ 1 + offset((5/3)*log(full_m)), data = dat_comp)
Izz_iso_model <- lm(log(max_Izz) ~ 1 + offset((5/3)*log(full_m)), data = dat_comp)

I_isometry <- ggplot() +
  geom_point(data = dat_comp, aes (x = full_m, y = max_Ixx), col = col_x) +
  geom_point(data = dat_comp, aes (x = full_m, y = max_Iyy), col = "#479030")+
  geom_point(data = dat_comp, aes (x = full_m, y = max_Izz), col = col_z)+
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(Ixx_iso_model)[1])*full_m^(5/3)), col = col_x, linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(Iyy_iso_model)[1])*full_m^(5/3)), col = "#479030", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(Izz_iso_model)[1])*full_m^(5/3)), col = col_z, linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Maximum moment of inertia (kg-m"^2,")", sep = "")),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)



# export as 5x5
CGx_iso_model <- lm(log(-min_CGx) ~ 1 + offset((1/3)*log(full_m)), data = dat_comp)
CGy_iso_model <- lm(log(max_wing_CGy) ~ 1 + offset((1/3)*log(full_m)), data = dat_comp)
CGz_iso_model <- lm(log(max_CGz) ~ 1 + offset((1/3)*log(full_m)), data = dat_comp)



CG_isometry <- ggplot() +
  geom_point(data = dat_comp, aes (x = full_m, y = -min_CGx), col = col_x)+
  geom_point(data = dat_comp, aes (x = full_m, y = max_wing_CGy), col = "#479030")+
  geom_point(data = dat_comp, aes (x = full_m, y = max_CGz), col = col_z)+
  # verify the fits match PGLSS model
  geom_line(data = test, aes(x = full_m, y = exp(fit)), col = col_x) +
  geom_ribbon(data = test,
              aes(x = full_m, ymin = exp(lwr), ymax = exp(upr)),
              col = col_x, fill = col_x, alpha = 0.4) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(CGx_iso_model)[1])*full_m^(1/3)), col = col_x, linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(CGy_iso_model)[1])*full_m^(1/3)), col = "#479030", linetype = 2) +
  geom_line(data = dat_comp, aes(x = full_m, y = exp(coef(CGz_iso_model)[1])*full_m^(1/3)), col = col_z, linetype = 2) +
  th +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Center of gravity (m)", sep = "")),
                     breaks = c(1E-2,1E-1,1E0), labels = c(expression(10^-2),expression(10^-1),expression(10^0)))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-2, yend = 1E0) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

bottomrow <- plot_grid(CG_isometry,I_isometry,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(1,1),
                       #labels
                       labels = c("",""),
                       label_size = 10,
                       label_fontfamily = "sans")

### ----------- Gull only data ---------------
gull_model_Ixx <- lm(full_Ixx~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_model_Iyy <- lm(full_Iyy~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_model_Izz <- lm(full_Izz~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_model_Ixz <- lm(full_Ixz~elbow*manus, data = subset(dat_final,species == "lar_gla"))
gull_Ixx <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Ixx)) + th

gull_Iyy <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Iyy)) + th

gull_Izz <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Izz)) + th

gull_Ixz <- ggplot()+
  geom_point(data = subset(dat_final,species == "lar_gla"), aes(x = elbow, y = manus, col = full_Ixz)) + th

## ----------------- Plots for sharing with others -----------------
# Species rainbow colours
cc_rain1  <- scales::seq_gradient_pal("#DA9101", "#CA302F", "Lab")(seq(0,1,length.out=6))
cc_rain2  <- scales::seq_gradient_pal("#FCC201", "#DA9101", "Lab")(seq(0,1,length.out=3))
cc_rain3  <- scales::seq_gradient_pal("#138715", "#FCC201", "Lab")(seq(0,1,length.out=6))
cc_rain4  <- scales::seq_gradient_pal("#1FC3CD", "#138715", "Lab")(seq(0,1,length.out=4))
cc_rain5  <- scales::seq_gradient_pal("#304CC8", "#1FC3CD", "Lab")(seq(0,1,length.out=5))
cc_rain6  <- scales::seq_gradient_pal("#5F2CC8", "#304CC8", "Lab")(seq(0,1,length.out=3))
cc_rain   <- c(cc_rain6, cc_rain5[2:5], cc_rain4[2:4], cc_rain3[2:6],cc_rain2[2:3],cc_rain1[2:6])

validation_Ixx_plot <- ggplot()+
  geom_errorbar(data = dat_comp, aes(x = full_m, ymax = max_wing_Ixx, ymin = sachs_pred_Ixx, col = species_order), alpha = 0.5) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "gray") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray", linetype = 2) + # Berg and Rayner, 1995
  geom_point(data = dat_comp, aes(x = full_m, y = sachs_pred_Ixx, fill = species_order, col = species_order), pch = 22, size = 1.5) + # Sachs, 2005
  geom_point(data = dat_comp, aes(x = full_m, y = max_wing_Ixx, fill = species_order), pch = 21, col = "black", size = 2)+
  th  +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  scale_colour_manual(values = rev(cc_rain), name = "species") +
  scale_x_continuous(trans='log10', name = "Mass (kg)",
                     breaks = c(0.01,0.1,1,10), labels = c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1)))+
  scale_y_continuous(trans='log10', name = expression(paste("Maximum wing I"["xx"]," (kg-m"^2,")", sep = "")),
                     breaks = c(1E-6,1E-5,1E-4,1E-3,1E-2,1E-1), labels = c(expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-6, yend = 1E-1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

chulls_xz <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                   function(df) df[chull(df$elbow, df$manus), ])

ROM_plot <- ggplot()+
  geom_polygon(data = chulls_xz, aes(x = elbow, y = manus, fill = species_order), col = NA) +
  #stat_contour_filled(data = data.fit, aes(x = elbow, y = manus, z = prop_q, colour = ..level..),
  #             breaks = quantile(data.fit$prop_q, seq(0, 1, 0.05)), size = 0.8) +
  facet_wrap(~species_order, nrow = 2) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.line =  element_line(colour = "black", linetype=1))+
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)') +
  scale_y_continuous(name='Wrist angle (°)') + theme(legend.position="none")

validation_Izz_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = wing_Izz, col = clade)) +
  geom_line(data = unique(dat_final[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_final[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_final[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Ixx_body_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = full_Ixx, col = clade)) +
  geom_point(aes(x = 0.0875, y = 928/(100^2*1000)), col = "black") + # Hedrick and Biewener
  geom_point(aes(x = 0.2893, y = 12889/(100^2*1000)), col = "black") + # Hedrick and Biewener
  th +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

q_plot <- ggplot() +
  geom_point(data = aggregate(list(prop_q_dot = dat_final$prop_q_dot, full_m = dat_final$full_m),  by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), max),
             aes(x = full_m, y = prop_q_dot, fill = species), col = "black", alpha = 0.8, pch = 21) +
  geom_line(data = unique(dat_final[,c("species_order","BirdID","full_m")]), aes(x = full_m, y = 0.66*full_m^-0.35),col = "gray") + # Kirkpatrick, 1994
  #geom_point(data = aggregate(list(prop_q_dot = dat_final$prop_q_dot, full_m = dat_final$full_m),  by=list(species = dat_final$species_order, BirdID = dat_final$BirdID), min),
  #           aes(x = full_m, y = prop_q_dot, col = species)) +
  th +
  scale_x_continuous(trans='log10', limits = c(0.01,10), breaks = c(0.01,0.1,1,10)) +
  scale_y_continuous(trans='log10', limits = c(0.1,10), breaks = c(0.1,1,10))+
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0, y = 1E-1, yend = 1E1) +
  annotate(geom = "segment", x = 0.01, xend = 10, y = 0, yend = 0)

test_plot <- ggplot() +
  geom_point(aes(x = dat_final$elbow[which(dat_final$species %in% test$species & dat_final$prop_q_dot %in% test$prop_q_dot)],
                 y = dat_final$manus[which(dat_final$species %in% test$species & dat_final$prop_q_dot %in% test$prop_q_dot)],
                 col = log(dat_final$prop_q_dot[which(dat_final$species %in% test$species & dat_final$prop_q_dot %in% test$prop_q_dot)]))) + th

ROM_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = elbow, y = manus, col = species_order), alpha = 0.5) +
  facet_wrap(~species_order, nrow = 2) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = ggplot2::element_rect(fill = "transparent"),
    # Background behind actual data points
    plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
    axis.line =  element_line(colour = "black", linetype=1))+
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)') +
  scale_y_continuous(name='Wrist angle (°)')

# Species rainbow colours
cc_rain1  <- scales::seq_gradient_pal("#DA9101", "#CA302F", "Lab")(seq(0,1,length.out=6))
cc_rain2  <- scales::seq_gradient_pal("#FCC201", "#DA9101", "Lab")(seq(0,1,length.out=3))
cc_rain3  <- scales::seq_gradient_pal("#138715", "#FCC201", "Lab")(seq(0,1,length.out=6))
cc_rain4  <- scales::seq_gradient_pal("#304CC8", "#138715", "Lab")(seq(0,1,length.out=6))
cc_rain5  <- scales::seq_gradient_pal("#4E2388", "#304CC8", "Lab")(seq(0,1,length.out=5))
cc_rain   <- c(cc_rain5, cc_rain4[2:6], cc_rain3[2:6],cc_rain2[2:3],cc_rain1[2:6])

chulls_xz <- ddply(dat_final[,c("species_order","BirdID","elbow","manus")], .(species_order, BirdID),
                   function(df) df[chull(df$elbow, df$manus), ])

ROM_plot <- ggplot()+
  geom_polygon(data = chulls_xz, aes(x = elbow, y = manus, fill = species_order), col = NA, alpha = 0.5) +
  #stat_contour_filled(data = data.fit, aes(x = elbow, y = manus, z = prop_q, colour = ..level..),
  #             breaks = quantile(data.fit$prop_q, seq(0, 1, 0.05)), size = 0.8) +
  facet_wrap(~species_order, nrow = 2) +
  # colour control
  scale_fill_manual(values = rev(cc_rain), name = "species") +
  theme(strip.background = element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        # Background behind actual data points
        plot.background  = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.line =  element_line(colour = "black", linetype=1))+
  coord_fixed() +
  scale_x_continuous(name='Elbow angle (°)') +
  scale_y_continuous(name='Wrist angle (°)') + theme(legend.position="none")


test <- lm(prop_q_dot_nd ~ elbow*manus + species_order, data = dat_final)
xgrid <-  seq(floor(min(dat_final$elbow)), ceiling(max(dat_final$elbow)), 1)
ygrid <-  seq(floor(min(dat_final$manus)), ceiling(max(dat_final$manus)), 1)
zgrid <-  phylo_order$species
data.fit       <- expand.grid(elbow = xgrid, manus = ygrid, species_order = zgrid)
data.fit$prop_q  <-  predict(test, newdata = data.fit)

# Visualize the wings as required - For each calibration verify that the axis is RH
m = 1:nrow(dat_wing_curr)
max = 0.4
plot(dat_wing_curr$pt2_Y[m],dat_wing_curr$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_wing_curr$pt3_Y[m],dat_wing_curr$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_wing_curr$pt4_Y[m], dat_wing_curr$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_wing_curr$pt1_Y[m],dat_wing_curr$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_wing_curr$pt8_Y[m],dat_wing_curr$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_wing_curr$pt10_Y[m],dat_wing_curr$pt10_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

plot(dat_clean$pt2_Y[m],dat_clean$pt2_Z[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_clean$pt3_Y[m],dat_clean$pt3_Z[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_clean$pt4_Y[m], dat_clean$pt4_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_clean$pt8_Y[m],dat_clean$pt8_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_clean$pt11_Y[m],dat_clean$pt11_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")
points(dat_clean$pt1_Y[m],dat_clean$pt1_Z[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")

plot(dat_clean$pt2_Z[m],dat_clean$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_clean$pt3_Z[m],dat_clean$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_clean$pt4_Z[m], dat_clean$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_clean$pt1_Z[m],dat_clean$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_clean$pt8_Z[m],dat_clean$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_clean$pt10_Z[m],dat_clean$pt10_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

x = c(dat_clean$pt1_X[m],dat_clean$pt2_X[m],dat_clean$pt3_X[m],dat_clean$pt4_X[m],
      dat_clean$pt6_X[m],dat_clean$pt7_X[m],dat_clean$pt8_X[m],dat_clean$pt9_X[m],
      dat_clean$pt10_X[m],dat_clean$pt11_X[m],dat_clean$pt12_X[m])
y = c(dat_clean$pt1_Y[m],dat_clean$pt2_Y[m],dat_clean$pt3_Y[m],dat_clean$pt4_Y[m],
      dat_clean$pt6_Y[m],dat_clean$pt7_Y[m],dat_clean$pt8_Y[m],dat_clean$pt9_Y[m],
      dat_clean$pt10_Y[m],dat_clean$pt11_Y[m],dat_clean$pt12_Y[m])
z = c(dat_clean$pt1_Z[m],dat_clean$pt2_Z[m],dat_clean$pt3_Z[m],dat_clean$pt4_Z[m],
      dat_clean$pt6_Z[m],dat_clean$pt7_Z[m],dat_clean$pt8_Z[m],dat_clean$pt9_Z[m],
      dat_clean$pt10_Z[m],dat_clean$pt11_Z[m],dat_clean$pt12_Z[m])
plot3d(x, y, z, col = c("black","gray20","gray50","gray70","red", "orange", "yellow", "green", "blue", "navy", "purple"))

plot(dat_clean$elbow,dat_clean$manus)

plot(dat_complete$pt2_Y[m],dat_complete$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_complete$pt3_Y[m],dat_complete$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_complete$pt4_Y[m], dat_complete$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_complete$pt1_Y[m],dat_complete$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_complete$pt8_Y[m],dat_complete$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_complete$pt9_Y[m],dat_complete$pt9_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

points(tmp$pt2_Y[m],tmp$pt2_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(tmp$pt3_Y[m],tmp$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(tmp$pt4_Y[m], tmp$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(tmp$pt1_Y[m],tmp$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(tmp$pt8_Y[m],tmp$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(tmp$pt9_Y[m],tmp$pt9_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")


# ## TO CHECK DATA
# x_vertices = c(dat_wing$pt1_X[i], dat_wing$pt6_X[i],dat_wing$pt7_X[i],dat_wing$pt8_X[i],dat_wing$pt9_X[i],dat_wing$pt10_X[i],dat_wing$pt11_X[i],dat_wing$pt12_X[i])
# y_vertices = c(dat_wing$pt1_Y[i], dat_wing$pt6_Y[i],dat_wing$pt7_Y[i],dat_wing$pt8_Y[i],dat_wing$pt9_Y[i],dat_wing$pt10_Y[i],dat_wing$pt11_Y[i],dat_wing$pt12_Y[i])
# z_vertices = c(dat_wing$pt1_Z[i], dat_wing$pt6_Z[i],dat_wing$pt7_Z[i],dat_wing$pt8_Z[i],dat_wing$pt9_Z[i],dat_wing$pt10_Z[i],dat_wing$pt11_Z[i],dat_wing$pt12_Z[i])
#
#
# x = c(dat_wing$pt1_X[m],dat_wing$pt2_X[m],dat_wing$pt3_X[m],dat_wing$pt4_X[m],
#       dat_wing$pt6_X[m],dat_wing$pt7_X[m],dat_wing$pt8_X[m],dat_wing$pt9_X[m],
#       dat_wing$pt10_X[m],dat_wing$pt11_X[m],dat_wing$pt12_X[m])
# y = c(dat_wing$pt1_Y[m],dat_wing$pt2_Y[m],dat_wing$pt3_Y[m],dat_wing$pt4_Y[m],
#       dat_wing$pt6_Y[m],dat_wing$pt7_Y[m],dat_wing$pt8_Y[m],dat_wing$pt9_Y[m],
#       dat_wing$pt10_Y[m],dat_wing$pt11_Y[m],dat_wing$pt12_Y[m])
# z = c(dat_wing$pt1_Z[m],dat_wing$pt2_Z[m],dat_wing$pt3_Z[m],dat_wing$pt4_Z[m],
#       dat_wing$pt6_Z[m],dat_wing$pt7_Z[m],dat_wing$pt8_Z[m],dat_wing$pt9_Z[m],
#       dat_wing$pt10_Z[m],dat_wing$pt11_Z[m],dat_wing$pt12_Z[m])
# plot3d(x, y, z, col = c("black","gray20","gray50","gray70","red", "orange", "yellow", "green", "blue", "navy", "purple"))




chulls_xz <- ddply(dat_final[,c("species","BirdID","full_CGx_specific","full_CGz_specific")], .(species, BirdID),
                   function(df) df[chull(df$full_CGx_specific, df$full_CGz_specific), ])

CG_side_view <- ggplot()+
  geom_polygon(data=chulls_xz, aes(x=-full_CGx_specific, y=-full_CGz_specific, group=interaction(species,BirdID), fill = species), alpha = 0.5) +
  geom_point(data=dat_comp, aes(x=-mean_CGx_specific, y=-mean_CGz_specific, col = species)) +

  scale_y_continuous(name = "z (% of full body length)",limits = c(-0.11,0.11), breaks = c(-0.1,0,0.1), labels = c(-10,0,10)) +
  scale_x_continuous(name = "x (% of full body length)",limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  coord_fixed() + th +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.1, yend = 0.1) +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))

chulls_xy <- ddply(dat_final[,c("species","BirdID","wing_CGx_specific","wing_CGy_specific")], .(species, BirdID),
                   function(df) df[chull(df$wing_CGx_specific, df$wing_CGy_specific), ])

CG_top_view <- ggplot()+
  geom_polygon(data=chulls_xy, aes(x=wing_CGy_specific, y=wing_CGx_specific, group=interaction(species,BirdID), fill = species), alpha = 0.15) +
  geom_point(data=dat_comp, aes(x=mean_wing_CGy_specific, y=mean_wing_CGx_specific, col = species)) +
  scale_y_continuous(name = "x (% of half wing span)",limits = c(-0.20,0.05), breaks = c(-0.2,-0.1,0), labels = c(20,10,0)) +
  scale_x_continuous(name = "y (% of half wing span)",limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  coord_fixed() + th +
  geom_rangeframe() +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -0.2, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0))





q_dot_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = wing_CGy_specific, col = species), pch = 15) + th

del_M_plot <- ggplot()+
  geom_point(data = subset(dat_all_e, BirdID != "21_0203" &species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade)) +
  geom_point(data = subset(dat_all_t, BirdID != "21_0203" &species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade), pch = 15) +
  th +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

cg_plot <- ggplot()+
  geom_point(data = dat_final, aes(x = full_m, y = full_CGz_specific, col = clade)) + th

CG_range_plot <- ggplot()+
  geom_point(data = dat_bird, aes(x = total_bird_mass, y = range_CGx_specific, col = clade)) + th +
  scale_x_continuous(trans='log10')


