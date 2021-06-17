#### ----- Script containing all aesthetics needed for plotting ----------


# ---------- Pre-define themes and base plots ----------
# base theme
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

# blank plot to fill space
blank_plot <- ggplot() + theme_void()


# ---------- Pre-define colours ----------
col_x = "#09015f"
col_y = "#55b3b1"
col_z = "#f6c065"
col_xz = "#af0069"

col_elb = "#55b3b1"
col_man = "#af0069"
col_elbman = "#5F0995"

# Species rainbow colours
cc_rain1  <- scales::seq_gradient_pal("#DA9101", "#CA302F", "Lab")(seq(0,1,length.out=6))
cc_rain2  <- scales::seq_gradient_pal("#FCC201", "#DA9101", "Lab")(seq(0,1,length.out=3))
cc_rain3  <- scales::seq_gradient_pal("#138715", "#FCC201", "Lab")(seq(0,1,length.out=6))
cc_rain4  <- scales::seq_gradient_pal("#1FC3CD", "#138715", "Lab")(seq(0,1,length.out=4))
cc_rain5  <- scales::seq_gradient_pal("#304CC8", "#1FC3CD", "Lab")(seq(0,1,length.out=5))
cc_rain6  <- scales::seq_gradient_pal("#5F2CC8", "#304CC8", "Lab")(seq(0,1,length.out=3))
cc_rain   <- c(cc_rain6, cc_rain5[2:5], cc_rain4[2:4], cc_rain3[2:6],cc_rain2[2:3],cc_rain1[2:6])

# Variation colours
#cc_var8  <- scales::seq_gradient_pal("#004237", "black", "Lab")(seq(0,1,length.out=8))
cc_var0  <- scales::seq_gradient_pal("#207567", "#004237", "Lab")(seq(0,1,length.out=10))
cc_var1  <- scales::seq_gradient_pal("#358873", "#207567", "Lab")(seq(0,1,length.out=5))
cc_var2  <- scales::seq_gradient_pal("#4E9C81", "#358873", "Lab")(seq(0,1,length.out=3))
cc_var3  <- scales::seq_gradient_pal("#6BAF92", "#4E9C81", "Lab")(seq(0,1,length.out=3))
cc_var4  <- scales::seq_gradient_pal("#8DC3A7", "#6BAF92", "Lab")(seq(0,1,length.out=4))
cc_var5  <- scales::seq_gradient_pal("#B4D6C1", "#8DC3A7", "Lab")(seq(0,1,length.out=6))
cc_var6  <- scales::seq_gradient_pal("#DFEAE2", "#B4D6C1", "Lab")(seq(0,1,length.out=10))
cc_var7  <- scales::seq_gradient_pal("white", "#DFEAE2", "Lab")(seq(0,1,length.out=10))
cc_var   <- c(cc_var7[5:10], cc_var6[2:10], cc_var5[2:6], cc_var4[2:4], cc_var3[2:3],cc_var2[2:3],cc_var1[2:5],cc_var0[2:10])

# -------------- Set up common labels -----------
lab_eta = expression(paste(eta^2))

# --------------- Set up Phylogeny plot -----------
phylo_label = c("Cooper's hawk","Sharp-shinned hawk", "Western grebe", "Mallard","Great blue heron","Canada goose","Common nighthawk",
                "Lady Amherst's pheasant","Northern flicker","Pigeon","Common raven","Steller's jay","Black swift",
                "Merlin","Peregrine falcon","Leach's storm petrel","Glaucous-winged gull","Himalayan monal","Silver pheasant",
                "Belted kingfisher","American white pelican","Barn owl")
dat_phylo_plot <- data.frame(label = pruned_mcc$tip.label, genus = phylo_label)
phylo_plot <- ggtree(pruned_mcc) %<+% dat_phylo_plot %>% rotate(23) +
  geom_tiplab(aes(label=genus)) +
  theme(plot.margin = unit(c(14,8,14,8), "mm")) +
  theme_tree2() + geom_rootedge(3)
phylo_plot_complete <- revts(phylo_plot) +
  scale_x_continuous(breaks=c(-80,-60,-40,-20,0), labels=c(80,60,40,20,0))

phylo_data <- ggtree::ggtree(pruned_mcc)[["data"]][,c("label","y")]
phylo_order <- merge(phylo_data[complete.cases(phylo_data),], unique(dat_bird[,c("binomial","species")]), by.x = "label", by.y = "binomial")
phylo_order <- arrange(phylo_order,y)
# add new column that is factored in the correct order
dat_comp$species_order        = factor(dat_comp$species, levels = phylo_order$species)
dat_final$species_order       = factor(dat_final$species, levels = phylo_order$species)
shoulder_motion$species_order = factor(shoulder_motion$species, levels = phylo_order$species)
# create the shading for everyother species
dat_model_out$species <- factor(dat_model_out$species, levels = phylo_order$species)
shading <- data.frame(col1 = levels(dat_model_out$species)[seq(from = 1, to = max(as.numeric(as.factor(dat_model_out$species)))-1, by = 2)],
                      col2 = levels(dat_model_out$species)[seq(from = 2, to = max(as.numeric(as.factor(dat_model_out$species))), by = 2)])
