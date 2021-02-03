# This script was written to analyse the data that is returned from the birdmoment package

# --------------------- Read in data -----------------------
# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

filename_wing = list.files(path = path_data_folder, pattern = paste("allspecimen_winginfo"))
dat_wing      = read.csv(file = paste(path_data_folder,filename_wing,sep= ""))
filename_feat = list.files(path = path_data_folder, pattern = paste("allspecimen_featherinfo"))
dat_feat      = read.csv(file = paste(path_data_folder,filename_feat,sep= ""))
dat_feat      = dat_feat[,-c(1,3:8)]
test          = reshape2::dcast(dat_feat, species + BirdID ~ feather, value.var=c("l_vane","l_cal","w_cal","w_vd","w_vp","m_f"))

filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))
filename_body = list.files(path = path_data_folder, pattern = paste("bodyCGandMOI"))
dat_body      = read.csv(file = paste(path_data_folder,filename_body,sep= ""))
dat_body <- reshape2::dcast(dat_body, species + BirdID + TestID + FrameID ~ component + object, value.var="value")
remove(filename_wing,filename_feat,filename_bird,filename_body)

filename_results = list.files(path = path_data_folder, pattern = paste("results"))
dat_results = read.csv(paste(path_data_folder,filename_results[1],sep=""))
for (i in 2:length(filename_results)){
  dat_results = rbind(dat_results,read.csv(paste(path_data_folder,filename_results[i],sep="")))
}

# Merge with all info we have about the individual
dat_results <- merge(dat_results,dat_wing[,c("species","BirdID","TestID","FrameID","elbow","manus","pt1_X","pt1_Z")], by = c("species","BirdID","TestID","FrameID"))
dat_results <- merge(dat_results,dat_body[,-c(3,4)], by = c("species","BirdID"))
dat_results <- merge(dat_results,dat_bird[,-c(1,4)], by = c("species","BirdID"))
dat_results$torso_length <- dat_results$torsotail_length - dat_results$tail_length

# ------------- Check how the elbow and wrist affects the mass properties of the wing -----------------
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
    plot.background  = ggplot2::element_rect(fill = "transparent", color = NA)
  )

ROM_plot <- ggplot()+
  geom_point(data = dat_results, aes(x = elbow, y = manus, col = species)) + th

cg_plot <- ggplot()+
  geom_point(data = dat_results, aes(x = full_m, y = full_CGx/torsotail_length, col = species)) + th

I_plot <- ggplot()+
  geom_point(data = dat_results, aes(x = elbow, y = wing_CGy/wing_span_cm, col = species)) + th

I_specific_plot <- ggplot()+
  geom_point(data = dat_results, aes(x = manus, y = full_Ixz/(full_m*torsotail_length^2), col = species)) + th

I_specific_plot <- ggplot()+
  geom_point(data = dat_results, aes(x = full_m, y = full_CGx, col = species)) + th

test <- lm(full_CGz/torsotail_length ~ elbow*manus+species, data = dat_results)
test <- lm(full_Ixz/(full_m*torsotail_length^2) ~ elbow*manus+species, data = dat_results)
