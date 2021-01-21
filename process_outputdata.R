# This script was written to analyse the data that is returned from the birdmoment package

# --------------------- Read in data -----------------------
# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

filename_wing = list.files(path = path_data_folder, pattern = paste("allspecimen_winginfo"))
dat_wing      = read.csv(file = paste(path_data_folder,filename_wing,sep= ""))
filename_feat = list.files(path = path_data_folder, pattern = paste("allspecimen_featherinfo"))
dat_feat      = read.csv(file = paste(path_data_folder,filename_feat,sep= ""))
filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))
remove(filename_wing,filename_feat,filename_bird)

filename_results = list.files(path = path_data_folder, pattern = paste("results"))
dat_results = read.csv(paste(path_data_folder,filename_results[1],sep=""))
for (i in 2:length(filename_results)){
  dat_results = rbind(dat_results,read.csv(paste(path_data_folder,filename_results[i],sep="")))
}

dat_results <- reshape2::dcast(dat_results, species + BirdID + TestID + FrameID ~ component + object, value.var="value") # Delete this eventually - will save the data differently
dat_results <- merge(dat_results,dat_wing[,c("species","BirdID","testid","frameID","elbow","manus")], by.x = c("species","BirdID","TestID","FrameID"), by.y = c("species","BirdID","testid","frameID"))

dat_results <- subset(dat_results, species != "meg_alc")

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
    # Legend
    legend.position  = 'none',
    # Background transparency
    # Background of panel
    panel.background = ggplot2::element_rect(fill = "transparent"),
    # Background behind actual data points
    plot.background  = ggplot2::element_rect(fill = "transparent", color = NA)
  )

cg_plot <- ggplot()+
  geom_point(data = dat_results, aes(x = bones_CGy, y = feathers_CGy, col = wing_CGy))
