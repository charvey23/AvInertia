# This script was written to analyse the data that is returned from the birdmoment package
library(pracma)
library(ggplot2)

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

## TO CHECK DATA
# x_vertices = c(dat_wing$pt1_X[i], dat_wing$pt6_X[i],dat_wing$pt7_X[i],dat_wing$pt8_X[i],dat_wing$pt9_X[i],dat_wing$pt10_X[i],dat_wing$pt11_X[i],dat_wing$pt12_X[i])
# y_vertices = c(dat_wing$pt1_Y[i], dat_wing$pt6_Y[i],dat_wing$pt7_Y[i],dat_wing$pt8_Y[i],dat_wing$pt9_Y[i],dat_wing$pt10_Y[i],dat_wing$pt11_Y[i],dat_wing$pt12_Y[i])
# z_vertices = c(dat_wing$pt1_Z[i], dat_wing$pt6_Z[i],dat_wing$pt7_Z[i],dat_wing$pt8_Z[i],dat_wing$pt9_Z[i],dat_wing$pt10_Z[i],dat_wing$pt11_Z[i],dat_wing$pt12_Z[i])


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

# --------------------- Read in data -----------------------
# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

filename_wing = list.files(path = path_data_folder, pattern = paste("allspecimen_winginfo"))
dat_wing      = read.csv(file = paste(path_data_folder,filename_wing,sep= ""))

filename_feat = list.files(path = path_data_folder, pattern = paste("allspecimen_featherinfo"))
dat_feat      = read.csv(file = paste(path_data_folder,filename_feat,sep= ""))
dat_feat      = dat_feat[,-c(1,3:7)]

filename_bird = list.files(path = path_data_folder, pattern = paste("allspecimen_birdinfo"))
dat_bird      = read.csv(file = paste(path_data_folder,filename_bird,sep= ""))
dat_bird$clade <- c("Accipitriformes","Accipitriformes","Galloanserae","Galloanserae","Aequorlitornithes","Aequorlitornithes","Strisores","Strisores",
                    "Galloanserae","Columbaves","Columbaves","Columbaves","Passeriformes","Passeriformes","Passeriformes","Passeriformes","Strisores",
                    "Australaves","Australaves","Australaves","Australaves","Galloanserae","Galloanserae","Coraciimorphae","Coraciimorphae",
                    "Coraciimorphae","Aequorlitornithes","Aequorlitornithes","Strigiformes","Strigiformes")

filename_body = list.files(path = path_data_folder, pattern = paste("bodyCGandMOI"))
dat_body      = read.csv(file = paste(path_data_folder,filename_body,sep= ""))
dat_body <- reshape2::dcast(dat_body, species + BirdID + TestID + FrameID ~ component + object, value.var="value")
remove(filename_wing,filename_feat,filename_bird,filename_body)

filename_results = list.files(path = path_data_folder, pattern = paste("results"))
dat_results = read.csv(paste(path_data_folder,filename_results[1],sep=""))
for (i in 2:length(filename_results)){
  dat_results = rbind(dat_results,read.csv(paste(path_data_folder,filename_results[i],sep="")))
}
# adjust the sharp-shinned hawk 20_1016 to a coopers hawk
dat_wing$species[which(dat_wing$species == "acc_str" & dat_wing$BirdID == "20_1016")] = "acc_coo"
dat_body$species[which(dat_body$species == "acc_str" & dat_body$BirdID == "20_1016")] = "acc_coo"

# remove wings where they accidentally got flipped around 44 from chr_amh and 7 from oce_leu CHECK THIS as you add more
# dat_results[which(dat_wing$pt2_X > 0 & dat_wing$species == "chr_amh")]

# Merge with all info we have about the individual
dat_all <- merge(dat_results,dat_wing[,c("species","BirdID","TestID","FrameID","elbow","manus","pt1_X","pt1_Z",
                                         "pt6_X","pt6_Y","pt7_X","pt7_Y","pt8_X","pt8_Y","pt9_X","pt9_Y","pt10_X","pt10_Y",
                                         "pt11_X","pt11_Y","pt11_Z","pt12_X","pt12_Y","pt12_Z")], by = c("species","BirdID","TestID","FrameID"))
dat_all <- merge(dat_all,dat_body[,-c(3,4)], by = c("species","BirdID"))
dat_all <- merge(dat_all,dat_bird[,-c(1,4)], by = c("species","BirdID"))
dat_all$torso_length <- dat_all$torsotail_length - dat_all$tail_length

# Adjust the calculations accordingly
dat_all$full_VRP_Ixx <- dat_all$head_Ixx+dat_all$neck_Ixx+dat_all$tail_Ixx+dat_all$torso_Ixx+2*dat_all$wing_Ixx
dat_all$full_VRP_Iyy <- dat_all$head_Iyy+dat_all$neck_Iyy+dat_all$tail_Iyy+dat_all$torso_Iyy+2*dat_all$wing_Iyy
dat_all$full_VRP_Izz <- dat_all$head_Izz+dat_all$neck_Izz+dat_all$tail_Izz+dat_all$torso_Izz+2*dat_all$wing_Izz
dat_all$full_VRP_Ixz <- dat_all$head_Ixz+dat_all$neck_Ixz+dat_all$tail_Ixz+dat_all$torso_Ixz+2*dat_all$wing_Ixz

dat_all$full_CGx <- ((dat_all$head_CGx*dat_all$head_m)+(dat_all$neck_CGx*dat_all$neck_m)+(dat_all$tail_CGx*dat_all$tail_m)+
                     (dat_all$torso_CGx*dat_all$torso_m)+2*(dat_all$wing_CGx*dat_all$wing_m))/dat_all$full_m
dat_all$full_CGz <- ((dat_all$head_CGz*dat_all$head_m)+(dat_all$neck_CGz*dat_all$neck_m)+(dat_all$tail_CGz*dat_all$tail_m)+
                       (dat_all$torso_CGz*dat_all$torso_m)+2*(dat_all$wing_CGz*dat_all$wing_m))/dat_all$full_m

dat_all$S_proj <- 0
for(i in 1:nrow(dat_all)){
  I_vrp = matrix(0, nrow = 3, ncol = 3)
  I_vrp[1,1] = dat_all$full_VRP_Ixx[i]
  I_vrp[2,2] = dat_all$full_VRP_Iyy[i]
  I_vrp[3,3] = dat_all$full_VRP_Izz[i]
  I_vrp[3,1] = dat_all$full_VRP_Ixz[i]
  I_vrp[1,3] = dat_all$full_VRP_Ixz[i]

  CG = c(dat_all$full_CGx[i], 0, dat_all$full_CGz[i])
  I = parallelaxis(I_vrp,-CG,dat_all$full_m[i],"A")

  dat_all$full_Ixx[i] = I[1,1]
  dat_all$full_Iyy[i] = I[2,2]
  dat_all$full_Izz[i] = I[3,3]
  dat_all$full_Ixz[i] = I[3,1]
  # Calculate the projected area for each wing - this is the correct order because X is negative and Y is positive
  x_vertices = c(dat_all$pt6_X[i],dat_all$pt7_X[i],dat_all$pt8_X[i],dat_all$pt9_X[i],dat_all$pt10_X[i],dat_all$pt11_X[i],dat_all$pt12_X[i])
  y_vertices = c(dat_all$pt6_Y[i],dat_all$pt7_Y[i],dat_all$pt8_Y[i],dat_all$pt9_Y[i],dat_all$pt10_Y[i],dat_all$pt11_Y[i],dat_all$pt12_Y[i])
  dat_all$S_proj[i] <- polyarea(x_vertices, y_vertices)
}


dat_all$full_CGx_specific <- (dat_all$full_CGx-dat_all$pt1_X)/dat_all$torso_length
dat_all$full_CGz_specific <- (dat_all$full_CGz-dat_all$pt1_Z)/dat_all$torso_length

### ------------- Compute summed quantities -----------------

# Maximum projected wing area
test     <- aggregate(list(S_proj_max = dat_all$S_proj),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), max)
dat_bird <- merge(dat_bird[,-c(1)],test, by = c("species","BirdID"))
dat_all  <- merge(dat_all,test, by = c("species","BirdID"))
# Range of CGx specific
test     <- aggregate(list(range_CGx_specific = dat_all$full_CGx_specific),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), FUN=function(x){max(x)-min(x)})
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
# Maximum of CGx specific
test     <- aggregate(list(max_CGx_specific = dat_all$full_CGx_specific),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
# Minimum of CGx specific
test     <- aggregate(list(min_CGx_specific = dat_all$full_CGx_specific),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), min)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
# Range of CGz specific
test     <- aggregate(list(range_CGz_specific = dat_all$full_CGz_specific),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), FUN=function(x){max(x)-min(x)})
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
# Maximum of CGz specific
test     <- aggregate(list(max_CGz_specific = dat_all$full_CGz_specific),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
# Minimum of CGz specific
test     <- aggregate(list(min_CGz_specific = dat_all$full_CGz_specific),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), min)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))

# Maximum of wing Ixx
test     <- aggregate(list(max_wing_Ixx = dat_all$wing_Ixx),  by=list(species = dat_all$species, BirdID = dat_all$BirdID), max)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
# Mean "root chord"
test     <- aggregate(list(mean_rc = sqrt((dat_all$pt11_X-dat_all$pt12_X)^2+(dat_all$pt11_Y-dat_all$pt12_Y)^2+(dat_all$pt11_Z-dat_all$pt12_Z)^2)),
                  by=list(species = dat_all$species, BirdID = dat_all$BirdID), mean)
dat_bird <- merge(dat_bird,test, by = c("species","BirdID"))
dat_all  <- merge(dat_all, test, by = c("species","BirdID"))
### ---------------------- Compute inertial metrics -------------------------

# uses scaling from:
# Alerstam, T., Rosén, M., Bäckman, J., Ericson, P. G., & Hellgren, O. (2007).
# Flight speeds among bird species: allometric and phylogenetic effects. PLoS Biol, 5(8), e197.
dat_all$prop_q_dot     <- (-(dat_all$full_CGx-dat_all$pt1_X)*dat_all$S_proj_max*dat_all$full_m^0.24)/(dat_all$full_Iyy)
dat_all$prop_q_dot_nd  <- (-(dat_all$full_CGx-dat_all$pt1_X)*dat_all$S_proj_max*dat_all$body_length^2)/(dat_all$full_Iyy)
dat_all$del_M_specific <- dat_all$prop_q_dot*dat_all$full_Iyy/(dat_all$full_m*dat_all$body_length)


filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_alldata.csv",sep="")
write.csv(dat_all,filename)

CG_range_plot <- ggplot()+
  geom_point(data = dat_bird, aes(x = clade, y = max_CGx_specific, col = species), pch = 17) +
  geom_point(data = dat_bird, aes(x = clade, y = min_CGx_specific, col = species), pch = 15) + th

ROM_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = elbow, y = manus, col = abs(full_Ixz/full_Izz))) + th + facet_wrap(~species)

q_dot_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = full_m, y = prop_q_dot, col = species)) + th +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

del_M_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = full_m, y = del_M_specific, col = clade)) + th +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

cg_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = torso_length, y = mean_rc, col = species)) + th  +
  scale_x_continuous(trans='log10')

CG_range_plot <- ggplot()+
  geom_point(data = dat_bird, aes(x = total_bird_mass, y = range_CGx_specific, col = clade)) + th +
  scale_x_continuous(trans='log10')

Ixz_plot <- ggplot()+
  geom_point(data = subset(dat_all, species != "oce_leu"), aes(x = full_m, y = wing_Ixz/(full_m*torso_length^2), col = species)) + th

test <- lmer(prop_q_dot ~ elbow*manus+species+full_m+(1|BirdID), data = dat_all)
test <- lm((full_Ixz/full_Izz) ~ elbow*manus+species+full_m, data = dat_all)

#Compare to Kirkpatrick
plot(as.factor(dat_results$species), dat_results$full_Ixx)
points(as.factor(dat_results$species),3.76*10^-3*dat_results$full_m^2.05, col = "red")


## ----------------- Plots for sharing with others -----------------

validation_Ixx_plot <- ggplot()+
  geom_point(data = dat_all, aes(x = full_m, y = wing_Ixx, col = clade)) +
  geom_line(data = unique(dat_all[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_all[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_all[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Izz_plot <- ggplot()+
  geom_point(data = dat_all, aes(x = full_m, y = wing_Izz, col = clade)) +
  geom_line(data = unique(dat_all[,c("species","BirdID","full_m")]), aes(x = full_m, y = 3.76*10^-3*full_m^2.05),col = "black") + # Kirkpatrick, 1994
  geom_line(data = unique(dat_all[,c("species","BirdID","full_m")]), aes(x = full_m, y = 1.94*10^-3*full_m^1.953),col = "gray") + # Berg and Rayner, 1995
  geom_point(data = unique(dat_all[,c("species","BirdID","full_m","wing_m","wing_span_cm")]), aes(x = full_m, y = full_m*(sqrt((0.14*wing_m/full_m))*wing_span_cm*0.5)^2),col = "black",pch = "-",size = 7) + # Sachs, 2005
  th  +
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')

validation_Ixx_body_plot <- ggplot()+
  geom_point(data = dat_all, aes(x = full_m, y = full_Ixx, col = clade)) +
  geom_point(aes(x = 0.0875, y = 928/(100^2*1000)), col = "black") + # Hedrick and Biewener
  geom_point(aes(x = 0.2893, y = 12889/(100^2*1000)), col = "black") + # Hedrick and Biewener
  th +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

