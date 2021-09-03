# This script is used to verify which estimation of the chord makes the best fit to the gull neutral point motion
library(pracma)
library(ggplot2)
library(ggthemes)  # need for geom_rangeframe
# This function must be run after code from our previous study
# https://github.com/charvey23/UMWTStaticStability - analyse_lltoutputs.R
source("/Users/christinaharvey/Documents/AvInertia/AnalysisFunctions/calc_np_loc_functions.R")
dat_stab$x_np = dat_stab$cmcl*dat_stab$c_max

# Calculate the project area for each wing
dat_all$S_wing  <- 0
for (i in 1:length(dat_all$frameID)){
  ## ----- Calculate the total wing area ------
  A1 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt11X","Pt11Y","Pt11Z")]-dat_all[i,c("Pt2X","Pt2Y","Pt2Z")])),
                          as.vector(t(dat_all[i,c("Pt11X","Pt11Y","Pt11Z")]-dat_all[i,c("Pt12X","Pt12Y","Pt12Z")]))), p = 2))
  A2 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt2X","Pt2Y","Pt2Z")]-dat_all[i,c("Pt6X","Pt6Y","Pt6Z")])),
                          as.vector(t(dat_all[i,c("Pt2X","Pt2Y","Pt2Z")]-dat_all[i,c("Pt12X","Pt12Y","Pt12Z")]))), p = 2))
  A3 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt11X","Pt11Y","Pt11Z")]-dat_all[i,c("Pt10X","Pt10Y","Pt10Z")])),
                          as.vector(t(dat_all[i,c("Pt11X","Pt11Y","Pt11Z")]-dat_all[i,c("Pt2X","Pt2Y","Pt2Z")]))), p = 2))
  A4 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt10X","Pt10Y","Pt10Z")]-dat_all[i,c("Pt6X","Pt6Y","Pt6Z")])),
                          as.vector(t(dat_all[i,c("Pt10X","Pt10Y","Pt10Z")]-dat_all[i,c("Pt2X","Pt2Y","Pt2Z")]))), p = 2))
  A5 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt10X","Pt10Y","Pt10Z")]-dat_all[i,c("Pt9X","Pt9Y","Pt9Z")])),
                          as.vector(t(dat_all[i,c("Pt10X","Pt10Y","Pt10Z")]-dat_all[i,c("Pt6X","Pt6Y","Pt6Z")]))), p = 2))
  A6 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt9X","Pt9Y","Pt9Z")]-dat_all[i,c("Pt7X","Pt7Y","Pt7Z")])),
                          as.vector(t(dat_all[i,c("Pt9X","Pt9Y","Pt9Z")]-dat_all[i,c("Pt6X","Pt6Y","Pt6Z")]))), p = 2))
  A7 = abs(0.5*Norm(cross(as.vector(t(dat_all[i,c("Pt9X","Pt9Y","Pt9Z")]-dat_all[i,c("Pt8X","Pt8Y","Pt8Z")])),
                          as.vector(t(dat_all[i,c("Pt9X","Pt9Y","Pt9Z")]-dat_all[i,c("Pt7X","Pt7Y","Pt7Z")]))), p = 2))

  dat_all$S_wing[i] = sum(A1,A2,A3,A4,A5,A6,A7)
}

dat_chord <- merge(dat_stab,dat_all, by =c("FrameID","WingID","TestID","elbow","manus"))
dat_chord1 <- aggregate(list(b = dat_num$b,
                        MAC = dat_num$MAC,
                        S_MX = dat_num$S),  by=list(species = dat_num$species,
                                                 WingID = dat_num$WingID,
                                                 TestID = dat_num$TestID,
                                                 FrameID = dat_num$FrameID), max)
dat_chord <- merge(dat_chord,dat_chord1, by =c("FrameID","WingID","TestID"), all = FALSE, all.y = FALSE)
colnames(dat_chord) <- make.unique(names(dat_chord))

dat_chord$c_root_max <- 0
dat_chord$c_root_max[which(dat_chord$WingID == "17_0285")] = max(dat_chord$c_root[which(dat_chord$WingID == "17_0285")])
dat_chord$c_root_max[which(dat_chord$WingID == "17_0243")] = max(dat_chord$c_root[which(dat_chord$WingID == "17_0243")])
dat_chord$c_root_max[which(dat_chord$WingID == "16_0048")] = max(dat_chord$c_root[which(dat_chord$WingID == "16_0048")])

dat_chord$mean_proj_chord_nobody = (dat_chord$S_proj/(dat_chord$b-(0.5*0.09481795)))
dat_chord$mean_chord_nobody = (dat_chord$S_wing/(dat_chord$b-(0.5*0.09481795)))

dat_chord$c_mean_max <- 0
dat_chord$c_mean_max[which(dat_chord$WingID == "17_0285")] = max(dat_chord$mean_chord_nobody[which(dat_chord$WingID == "17_0285")])
dat_chord$c_mean_max[which(dat_chord$WingID == "17_0243")] = max(dat_chord$mean_chord_nobody[which(dat_chord$WingID == "17_0243")])
dat_chord$c_mean_max[which(dat_chord$WingID == "16_0048")] = max(dat_chord$mean_chord_nobody[which(dat_chord$WingID == "16_0048")])

### ------------- Calculate the quarter chord position for key chord metrics ---------

dat_chord$mean_proj_ac_nobody <- 0
dat_chord$mean_ac_nobody      <- 0
dat_chord$mac_ac_nobody       <- 0
dat_chord$cen_ac_nobody       <- 0
dat_chord$st_mean_ac_nobody   <- 0
for (k in 1:nrow(dat_chord)){
  peri_pts = rbind(c(dat_chord$Pt12X[k],dat_chord$Pt12Y[k],dat_chord$Pt12Z[k]),
                   c(dat_chord$Pt6X[k],dat_chord$Pt6Y[k],dat_chord$Pt6Z[k]),
                   c(dat_chord$Pt7X[k],dat_chord$Pt7Y[k],dat_chord$Pt7Z[k]),
                   c(dat_chord$Pt8X[k],dat_chord$Pt8Y[k],dat_chord$Pt8Z[k]),
                   c(dat_chord$Pt9X[k],dat_chord$Pt9Y[k],dat_chord$Pt9Z[k]),
                   c(dat_chord$Pt10X[k],dat_chord$Pt10Y[k],dat_chord$Pt10Z[k]),
                   c(dat_chord$Pt11X[k],dat_chord$Pt11Y[k],dat_chord$Pt11Z[k]))
  if (dat_chord$Pt11Y[k] > dat_chord$Pt12Y[k]){
    peri_pts = rbind(peri_pts,c(dat_chord$Pt11X[k],dat_chord$Pt12Y[k],dat_chord$Pt12Z[k]))
  }
  dat_chord$mean_proj_ac_nobody[k] <- find_c4_x(peri_pts,dat_chord$mean_proj_chord_nobody[k])
  dat_chord$mean_ac_nobody[k] <- find_c4_x(peri_pts,dat_chord$mean_chord_nobody[k])
  dat_chord$mac_ac_nobody[k] <- find_c4_x(peri_pts,"MAC")
  dat_chord$cen_ac_nobody[k] <- find_c4_x(peri_pts,"area centroid")
  dat_chord$st_mean_ac_nobody[k] <- find_c4_x(peri_pts,"standard mean")
}
# Metric 1
dat_chord$root_ac_nondim      <- dat_chord$c_root/dat_chord$c_root_max
# Metric 2
dat_chord$mean_proj_ac_nondim <- dat_chord$mean_proj_ac_nobody/dat_chord$c_root_max
# Metric 3
dat_chord$mean_ac_nondim      <- dat_chord$mean_ac_nobody/dat_chord$c_root_max
# Metric 4
dat_chord$mac_ac_nondim       <- dat_chord$mac_ac_nobody/dat_chord$c_root_max
# Metric 5
dat_chord$cen_ac_nondim       <- dat_chord$cen_ac_nobody/dat_chord$c_root_max
# Metric 6
dat_chord$st_mean_ac_nondim   <- dat_chord$st_mean_ac_nobody/dat_chord$c_root_max

# Dependent value - true neutral point
dat_chord$x_np_nondim         <- dat_chord$x_np/dat_chord$c_root_max

# Evaluate metric 1
test_1 <- lm(log(-x_np_nondim)~log(root_ac_nondim), data = dat_chord)
summary(test_1)
confint(test_1)
# Evaluate metric 2
test_2 <- lm(log(-x_np_nondim)~log(-mean_proj_ac_nondim), data = dat_chord)
summary(test_2)
confint(test_2)
# Evaluate metric 3
test_3 <- lm(log(-x_np_nondim)~log(-mean_ac_nondim), data = dat_chord)
summary(test_3)
confint(test_3)
# Evaluate metric 4
test_4 <- lm(log(-x_np_nondim)~log(-mac_ac_nondim), data = dat_chord)
summary(test_4)
confint(test_4)
# Evaluate metric 5
test_5 <- lm(log(-x_np_nondim)~log(-cen_ac_nondim), data = dat_chord)
summary(test_5)
confint(test_5)
# Evaluate metric 6
test_6 <- lm(log(-x_np_nondim)~log(-st_mean_ac_nondim), data = dat_chord)
summary(test_6)
confint(test_6)
#
# ### Option #1: Use the root chord change to estimate the neutral point position
# dat_chord$c_root_max <- 0
# dat_chord$c_root_max[which(dat_chord$WingID == "17_0285")] = max(dat_chord$c_root[which(dat_chord$WingID == "17_0285")])
# dat_chord$c_root_max[which(dat_chord$WingID == "17_0243")] = max(dat_chord$c_root[which(dat_chord$WingID == "17_0243")])
# dat_chord$c_root_max[which(dat_chord$WingID == "16_0048")] = max(dat_chord$c_root[which(dat_chord$WingID == "16_0048")])
#
# dat_chord$x_root_ac = -0.25*dat_chord$c_root
# # center at the mean
# dat_chord$x_root_ac = dat_chord$x_root_ac-mean(dat_chord$x_root_ac)
# mod_root <- lm(x_np ~ x_root_ac, data = subset(dat_chord,WingID=="17_0285"))
# ### Option #2: Use the mean chord change to estimate the neutral point position
# # Option #2a: the project wing area divided by the half span
# #this is the one wing projected area divided by the half span which includes the body width
# # substract half the constant body width from this one
# # Option #2a - no body area
# dat_chord$x_mean_proj_ac1 = -0.25*(dat_chord$S_proj/(dat_chord$b-(0.5*0.09481795)))
# # center at the mean
# dat_chord$x_mean_proj_ac1 = dat_chord$x_mean_proj_ac1-mean(dat_chord$x_mean_proj_ac1)
# #fit model
# mod_mean_proj1 <- lm(x_np ~ x_mean_proj_ac1, data = dat_chord)
# # Option #2b -include the body area
# dat_chord$x_mean_proj_ac2 = -0.25*((dat_chord$S_proj+(0.5*0.00305))/(dat_chord$b))
# # center at the mean
# dat_chord$x_mean_proj_ac2 = dat_chord$x_mean_proj_ac2-mean(dat_chord$x_mean_proj_ac2)
# #fit model
# mod_mean_proj2 <- lm(x_np ~ x_mean_proj_ac2, data = dat_chord)
#
# # Option #2c - no body area
# dat_chord$x_mean_ac1 = -0.25*(dat_chord$S_wing/(dat_chord$b-(0.5*0.09481795)))
# # center at the mean
# dat_chord$x_mean_ac1 = dat_chord$x_mean_ac1-mean(dat_chord$x_mean_ac1)
# #fit model
# mod_mean1 <- lm(x_np ~ x_mean_ac1, data = dat_chord)
# # Option #2d - include the body area
# dat_chord$x_mean_ac2 = -0.25*((dat_chord$S_wing+(0.5*0.00305))/(dat_chord$b))
# # center at the mean
# dat_chord$x_mean_ac2 = dat_chord$x_mean_ac2-mean(dat_chord$x_mean_ac2)
# #fit model
# mod_mean2 <- lm(x_np ~ x_mean_ac2, data = dat_chord)
#
# ### Option #3: Use the mean aerodynamic chord change to estimate the neutral point position
# dat_chord$x_mac_ac = -0.25*dat_chord$MAC
# # center at the mean
# dat_chord$x_mac_ac = dat_chord$x_mac_ac-mean(dat_chord$x_mac_ac)
# #fit model
# mod_MAC <- lm(x_np ~ x_mac_ac, data = dat_chord)
#
# th <- theme_classic() +
#   theme(
#     # Text
#     axis.title = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5),
#     axis.text = element_text(size = 8, colour = "black"),
#     axis.text.x = element_text(margin = margin(t = 10, unit = "pt")),
#     axis.text.y = element_text(margin = margin(r = 10)),
#     # Axis line
#     axis.line = element_blank(),
#     axis.ticks.length = unit(-5,"pt"),
#     # Legend
#     legend.position = 'none',
#     # Background transparency
#     # Background of panel
#     panel.background = element_rect(fill = "transparent"),
#     # Background behind actual data points
#     plot.background = element_rect(fill = "transparent", color = NA)
#   )
#
# plot_chord <- ggplot()+
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   geom_smooth(data = dat_chord, aes(x = x_mean_proj_ac1, y = x_np), col = "red", linetype = 3, method = "lm", formula = "y~x", size = 0.5) +
#   geom_smooth(data = dat_chord, aes(x = x_mean_proj_ac2, y = x_np), col = "red",  method = "lm", formula = "y~x",size = 0.5) +
#   geom_smooth(data = dat_chord, aes(x = x_mean_ac1, y = x_np), method = "lm", linetype = 3, formula = "y~x",size = 0.5) +
#   geom_smooth(data = dat_chord, aes(x = x_mean_ac2, y = x_np), method = "lm", formula = "y~x",size = 0.5) +
#   geom_smooth(data = dat_chord, aes(x = x_root_ac, y = x_np), method = "lm", formula = "y~x",size = 0.5) +
#   geom_smooth(data = dat_chord, aes(x = x_mac_ac, y = x_np), col = "green", fill= "green", alpha = 0.3, method = "lm", formula = "y~x",size = 0.5) + th +
#   coord_fixed()+
#   scale_y_continuous(limits = c(-0.03,0.03), breaks = c(-0.03,-0.015,0,0.015,0.03), name = "Neutral point - centered (m)") +
#   scale_x_continuous(limits = c(-0.03,0.03), breaks = c(-0.03,-0.015,0,0.015,0.03), name = "Proxy chord length - centered (m)") +
#   geom_rangeframe() +
#   annotate(geom = "segment", x = log(0), xend = log(0), y = -0.03, yend = 0.03) +
#   annotate(geom = "segment", x = -0.03, xend = 0.03, y = log(0), yend = log(0))




