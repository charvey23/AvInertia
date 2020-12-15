# ---- This code is to be used with the LLT takes raw digitized points and prepares them for testing
#When re-running for new data double check that the transformations are working correctly.

# CAUTION: This code is now set up to orient wings as follows:
#          1. Pt3 (Wrist) will be inline with Pt1 (humeral head) along the y axis
#          2. Pt3 (Wrist) will be inline with Pt1 (humeral head) along the z axis
#          3. Pt3 (Wrist) will be inline with Pt10 (S1) along the x axis

# - NOTE: For eagles we selected the point 8 to be longest primary P2 OR P3 and then point 9 to be P5
#-        Read in info and set working directory

# Created: Christina Harvey
# Last updated: 25-May-2020

# Load required libraries
library(reshape2)
library(pracma)   #needed for the math in 3D angles

#setwd("/Users/christinaharvey/Google Drive/ComparativeAvianStabilityStudy/CompStability(LLT&Digitize)/Code") #For Mac
setwd("/Users/christinaharvey/Google Drive/DoctoralThesis/WingMorphology/") #For Windows
source('jointangles.R')


# ------------------------- Set file directory -------------------------
setwd("/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/05_subsampled_optitrack")

# ------------------------- Read in data -------------------------
dat_info <- read.csv(file = "2020_12_15_IDfile.csv")

for (i in 5:nrow(dat_info)){
  # Read in specific wing data
  filename <- paste("2020_12_15_",dat_info$species[i],"_",dat_info$birdid[i],"_00",dat_info$testid[i],"_subsampled.csv",sep = "")
  rawdat   <- read.csv(filename, stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA") )

  if (i < 12){
    rawdat[,3:35] <- rawdat[,3:35]/1000  #First few are in mm
    # check humerus length
    mean(sqrt((rawdat$pt3_X - rawdat$pt2_X^2)^2 + (rawdat$pt3_Y - rawdat$pt2_Y^2)^2 + (rawdat$pt3_Z - rawdat$pt2_Z^2)^2))
  }

}


iterationcount <- nrow(rawdat)

# Visualize the wings as required - For each calibration verify that the axis is RH
plot(rawdat$pt2_X,rawdat$pt2_Y, xlim =c(-0.4,0.4), ylim = c(-0.4,0.4))
points(rawdat$pt3_X,rawdat$pt3_Y,xlim =c(-0.4,0.4), ylim = c(-0.4,0.4),col = "blue")
points(rawdat$pt4_X,rawdat$pt4_Y, xlim =c(-0.4,0.4), ylim = c(-0.4,0.4), col = "green")
points(rawdat$pt1_X,rawdat$pt1_Y, xlim =c(-0.4,0.4), ylim = c(-0.4,0.4), col = "red")


# ------------------------- Initialize all matrices -------------------------
clean        <- data.frame(matrix(nrow = nrow(rawdat)*49, ncol = 38))
names(clean) <- c('frameID','elbow','manus','sweep','dihedral',
                  'Pt1X','Pt1Y','Pt1Z','Pt2X','Pt2Y','Pt2Z','Pt3X','Pt3Y','Pt3Z',
                  'Pt4X','Pt4Y','Pt4Z','Pt6X','Pt6Y','Pt6Z',
                  'Pt7X','Pt7Y','Pt7Z','Pt8X','Pt8Y','Pt8Z','Pt9X','Pt9Y','Pt9Z',
                  'Pt10X','Pt10Y','Pt10Z','Pt11X','Pt11Y','Pt11Z','Pt12X','Pt12Y','Pt12Z')

pts    <- data.frame(matrix(nrow = nrow(rawdat) * 11, ncol = 5))
names(pts) <- c("ID", "Pt.No",
                "X", "Y", "Z")
pts.sh    <- data.frame(matrix(nrow = nrow(pts) * 49, ncol = 2))
names(pts.sh) <- c("ID", "Pt.No")
xaxis <- c(1, 0, 0)
yaxis <- c(0, 1, 0)
zaxis <- c(0, 0, 1)
# ------------------------- make raw data into long form -------------------------
count  = 1
for (j in 1:iterationcount) {
  iter   = 1
  for (k in 1:12) {
    if (k == 5){next} # no point 5 was saved
    pts[count, 1]   <- j
    pts[count, 2]   <- k
    pts[count, 3]   <- rawdat[j, iter]    #pt x
    pts[count, 4]   <- rawdat[j, 1 + iter]  #pt y
    pts[count, 5]   <- rawdat[j, 2 + iter]  #pt z
    count = count + 1
    iter  = iter  + 3
  }
}

#------------------------------- Step 1 -------------------------------
#Make pt 1 on humerus the beginning location
for (i in 1:iterationcount) {
  pts$pt1x[pts$ID == i] = pts$X[pts$ID == i & pts$Pt.No == 1]
  pts$pt1y[pts$ID == i] = pts$Y[pts$ID == i & pts$Pt.No == 1]
  pts$pt1z[pts$ID == i] = pts$Z[pts$ID == i & pts$Pt.No == 1]
}

pts$x_adj = pts$X - pts$pt1x
pts$y_adj = pts$Y - pts$pt1y
pts$z_adj = pts$Z - pts$pt1z

#-------------------------------  Step 2 -------------------------------
# Set pt 3 in line with Pt 1 along the wingspan
#Calculate angle between axis that runs along wing length and pt 3 (rotate about z)
#Project onto xaxis & make negative if j2 > 0: Gull 16-0048
for (i in 1:iterationcount) {
  i1 = pts$x_adj[pts$ID == unique(pts$ID)[i] & pts$Pt.No == 3]
  j1 = pts$y_adj[pts$ID == unique(pts$ID)[i] & pts$Pt.No == 3]
  k1 = 0
  vec      = c(i1, j1, k1)
  interior = dot(vec, xaxis) / (Norm(vec))

  pts$thetaz[pts$ID == unique(pts$ID)[i]]   = acos(interior) + pi / 2

  if (j1 > 0) {
    pts$thetaz[pts$ID == unique(pts$ID)[i]]   = -acos(interior) + pi / 2
  }

}
pts$x_adj1 = cos(pts$thetaz) * pts$x_adj - sin(pts$thetaz) * pts$y_adj
pts$y_adj1 = sin(pts$thetaz) * pts$x_adj + cos(pts$thetaz) * pts$y_adj
pts$z_adj1 = pts$z_adj

# CAUTION: Verify that Pt3X=x_adj1=0 after this step - View(subset(pts,Pt.No == 3))

#------------------------------- Step 3 -------------------------------
# Set pt 3 in line with Pt 1 looking down the length of wing
# Calculate angle between x axis and pt 3 (rotate about y)
# Project onto xaxis & make negative if k1 < 0: Gull 16-0048
for (i in 1:iterationcount) {
  i2 = 0
  j2 = pts$y_adj1[pts$ID == unique(pts$ID)[i] & pts$Pt.No == 3]
  k2 = pts$z_adj1[pts$ID == unique(pts$ID)[i] & pts$Pt.No == 3]

  vec      = c(i2, j2, k2)
  interior = dot(vec, yaxis) / (Norm(vec))

  pts$thetax[pts$ID == unique(pts$ID)[i]]   =  acos(interior)

  if (k2 > 0) {
    pts$thetax[pts$ID == unique(pts$ID)[i]]   = -acos(interior)
  }

}
pts$x_adj2 = pts$x_adj1
pts$y_adj2 = cos(pts$thetax) * pts$y_adj1 - sin(pts$thetax) * pts$z_adj1
pts$z_adj2 = sin(pts$thetax) * pts$y_adj1 + cos(pts$thetax) * pts$z_adj1


# CAUTION: Verify that Pt3Z=z_adj2=0 after this step - View(subset(pts,Pt.No == 3))

#------------------------------- Step 4 -------------------------------
# Set pt 10 in line with Pt 3 looking from back of wing
# Calculate angle between x axis and pt 10 (rotate about y)
# Project onto yaxis & make negative if k2 > 0: Gull 16-0048

for (i in 1:iterationcount) {
  i3 = pts$x_adj2[pts$ID == unique(pts$ID)[i] & pts$Pt.No == 10]
  j3 = 0
  k3 = pts$z_adj2[pts$ID == unique(pts$ID)[i] & pts$Pt.No == 10]

  vec      = c(i3, j3, k3)
  interior = dot(vec, xaxis) / (Norm(vec))

  pts$thetay[pts$ID == unique(pts$ID)[i]]   = -acos(interior)

  if (k3 > 0) {
    pts$thetay[pts$ID == unique(pts$ID)[i]]   = acos(interior)
  }

}

#- Step 6: Rotate about x by the calculated angle
  pts$x_adj3 = cos(pts$thetay)*pts$x_adj2+sin(pts$thetay)*pts$z_adj2
  pts$y_adj3 = pts$y_adj2
  pts$z_adj3 = -sin(pts$thetay)*pts$x_adj2+cos(pts$thetay)*pts$z_adj2
# CAUTION: Verify that Pt10Z=z_adj3=0 after this step - View(subset(pts,Pt.No == 10))


#------------------------------- Step 5 -------------------------------
# To be correct this will be a R wing
# y positive along wing, x positive towards nose, z positive down
# All   Gulls,  Hawk 18-1R & 18-2, Eagle 14-0324, Nighthawk 18-3 & 18-2:     xbody = -x; ybody = y; zbody = z
# Hawk 18-3, Eagle 17-113 & 17-193, Nighthawk 18-1L:                         xbody = -x; ybody = y; zbody = -z

pts$x_final = - pts$x_adj3
pts$y_final =   pts$y_adj3
pts$z_final = - pts$z_adj3


# --------------------- Shoulder sweep --------------------------------------------
# for each wing configuration we will sweep and dihedral about pt1
#------------------------------- Step 1 -------------------------------
sweep_inc    = c(-30,-20,-10,0,10,20,30)*pi/180
dihedral_inc = c(-30,-20,-10,0,10,20,30)*pi/180
# to sweep backwards rotate positive about z
pts.sh$ID[1:nrow(pts)]       = pts$ID
pts.sh$Pt.No[1:nrow(pts)]    = pts$Pt.No
pts.sh$sweep[1:nrow(pts)]    = 0
pts.sh$dihedral[1:nrow(pts)] = 0
pts.sh$x_sw[1:nrow(pts)]     = pts$x_final
pts.sh$y_sw[1:nrow(pts)]     = pts$y_final
pts.sh$z_sw[1:nrow(pts)]     = pts$z_final
pts.sh$x_swdi[1:nrow(pts)]   = pts$x_final
pts.sh$y_swdi[1:nrow(pts)]   = pts$y_final
pts.sh$z_swdi[1:nrow(pts)]   = pts$z_final

count = nrow(pts)+1
for (sw in 1:7){
  for (di in 1:7){
    if(sw == 4 & di == 4){next}
    pts.sh$ID[count:(count+nrow(pts)-1)]        = pts$ID
    pts.sh$Pt.No[count:(count+nrow(pts)-1)]     = pts$Pt.No
    pts.sh$sweep[count:(count+nrow(pts)-1)]     = sweep_inc[sw]*180/pi # to sweep backwards rotate positive about z
    pts.sh$x_sw[count:(count+nrow(pts)-1)]      = cos(sweep_inc[sw])*pts$x_final-sin(sweep_inc[sw])*pts$y_final
    pts.sh$y_sw[count:(count+nrow(pts)-1)]      = sin(sweep_inc[sw])*pts$x_final+cos(sweep_inc[sw])*pts$y_final
    pts.sh$z_sw[count:(count+nrow(pts)-1)]      = pts$z_final

      pts.sh$dihedral[count:(count+nrow(pts)-1)]  = -dihedral_inc[di]*180/pi  # positive rotation about x gives anhedral so flip sign s.t. positive dihedral indicates upwards

      pts.sh$x_swdi[count:(count+nrow(pts)-1)] = pts.sh$x_sw[count:(count+nrow(pts)-1)]
      pts.sh$y_swdi[count:(count+nrow(pts)-1)] = cos(dihedral_inc[di])*pts.sh$y_sw[count:(count+nrow(pts)-1)]-sin(dihedral_inc[di])*pts.sh$z_sw[count:(count+nrow(pts)-1)]
      pts.sh$z_swdi[count:(count+nrow(pts)-1)] = sin(dihedral_inc[di])*pts.sh$y_sw[count:(count+nrow(pts)-1)]+cos(dihedral_inc[di])*pts.sh$z_sw[count:(count+nrow(pts)-1)]
      count = count + nrow(pts)
  }
}

# --------------------------------- clean up the data ---------------------------------
clean$sweep    = pts.sh$sweep[pts.sh$Pt.No == 12]
clean$dihedral = pts.sh$dihedral[pts.sh$Pt.No == 12]

count2 = 1
for (ptno in 1:12){
  clean[,count2+5] = pts.sh$x_swdi[pts.sh$Pt.No == ptno]
  clean[,count2+6] = pts.sh$y_swdi[pts.sh$Pt.No == ptno]
  clean[,count2+7] = pts.sh$z_swdi[pts.sh$Pt.No == ptno]
  count2 = count2 + 3
}

#-------------------------- Calculate elbow and manus angles from the unrotated data --------------------------
elbow = 0
manus = 0
for (j in 1:nrow(rawdat)){
  elbow[j] <- jointangles(rawdat$pt1_X[j],rawdat$pt1_Y[j],rawdat$pt1_Z[j],
                          rawdat$pt2_X[j],rawdat$pt2_Y[j],rawdat$pt2_Z[j],
                          rawdat$pt3_X[j],rawdat$pt3_Y[j],rawdat$pt3_Z[j])

  manus[j] <- jointangles(rawdat$pt2_X[j],rawdat$pt2_Y[j],rawdat$pt2_Z[j],
                          rawdat$pt3_X[j],rawdat$pt3_Y[j],rawdat$pt3_Z[j],
                          rawdat$pt4_X[j],rawdat$pt4_Y[j],rawdat$pt4_Z[j])
}
#will autorepeat the vectors
clean$frameID  <- rawdat$frameID
clean$elbow    <- elbow
clean$manus    <- manus

# Remove unphysical configurations
removerow = which(clean$Pt2Y < dis_bodyedge | clean$Pt3Y < dis_bodyedge | clean$Pt4Y < dis_bodyedge | clean$Pt6Y < dis_bodyedge  | clean$Pt7Y < dis_bodyedge  | clean$Pt8Y < dis_bodyedge | clean$Pt9Y < dis_bodyedge | clean$Pt10Y < dis_bodyedge)
clean <- clean[-removerow,]

#For Windows
write.csv(clean,'/Users/Inman PC/Google Drive/DoctoralThesis/LLT_AIAAPaper/AvianWingLLT/ReorientedWings/2020_05_25_nighthawk_18_1L_6bar.csv')
#write.csv(clean,'/Users/christinaharvey/Google Drive/ComparativeAvianStabilityStudy/CompStability(LLT&Digitize)/LLTcode/AvianWingData/2020_02_21_gull_16_0048_Test6.csv')

m = 1:nrow(rawdat)
x = c(clean$Pt1X[m],clean$Pt2X[m],clean$Pt3X[m],clean$Pt4X[m],clean$Pt6X[m],clean$Pt7X[m],clean$Pt8X[m],clean$Pt9X[m],clean$Pt10X[m],clean$Pt11X[m],clean$Pt12X[m])
y = c(clean$Pt1Y[m],clean$Pt2Y[m],clean$Pt3Y[m],clean$Pt4Y[m],clean$Pt6Y[m],clean$Pt7Y[m],clean$Pt8Y[m],clean$Pt9Y[m],clean$Pt10Y[m],clean$Pt11Y[m],clean$Pt12Y[m])
z = c(clean$Pt1Z[m],clean$Pt2Z[m],clean$Pt3Z[m],clean$Pt4Z[m],clean$Pt6Z[m],clean$Pt7Z[m],clean$Pt8Z[m],clean$Pt9Z[m],clean$Pt10Z[m],clean$Pt11Z[m],clean$Pt12Z[m])

#Concatenate all final outputs
setwd("/Users/Inman PC/Google Drive/DoctoralThesis/LLT_AIAAPaper/AvianWingLLT/ReorientedWings") #For Windows

# ------------------------- Read in data -------------------------
rawdat1 <- read.csv('2020_05_25_eagle_14_0324_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat1$species = "hal_leu"
rawdat1$WingID = "14_0324"
rawdat1$TestID = "6bar"
rawdat2 <- read.csv('2020_05_25_eagle_14_0324_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat2$species = "hal_leu"
rawdat2$WingID = "14_0324"
rawdat2$TestID = "sweep"
rawdat3 <- read.csv('2020_05_25_eagle_17_113_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat3$species = "hal_leu"
rawdat3$WingID = "17_113"
rawdat3$TestID = "6bar"
rawdat4 <- read.csv('2020_05_25_eagle_17_113_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat4$species = "hal_leu"
rawdat4$WingID = "17_113"
rawdat4$TestID = "sweep"
rawdat5 <- read.csv('2020_05_25_eagle_17_193_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat5$species = "hal_leu"
rawdat5$WingID = "17_193"
rawdat5$TestID = "6bar"
rawdat6 <- read.csv('2020_05_25_eagle_17_193_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat6$species = "hal_leu"
rawdat6$WingID = "17_193"
rawdat6$TestID = "sweep"

rawdat7  <- read.csv('2020_05_25_hawk_18_1R_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat7$species = "but_jam"
rawdat7$WingID = "18_1R"
rawdat7$TestID = "6bar"
rawdat8  <- read.csv('2020_05_25_hawk_18_1R_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat8$species = "but_jam"
rawdat8$WingID = "18_1R"
rawdat8$TestID = "sweep"
rawdat9  <- read.csv('2020_05_25_hawk_18_2_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat9$species = "but_jam"
rawdat9$WingID = "18_2"
rawdat9$TestID = "6bar"
rawdat10 <- read.csv('2020_05_25_hawk_18_2_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat10$species = "but_jam"
rawdat10$WingID = "18_2"
rawdat10$TestID = "sweep"
rawdat11 <- read.csv('2020_05_25_hawk_18_3_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat11$species = "but_jam"
rawdat11$WingID = "18_3"
rawdat11$TestID = "6bar"
rawdat12 <- read.csv('2020_05_25_hawk_18_3_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat12$species = "but_jam"
rawdat12$WingID = "18_3"
rawdat12$TestID = "sweep"

rawdat13 <- read.csv('2020_05_25_nighthawk_18_1L_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat13$species = "cho_min"
rawdat13$WingID = "18_1L"
rawdat13$TestID = "6bar"
rawdat14 <- read.csv('2020_05_25_nighthawk_18_1L_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat14$species = "cho_min"
rawdat14$WingID = "18_1L"
rawdat14$TestID = "sweep"
rawdat15 <- read.csv('2020_05_25_nighthawk_18_2R_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat15$species = "cho_min"
rawdat15$WingID = "18_2R"
rawdat15$TestID = "6bar"
rawdat16 <- read.csv('2020_05_25_nighthawk_18_2R_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat16$species = "cho_min"
rawdat16$WingID = "18_2R"
rawdat16$TestID = "sweep"
rawdat17 <- read.csv('2020_05_25_nighthawk_18_3R_6bar.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat17$species = "cho_min"
rawdat17$WingID = "18_3R"
rawdat17$TestID = "6bar"
rawdat18 <- read.csv('2020_05_25_nighthawk_18_3R_sweep.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat18$species = "cho_min"
rawdat18$WingID = "18_3R"
rawdat18$TestID = "sweep"

rawdat19 <- read.csv('2020_05_25_gull_16_0048_Test5.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat19$species = "lar_gla"
rawdat19$WingID = "16_0048"
rawdat19$TestID = "Test5"
rawdat20 <- read.csv('2020_05_25_gull_16_0048_Test6.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat20$species = "lar_gla"
rawdat20$WingID = "16_0048"
rawdat20$TestID = "Test6"
rawdat21 <- read.csv('2020_05_25_gull_16_0048_Test7.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat21$species = "lar_gla"
rawdat21$WingID = "16_0048"
rawdat21$TestID = "Test7"
rawdat22 <- read.csv('2020_05_25_gull_17_0243_Test5.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat22$species = "lar_gla"
rawdat22$WingID = "17_0243"
rawdat22$TestID = "Test5"
rawdat23 <- read.csv('2020_05_25_gull_17_0243_Test8.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat23$species = "lar_gla"
rawdat23$WingID = "17_0243"
rawdat23$TestID = "Test8"
rawdat24 <- read.csv('2020_05_25_gull_17_0285_Test1.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat24$species = "lar_gla"
rawdat24$WingID = "17_0285"
rawdat24$TestID = "Test1"
rawdat25 <- read.csv('2020_05_25_gull_17_0285_Test2.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
rawdat25$species = "lar_gla"
rawdat25$WingID = "17_0285"
rawdat25$TestID = "Test2"



finaldat <- rbind(rawdat19,rawdat20,rawdat21,rawdat22,rawdat23,rawdat24,rawdat25,rawdat1,rawdat2,rawdat3,rawdat4,rawdat5,rawdat6,rawdat7,rawdat8,rawdat9,rawdat10,rawdat11,rawdat12,rawdat13,rawdat14,rawdat15,rawdat16,rawdat17,rawdat18)
write.csv(finaldat,'/Users/Inman PC/Google Drive/DoctoralThesis/LLT_AIAAPaper/AvianWingLLT/ReorientedWings/2020_05_25_OrientedWings.csv')
