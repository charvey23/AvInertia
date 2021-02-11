# ---- This code is to be used with the LLT takes raw digitized points and prepares them for testing
#When re-running for new data double check that the transformations are working correctly.
library(stringr)
# CAUTION: This code is now set up to orient wings as follows:
#          1. Pt3 (Wrist) will be inline with Pt1 (humeral head) along the y axis
#          2. Pt3 (Wrist) will be inline with Pt1 (humeral head) along the z axis
#          3. Pt3 (Wrist) will be inline with Pt10 (S1) along the x axis

# Created: Christina Harvey
# Last updated: 16-Dec-2020

setwd("/Users/christinaharvey/Google Drive/DoctoralThesis/WingMorphology/")
source('jointangles.R')
library(stringr)

# ------------------------- Set file directory -------------------------
setwd("/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/05_subsampled_optitrack")

# ----------------- Define set variables ------------
point_list <- c("pt1","pt2","pt3","pt4","pt6","pt7","pt8","pt9","pt10","pt11","pt12")
dim_list   <- c("X","Y","Z")
xaxis <- c(1, 0, 0)
yaxis <- c(0, 1, 0)
zaxis <- c(0, 0, 1)

# ------------------------- Read in data -------------------------
dat_info   <- read.csv(file = "2020_12_15_IDfile.csv")

for (i in 1:nrow(dat_info)){
  # Read in specific wing data
  filename <- paste("2020_12_15_",dat_info$species[i],"_",dat_info$birdid[i],"_00",dat_info$testid[i],"_subsampled.csv",sep = "")
  wing = str_sub(dat_info$birdid[i],-1)

  rawdat   <- read.csv(filename, stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA") )

  if (i < 12){
    rawdat[,3:35] <- rawdat[,3:35]/1000}  # First few are in mm

  # check ulna/radius length to be sure on units
  dat_info$ur_length[i] = mean(sqrt((rawdat$pt3_X - rawdat$pt2_X)^2 + (rawdat$pt3_Y - rawdat$pt2_Y)^2 + (rawdat$pt3_Z - rawdat$pt2_Z)^2))

  ## ------- Reorient Wings ---------
  dat_clean <- rawdat
  #------------------------------- Step 1 -------------------------------
  ##### error doesn't adjust all points ######
  #Make pt 1 on humerus the beginning location
  for (j in 2:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      pt1_name = paste(point_list[1],dim_list[k],sep = "_")
      dat_clean[,col_name] <- dat_clean[,col_name] - dat_clean[,pt1_name]
    }
  }
  dat_clean[,c("pt1_X","pt1_Y","pt1_Z")] <- 0

  #-------------------------------  Step 2 -------------------------------
  # Set pt 3 in line with Pt 1 along the wingspan
  #Calculate angle between axis that runs along wing length and pt 3 (rotate about z)
  #Project onto xaxis & make negative if j1 > 0:
  i1 = dat_clean$pt3_X
  j1 = dat_clean$pt3_Y

  dot_xaxis = i1
  vec_norm  = sqrt(i1^2 + j1^2)
  interior  = dot_xaxis/vec_norm
  thetaz    = acos(interior) + pi / 2
  thetaz[which(j1 > 0)] = -acos(interior[which(j1 > 0)]) + pi / 2

  # Rotate about the z axis
  tmp <- dat_clean
  for (j in 1:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      if (k == 1){
        dat_clean[,col_name] = cos(thetaz) * tmp[,col_name] - sin(thetaz) * tmp[,paste(point_list[j],dim_list[2],sep = "_")]
      } else{
        if(k == 2){
          dat_clean[,col_name] = sin(thetaz) * tmp[,paste(point_list[j],dim_list[1],sep = "_")] + cos(thetaz) * tmp[,col_name]
        }
      }
    }
  }
  # CAUTION: Verify that Pt3X=x_adj1=0 after this step - View(subset(pts,Pt.No == 3))

  #------------------------------- Step 3 -------------------------------
  # Set pt 3 in line with Pt 1 looking down the length of wing
  # Calculate angle between x axis and pt 3 (rotate about y)
  # Project onto xaxis & make negative if k2 > 0
  j2 = dat_clean$pt3_Y
  k2 = dat_clean$pt3_Z

  dot_yaxis = j2
  vec_norm  = sqrt(j2^2 + k2^2)
  interior = dot_yaxis/vec_norm

  thetax = acos(interior)
  thetax[which(k2 > 0)] = -acos(interior[which(k2 > 0)])

  # Rotate about the x axis
  tmp <- dat_clean
  for (j in 1:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      if (k == 2){
        dat_clean[,col_name] = cos(thetax) * tmp[,col_name] - sin(thetax) * tmp[,paste(point_list[j],dim_list[3],sep = "_")]
      } else{
        if(k == 3){
          dat_clean[,col_name] = sin(thetax) * tmp[,paste(point_list[j],dim_list[2],sep = "_")] + cos(thetax) * tmp[,col_name]
        }
      }
    }
  }

  # CAUTION: Verify that Pt3Z=z_adj2=0 after this step - View(subset(pts,Pt.No == 3))

  #------------------------------- Step 4 -------------------------------
  # Set pt 10 in line with Pt 3 looking from back of wing
  # Calculate angle between x axis and pt 10 (rotate about y)
  # Project onto yaxis & make negative if k2 > 0
  i3 = dat_clean$pt10_X
  k3 = dat_clean$pt10_Z

  dot_xaxis = i3
  vec_norm  = sqrt(i3^2 + k3^2)
  interior  = dot_xaxis/vec_norm

  thetay = -acos(interior)
  thetay[which(k3 > 0)] = acos(interior[which(k3 > 0)])
  #thetay[which(dat_clean$pt2_X < 0)] = thetay[which(dat_clean$pt2_X < 0)] + pi

  # Rotate about the y axis
  tmp <- dat_clean
  for (j in 1:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      if (k == 1){
        dat_clean[,col_name] = cos(thetay) * tmp[,col_name] + sin(thetay) * tmp[,paste(point_list[j],dim_list[3],sep = "_")]
      } else{
        if(k == 3){
          dat_clean[,col_name] = -sin(thetay) * tmp[,paste(point_list[j],dim_list[1],sep = "_")] + cos(thetay) * tmp[,col_name]
        }
      }
    }
  }
  # CAUTION: Verify that Pt10Z=z_adj3=0 after this step - View(subset(pts,Pt.No == 10))


  #------------------------------- Step 5 -------------------------------
  for (j in 1:11) {
    col_name = paste(point_list[j], dim_list[1], sep = "_")
    dat_clean[, col_name] <- -dat_clean[, col_name]
    if(wing == "R"){
      col_name = paste(point_list[j], dim_list[3], sep = "_")
      dat_clean[, col_name] <- -dat_clean[, col_name]
    }
  }

  # ------------------------------ Save data
  filename_new <- paste("2020_01_14_",dat_info$species[i],"_",dat_info$birdid[i],"_00",dat_info$testid[i],"_aligned.csv",sep = "")
  write.csv(dat_clean,filename_new)
}

# Visualize the wings as required - For each calibration verify that the axis is RH
m = 1:nrow(dat_clean)
max = 0.5
plot(dat_clean$pt2_Y[m],dat_clean$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(dat_clean$pt3_Y[m],dat_clean$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(dat_clean$pt4_Y[m], dat_clean$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(dat_clean$pt1_Y[m],dat_clean$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(dat_clean$pt8_Y[m],dat_clean$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(dat_clean$pt11_Y[m],dat_clean$pt11_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")

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

plot(rawdat$pt2_Z[m],rawdat$pt2_X[m], xlim =c(-max,max), ylim = c(-max,max))
points(rawdat$pt3_Z[m],rawdat$pt3_X[m],xlim =c(-max,max), ylim = c(-max,max),col = "blue")
points(rawdat$pt4_Z[m], rawdat$pt4_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "green")
points(rawdat$pt1_Z[m],rawdat$pt1_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "red")
points(rawdat$pt8_Z[m],rawdat$pt8_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "purple")
points(rawdat$pt9_Z[m],rawdat$pt9_X[m], xlim =c(-max,max), ylim = c(-max,max), col = "yellow")
