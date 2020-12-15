## This code:
## - reads in digitized points from a full ROM study
## - computes the elbow and manus angles
## - based on bins rounded to the nearest degree

library(tidyverse) # needed for renaming

#------------ Define locations and source functions necessary --------
setwd("/Users/christinaharvey/Google Drive/DoctoralThesis/WingMorphology/")
source('jointangles.R')

allDup <- function (value)# found online by jholtman at gmail.com
{duplicated(value) | duplicated(value, fromLast = TRUE)}

#-------- Find all file names necessary -------
setwd("/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/03_processed_optitrack")
filename.3d   <- list.files(pattern = paste('*',".csv",sep=""))
max_sample_no = 3 # maximum amount of samples to have within one bin
bin_size      = 2 # #degx#deg bins that will have the max amount of samples

### ----- Step 1: Compute elbow and manus angles -------------
for (i in 1:length(filename.3d)){
  #-- Read in data
  setwd("/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/03_processed_optitrack")
  xyz_one <- read.csv(filename.3d[i], stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA") )
  xyz_one <- as_tibble(xyz_one)
  if (ncol(xyz_one)== 35){
    xyz_one <- rename(xyz_one,
                      # new name = old name
                      frameID = frame,
                      pt1_X = humerus_position_x,          pt1_Y = humerus_position_y,          pt1_Z = humerus_position_z,
                      pt2_X = elbow_position_x,            pt2_Y = elbow_position_y,            pt2_Z = elbow_position_z,
                      pt3_X = wrist_position_x,            pt3_Y = wrist_position_y,            pt3_Z = wrist_position_z,
                      pt4_X = cmc_position_x,              pt4_Y = cmc_position_y,              pt4_Z = cmc_position_z,
                      pt6_X = wrist_leading_position_x,    pt6_Y = wrist_leading_position_y,    pt6_Z = wrist_leading_position_z,
                      pt7_X = cmc_leading_position_x,      pt7_Y = cmc_leading_position_y,      pt7_Z = cmc_leading_position_z,
                      pt8_X = p10_position_x,              pt8_Y = p10_position_y,              pt8_Z = p10_position_z,
                      pt9_X = p7_position_x,               pt9_Y = p7_position_y,               pt9_Z = p7_position_z,
                      pt10_X = s1_position_x,              pt10_Y = s1_position_y,              pt10_Z = s1_position_z,
                      pt11_X = s10_position_x,             pt11_Y = s10_position_y,             pt11_Z = s10_position_z,
                      pt12_X = humerus_leading_position_x, pt12_Y = humerus_leading_position_y, pt12_Z = humerus_leading_position_z)
  }else{
    # For wings where the leading positions were too close to the joints to differentiate
    xyz_one <- rename(xyz_one,
                      # new name = old name
                      frameID = frame,
                      pt1_X = humerus_position_x,          pt1_Y = humerus_position_y,          pt1_Z = humerus_position_z,
                      pt2_X = elbow_position_x,            pt2_Y = elbow_position_y,            pt2_Z = elbow_position_z,
                      pt3_X = wrist_position_x,            pt3_Y = wrist_position_y,            pt3_Z = wrist_position_z,
                      pt4_X = cmc_position_x,              pt4_Y = cmc_position_y,              pt4_Z = cmc_position_z,
                      pt8_X = p10_position_x,              pt8_Y = p10_position_y,              pt8_Z = p10_position_z,
                      pt9_X = p7_position_x,               pt9_Y = p7_position_y,               pt9_Z = p7_position_z,
                      pt10_X = s1_position_x,              pt10_Y = s1_position_y,              pt10_Z = s1_position_z,
                      pt11_X = s10_position_x,             pt11_Y = s10_position_y,             pt11_Z = s10_position_z)

    xyz_one$pt12_X = xyz_one$pt1_X
    xyz_one$pt12_Y = xyz_one$pt1_Y
    xyz_one$pt12_Z = xyz_one$pt1_Z

    xyz_one$pt6_X = xyz_one$pt3_X
    xyz_one$pt6_Y = xyz_one$pt3_Y
    xyz_one$pt6_Z = xyz_one$pt3_Z

    xyz_one$pt7_X = xyz_one$pt4_X
    xyz_one$pt7_Y = xyz_one$pt4_Y
    xyz_one$pt7_Z = xyz_one$pt4_Z
  }

  xyz_one         <- xyz_one[complete.cases(xyz_one[,3:ncol(xyz_one)]),]  #Remove NaN rows

  if (nrow(xyz_one) == 0){
    print(i)
    next
  }
  xyz_one$species = paste(tolower(strsplit(filename.3d[i], "_")[[1]][2]),"_",strsplit(filename.3d[i], "_")[[1]][3], sep = "")
  xyz_one$birdid  = paste(strsplit(filename.3d[i], "_")[[1]][4],"_",strsplit(filename.3d[i], "_")[[1]][5], sep = "")
  xyz_one$testid  = strsplit(filename.3d[i], "_")[[1]][7]

  if (i == 1){
    dat_all <- xyz_one[1,c("species","birdid","testid")]
  } else{
    dat_all <- rbind(dat_all,xyz_one[1,c("species","birdid","testid")])
  }

  xyz_one$elbow <- NA
  xyz_one$manus <- NA
  #-- Loop through every row to compute elbow and manus angle
  for (j in 1:nrow(xyz_one)){
    xyz_one$elbow[j]     <- jointangles(xyz_one$pt1_X[j],xyz_one$pt1_Y[j],xyz_one$pt1_Z[j],
                                             xyz_one$pt2_X[j],xyz_one$pt2_Y[j],xyz_one$pt2_Z[j],
                                             xyz_one$pt3_X[j],xyz_one$pt3_Y[j],xyz_one$pt3_Z[j])
    xyz_one$manus[j]     <- jointangles(xyz_one$pt2_X[j],xyz_one$pt2_Y[j],xyz_one$pt2_Z[j],
                                             xyz_one$pt3_X[j],xyz_one$pt3_Y[j],xyz_one$pt3_Z[j],
                                             xyz_one$pt4_X[j],xyz_one$pt4_Y[j],xyz_one$pt4_Z[j])

  }
  # Bin the data by 1 deg increments
  xyz_one$elbow_round = bin_size*round(xyz_one$elbow/bin_size)
  xyz_one$manus_round = bin_size*round(xyz_one$manus/bin_size)

  # only keep the unique values
  dat_subsample <- xyz_one[!allDup(xyz_one[,c("elbow_round","manus_round")]),]
  # save all duplicated rows so that we can randomly select one value for each bin and add back
  dat_dup       <- xyz_one[allDup(xyz_one[,c("elbow_round","manus_round")]),]

  if (nrow(dat_dup) > 0 ){
    dat_dup_uni   <- unique(dat_dup[,c("elbow_round","manus_round")])
    # Loop through each duplicated values and save one random sample
    for (j in 1:nrow(dat_dup_uni)){
      # limit to the current unique bin
      tmp = subset(dat_dup, elbow_round == dat_dup_uni$elbow_round[j] &
                     manus_round == dat_dup_uni$manus_round[j])
      # select at maximum five samples per bin
      if (nrow(tmp) > max_sample_no){
        tmp = subset(tmp, frameID %in% sample(tmp$frameID, max_sample_no))
      }

      dat_subsample <- rbind(dat_subsample,tmp)
    }
  }


  # ---- Save output data ----
  setwd("/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/05_subsampled_optitrack")
  filename <- paste("2020_12_15_",xyz_one$species[1],"_",xyz_one$birdid[1],"_",xyz_one$testid[1],"_subsampled.csv",sep = "")
  write.csv(file = filename,dat_subsample, row.names = FALSE)
}

write.csv(file = "2020_12_15_IDfile.csv",dat_all, row.names = FALSE)
