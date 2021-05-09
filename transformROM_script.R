##### This scripts reads in all range of motion tracks for a specific species.
##### Resizes down to each individual, subsamples and finally reorients each wing shape

## Load required libraries
library(tidyverse)
library(Morpho)
library(geomorph)
library(stringr)
library(gtools)
library(pracma)

source('/Users/christinaharvey/Documents/birdmoment/transformROM_supportingfunctions.R')

### -------------- Read in all morphological measurements ------------
run_data_path = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/rundata/"
output_path   = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/05_transformed_optitrack/"

# all of the non-wing based measurements for all specimens
dat_bird     = readxl::read_xlsx(paste(run_data_path,"bird_measurements_readytorun.xlsx",sep= ""), sheet = 'Major body parts')
# all of the wing based measurements for all specimens
dat_wingspec = readxl::read_xlsx(paste(run_data_path,"bird_measurements_readytorun.xlsx",sep= ""), sheet = 'Within wings')

dat_bird$species = NA
dat_bird$BirdID  = NA
for (i in 1:nrow(dat_bird)){
  if (strsplit(dat_bird$bird_id, "_")[[i]][1] == "COLLI"){
    dat_bird$species[i] = "col_liv"
    dat_bird$BirdID[i]  = paste(strsplit(dat_bird$bird_id, "_")[[i]][2],strsplit(dat_bird$bird_id, "_")[[i]][3], sep = "_")
  }else{
    dat_bird$species[i] = paste(tolower(strsplit(dat_bird$bird_id, "_")[[i]][1]),strsplit(dat_bird$bird_id, "_")[[i]][2], sep = "_")
    dat_bird$BirdID[i]  = paste(strsplit(dat_bird$bird_id, "_")[[i]][3],strsplit(dat_bird$bird_id, "_")[[i]][4], sep = "_")
  }
}

dat_wingspec$species = NA
dat_wingspec$BirdID  = NA
for (i in 1:nrow(dat_wingspec)){
  if (strsplit(dat_wingspec$bird_id, "_")[[i]][1] == "COLLI"){
    dat_wingspec$species[i] = "col_liv"
    dat_wingspec$BirdID[i]  = paste(strsplit(dat_wingspec$bird_id, "_")[[i]][2],strsplit(dat_wingspec$bird_id, "_")[[i]][3], sep = "_")
  }else{
    dat_wingspec$species[i] = paste(tolower(strsplit(dat_wingspec$bird_id, "_")[[i]][1]),strsplit(dat_wingspec$bird_id, "_")[[i]][2], sep = "_")
    dat_wingspec$BirdID[i]  = paste(strsplit(dat_wingspec$bird_id, "_")[[i]][3],strsplit(dat_wingspec$bird_id, "_")[[i]][4], sep = "_")
  }
}

#-------- Find all file names necessary -------
working_path  = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/03_processed_optitrack/"
filename.3d   = list.files(path = working_path, pattern = paste('*',".csv",sep=""))
max_sample_no = 1 # maximum amount of samples to have within one bin
bin_size      = 2 # #deg x #deg bins that will have the max amount of samples

# create blank data frame
dat_ID = as.data.frame(matrix(nrow = length(filename.3d), ncol = 4))
colnames(dat_ID) = c("species","BirdID","TestID","wing")

## ----------------------------------------------------
## --------- Save identifying information -------------
## ----------------------------------------------------

# Save all the identifying information
for (i in 1:length(filename.3d)){
  dat_ID$species[i] = paste(tolower(strsplit(filename.3d[i], "_")[[1]][2]),"_",strsplit(filename.3d[i], "_")[[1]][3], sep = "")
  dat_ID$BirdID[i]  = paste(strsplit(filename.3d[i], "_")[[1]][4],"_",str_sub(strsplit(filename.3d[i], "_")[[1]][5],start = 1, end = -2), sep = "")
  dat_ID$TestID[i]  = as.numeric(strsplit(filename.3d[i], "_")[[1]][7])
  dat_ID$wing[i]    = str_sub(strsplit(filename.3d[i], "_")[[1]][5],-1)
}

no_species = length(unique(dat_ID$species))

## ----------------------------------------------------
## --------- Iterate through each species -------------
## ----------------------------------------------------

for (i in 1:no_species){
  curr_species  = unique(dat_ID$species)[i]

  # identify which filenames belong to the current species
  files_to_read = which(dat_ID$species == curr_species)
  # identify how many body measurements we have for the current species
  dat_body_curr = subset(dat_bird, use_body == "Y" & species == curr_species)

    # Read in all optitrack data from the given species
  dat_raw         = read.csv(paste(working_path,filename.3d[files_to_read[1]],sep = ""), stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA"))
  dat_raw$BirdID  = dat_ID$BirdID[files_to_read[1]]
  dat_raw$TestID  = dat_ID$TestID[files_to_read[1]]
  dat_raw$wing    = dat_ID$wing[files_to_read[1]]
  # First few are in mm adjust to cm
  if (files_to_read[1] < 13){dat_raw[,3:35] <- dat_raw[,3:35]*0.001}

  for (j in 2:length(files_to_read)){
    tmp = read.csv(paste(working_path,filename.3d[files_to_read[j]],sep = ""), stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA"))
    tmp$BirdID  = dat_ID$BirdID[files_to_read[j]]
    tmp$TestID  = dat_ID$TestID[files_to_read[j]]
    tmp$wing    = dat_ID$wing[files_to_read[j]]
    # First few are in mm adjust to cm
    if (files_to_read[j] < 13){tmp[,3:35] <- tmp[,3:35]*0.001}
    # The storm petrel had one wing track all points but not the other
    if (ncol(tmp) > ncol(dat_raw)){
      dat_raw$humerus_leading_position_x = dat_raw$humerus_position_x
      dat_raw$humerus_leading_position_y = dat_raw$humerus_position_y
      dat_raw$humerus_leading_position_z = dat_raw$humerus_position_z

      dat_raw$wrist_leading_position_x = dat_raw$wrist_position_x
      dat_raw$wrist_leading_position_y = dat_raw$wrist_position_y
      dat_raw$wrist_leading_position_z = dat_raw$wrist_position_z

      dat_raw$cmc_leading_position_x = dat_raw$cmc_position_x
      dat_raw$cmc_leading_position_y = dat_raw$cmc_position_y
      dat_raw$cmc_leading_position_z = dat_raw$cmc_position_z
    }
    dat_raw = smartbind(dat_raw, tmp)
  }

  dat_raw$species = curr_species
  col_dat_raw = c("FrameID","time_sec","pt4_X","pt4_Y","pt4_Z" , "pt7_X","pt7_Y","pt7_Z","pt2_X","pt2_Y","pt2_Z","pt1_X","pt1_Y","pt1_Z",
                  "pt12_X","pt12_Y","pt12_Z","pt8_X","pt8_Y","pt8_Z","pt9_X","pt9_Y","pt9_Z","pt10_X","pt10_Y","pt10_Z","pt11_X","pt11_Y",
                  "pt11_Z","pt3_X","pt3_Y","pt3_Z","pt6_X","pt6_Y","pt6_Z","BirdID","TestID","wing","species" )


  ### ---------------- Rename the columns appropriately -----------------
  dat_raw <- as_tibble(dat_raw)
  dat_raw <- rename(dat_raw,
                    FrameID = frame,
                    pt1_X = humerus_position_x,          pt1_Y = humerus_position_y,          pt1_Z = humerus_position_z,
                    pt2_X = elbow_position_x,            pt2_Y = elbow_position_y,            pt2_Z = elbow_position_z,
                    pt3_X = wrist_position_x,            pt3_Y = wrist_position_y,            pt3_Z = wrist_position_z,
                    pt4_X = cmc_position_x,              pt4_Y = cmc_position_y,              pt4_Z = cmc_position_z,
                    pt8_X = p10_position_x,              pt8_Y = p10_position_y,              pt8_Z = p10_position_z,
                    pt9_X = p7_position_x,               pt9_Y = p7_position_y,               pt9_Z = p7_position_z,
                    pt10_X = s1_position_x,              pt10_Y = s1_position_y,              pt10_Z = s1_position_z,
                    pt11_X = s10_position_x,             pt11_Y = s10_position_y,             pt11_Z = s10_position_z)
  # For wings where the leading positions were too close to the joints to differentiate set the joint positions as the edge points as well
  if (ncol(dat_raw) == 39){
    dat_raw <- rename(dat_raw,
                      pt6_X = wrist_leading_position_x,    pt6_Y = wrist_leading_position_y,    pt6_Z = wrist_leading_position_z,
                      pt7_X = cmc_leading_position_x,      pt7_Y = cmc_leading_position_y,      pt7_Z = cmc_leading_position_z,
                      pt12_X = humerus_leading_position_x, pt12_Y = humerus_leading_position_y, pt12_Z = humerus_leading_position_z)
  }else{
    dat_raw$pt12_X = dat_raw$pt1_X
    dat_raw$pt12_Y = dat_raw$pt1_Y
    dat_raw$pt12_Z = dat_raw$pt1_Z

    dat_raw$pt6_X = dat_raw$pt3_X
    dat_raw$pt6_Y = dat_raw$pt3_Y
    dat_raw$pt6_Z = dat_raw$pt3_Z

    dat_raw$pt7_X = dat_raw$pt4_X
    dat_raw$pt7_Y = dat_raw$pt4_Y
    dat_raw$pt7_Z = dat_raw$pt4_Z

  }

  dat_raw       = as.data.frame(dat_raw[complete.cases(dat_raw[,3:ncol(dat_raw)]),])  #Remove NaN rows
  dat_raw       = dat_raw[,c(col_dat_raw)] # re-order accordingly will work as long as the first species tracks all points

  # ---- Calculate the elbow and wrist angle for each frame
  dat_raw$elbow <- NA
  dat_raw$manus <- NA
  #-- Loop through every row to compute elbow and manus angle
  for (j in 1:nrow(dat_raw)){
    dat_raw$elbow[j]     <- jointangles(dat_raw$pt1_X[j],dat_raw$pt1_Y[j],dat_raw$pt1_Z[j],
                                        dat_raw$pt2_X[j],dat_raw$pt2_Y[j],dat_raw$pt2_Z[j],
                                        dat_raw$pt3_X[j],dat_raw$pt3_Y[j],dat_raw$pt3_Z[j])
    dat_raw$manus[j]     <- jointangles(dat_raw$pt2_X[j],dat_raw$pt2_Y[j],dat_raw$pt2_Z[j],
                                        dat_raw$pt3_X[j],dat_raw$pt3_Y[j],dat_raw$pt3_Z[j],
                                        dat_raw$pt4_X[j],dat_raw$pt4_Y[j],dat_raw$pt4_Z[j])

  }
  # Bin the data
  dat_raw$elbow_round = bin_size*round(dat_raw$elbow/bin_size)
  dat_raw$manus_round = bin_size*round(dat_raw$manus/bin_size)

  dat_raw$BirdID_FrameSpec = dat_raw$BirdID

  ## ----------------------------------------------------
  ## --------- Iterate through each specimen -------------
  ## ----------------------------------------------------

  for (ind in 1:nrow(dat_body_curr)){

    curr_BirdID = dat_body_curr$BirdID[ind] # save the current BirdID

    col_char = c("species","BirdID","TestID","BirdID_FrameSpec","FrameID","time_sec","wing","elbow","manus","elbow_round","manus_round")
    col_all  = c(col_char, colnames(dat_raw[,3:35]))

    #### ------- Resize the applicable wings --------
    for (j in 1:nrow(dat_body_curr)){

      # there are gull specific adjustments since wing doesn't match any body
      if (dat_body_curr$BirdID[j] != curr_BirdID | curr_species == "lar_gla"){

        if (curr_species == "lar_gla" & dat_body_curr$BirdID[j]== "20_0341"){
          adjust = subset(dat_raw, species == curr_species & BirdID == "21_0310")

          target_bone_len    = subset(dat_wingspec, species == curr_species & BirdID == "20_0341")$ulna_length_mm*0.001
          adjust_bone_length = mean(calc_dist(adjust[,c(9:11,30:32)]))
        } else{
          target = subset(dat_raw, species == curr_species & BirdID == curr_BirdID)
          adjust = subset(dat_raw, species == curr_species & BirdID == dat_body_curr$BirdID[j])
          # if there is no data for the wing ROM move onto the next wing
          if(nrow(adjust) == 0){next}
          target_bone_len    = subset(dat_wingspec, species == curr_species & BirdID == curr_BirdID)$humerus_length_mm*0.001
          adjust_bone_length = subset(dat_wingspec, species == curr_species & BirdID == dat_body_curr$BirdID[j])$humerus_length_mm*0.001
        }

        tmp = resize_to(specimen_to_adjust = adjust, char_colnames = col_char,
                        adjust_length = adjust_bone_length, target_length = target_bone_len)

        colnames(tmp) = col_all
        # Check:
        # scale_factor = target_bone_len/adjust_bone_length # length of target bone/length of bone in wing that will be resized
        # mean(calc_dist(tmp[,c(18:20,39:41)]))/mean(calc_dist(adjust[,c(9:11,30:32)])) == scale_factor
      }else {
        tmp = subset(dat_raw, species == curr_species & BirdID == curr_BirdID)
        tmp = tmp[,col_all]
      }

      if(!exists("dat_resized")){
        dat_resized = tmp} else{
        dat_resized = rbind(dat_resized,tmp)}
    }
    # Adjust the BirdID to reflect the current bird - BirdID_FrameSpec indicates which bird the frame originally came from
    dat_resized$BirdID = curr_BirdID

    #### ------- Subsample the current data frame --------
    # only keep the unique values
    dat_keep <- dat_resized[!allDup(dat_resized[,c("elbow_round","manus_round")]),]
    # save all duplicated rows so that we can randomly select one value for each bin and add back
    dat_dup  <- dat_resized[allDup(dat_resized[,c("elbow_round","manus_round")]),]

    dat_subsample = subsample_frames(dat_keep,dat_dup,curr_BirdID,max_sample_no)

    #### ------- Reorient the wings --------
    dat_complete = reorient_wings(dat_subsample)

    dat_complete$S_proj = NA
    dat_complete$S = NA
    for (i in 1:nrow(dat_complete)){
      ## ----- Calculate the projected wing area ------
      # this is the correct order because X is negative and Y is positive
      x_vertices = c(dat_complete$pt6_X[i],dat_complete$pt7_X[i],dat_complete$pt8_X[i],dat_complete$pt9_X[i],dat_complete$pt10_X[i],dat_complete$pt11_X[i],dat_complete$pt12_X[i])
      y_vertices = c(dat_complete$pt6_Y[i],dat_complete$pt7_Y[i],dat_complete$pt8_Y[i],dat_complete$pt9_Y[i],dat_complete$pt10_Y[i],dat_complete$pt11_Y[i],dat_complete$pt12_Y[i])
      dat_complete$S_proj[i] <- polyarea(x_vertices, y_vertices)
      ## ----- Calculate the total wing area ------
      A1 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt11_X","pt11_Y","pt11_Z")]-dat_complete[i,c("pt2_X","pt2_Y","pt2_Z")])),
                          as.vector(t(dat_complete[i,c("pt11_X","pt11_Y","pt11_Z")]-dat_complete[i,c("pt12_X","pt12_Y","pt12_Z")]))), p = 2)
      A2 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt2_X","pt2_Y","pt2_Z")]-dat_complete[i,c("pt6_X","pt6_Y","pt6_Z")])),
                          as.vector(t(dat_complete[i,c("pt2_X","pt2_Y","pt2_Z")]-dat_complete[i,c("pt12_X","pt12_Y","pt12_Z")]))), p = 2)
      A3 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt11_X","pt11_Y","pt11_Z")]-dat_complete[i,c("pt10_X","pt10_Y","pt10_Z")])),
                          as.vector(t(dat_complete[i,c("pt11_X","pt11_Y","pt11_Z")]-dat_complete[i,c("pt2_X","pt2_Y","pt2_Z")]))), p = 2)
      A4 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt10_X","pt10_Y","pt10_Z")]-dat_complete[i,c("pt6_X","pt6_Y","pt6_Z")])),
                          as.vector(t(dat_complete[i,c("pt10_X","pt10_Y","pt10_Z")]-dat_complete[i,c("pt2_X","pt2_Y","pt2_Z")]))), p = 2)
      A5 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt10_X","pt10_Y","pt10_Z")]-dat_complete[i,c("pt9_X","pt9_Y","pt9_Z")])),
                          as.vector(t(dat_complete[i,c("pt10_X","pt10_Y","pt10_Z")]-dat_complete[i,c("pt6_X","pt6_Y","pt6_Z")]))), p = 2)
      A6 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt9_X","pt9_Y","pt9_Z")]-dat_complete[i,c("pt7_X","pt7_Y","pt7_Z")])),
                          as.vector(t(dat_complete[i,c("pt9_X","pt9_Y","pt9_Z")]-dat_complete[i,c("pt6_X","pt6_Y","pt6_Z")]))), p = 2)
      A7 = 0.5*Norm(cross(as.vector(t(dat_complete[i,c("pt9_X","pt9_Y","pt9_Z")]-dat_complete[i,c("pt8_X","pt8_Y","pt8_Z")])),
                          as.vector(t(dat_complete[i,c("pt9_X","pt9_Y","pt9_Z")]-dat_complete[i,c("pt7_X","pt7_Y","pt7_Z")]))), p = 2)

      dat_complete$S[i] = sum(A1,A2,A3,A4,A5,A6,A7)
    }
    # ------------------------------ Save data
    filename_new <- paste(output_path,format(Sys.Date(), "%Y_%m_%d"),"_",curr_species,"_",curr_BirdID,"_transformed.csv",sep = "")
    write.csv(dat_complete,filename_new)
    remove(dat_resized) # must be remove to avoid saving from different species
  } # end of the specimen loop
} # end of the species loop
filename_new <- paste(run_data_path,format(Sys.Date(), "%Y_%m_%d"),"_IDfile.csv",sep = "")
write.csv(dat_ID,filename_new)

