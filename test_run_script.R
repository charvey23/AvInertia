devtools::load_all()

# --------------------- Read in data -----------------------
# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/rundata/"
# identification info on each of the individual specimens
dat_ind      = read.csv(file = paste(path_data_folder,"2020_12_15_IDfile.csv",sep= ""))
# all of the non-wing based measurements for all specimens - NEED TO RESAVE THIS - ALSO NEED TO CORRECT ONE ID NUMBER TO MATCH FEATHERS
dat_bird     = readxl::read_xlsx(paste(path_data_folder,"2020_12_19_bird_measurements.xlsx",sep= ""), sheet = 'Major body parts')

# all of the wing based measurements for all specimens
dat_wingspec = readxl::read_xlsx(paste(path_data_folder,"2020_12_19_bird_measurements.xlsx",sep= ""), sheet = 'Within wings')

# all of the feather measurements for one specimen of each species
dat_feat     = readxl::read_xlsx(paste(path_data_folder,"2021_01_13_trackingsheet.xlsx",sep= ""), sheet = 'FeatherData')

# all of the barb measurements for each species baed on measurements from previous studies
dat_barb     = readxl::read_xlsx(paste(path_data_folder,"2021_01_13_trackingsheet.xlsx",sep= ""), sheet = 'BarbInfo')

# ----------------- Clean Data -----------------------
# -- Clean up the data set relating to each individual specimen --

# Save info about whether this is a right or left wing
dat_ind$wing_side <- stringr::str_sub(dat_ind$birdid, -1)
dat_ind$BirdID    <- stringr::str_sub(dat_ind$birdid, 1, stringr::str_length(dat_ind$birdid)-1)

# -- Clean up the data set relating to full bird measurements of each individual specimen--
dat_bird$species <- NA
dat_bird$BirdID  <- NA
for (i in 1:nrow(dat_bird)){
  if (strsplit(dat_bird$bird_id, "_")[[i]][1] == "COLLI"){
    dat_bird$species[i] = "col_liv"
    dat_bird$BirdID[i]  = paste(strsplit(dat_bird$bird_id, "_")[[i]][2],strsplit(dat_bird$bird_id, "_")[[i]][3], sep = "_")
  }else{
    dat_bird$species[i] = paste(tolower(strsplit(dat_bird$bird_id, "_")[[i]][1]),strsplit(dat_bird$bird_id, "_")[[i]][2], sep = "_")
    dat_bird$BirdID[i]  = paste(strsplit(dat_bird$bird_id, "_")[[i]][3],strsplit(dat_bird$bird_id, "_")[[i]][4], sep = "_")
  }
}
# -- Clean up the data set relating to the wing specific measurements for each individual specimen --
dat_wingspec$species <- NA
dat_wingspec$BirdID  <- NA
for (i in 1:nrow(dat_wingspec)){
  if (strsplit(dat_wingspec$bird_id, "_")[[i]][1] == "COLLI"){
    dat_wingspec$species[i] = "col_liv"
    dat_wingspec$BirdID[i]  = paste(strsplit(dat_wingspec$bird_id, "_")[[i]][2],strsplit(dat_wingspec$bird_id, "_")[[i]][3], sep = "_")
  }else{
    dat_wingspec$species[i] = paste(tolower(strsplit(dat_wingspec$bird_id, "_")[[i]][1]),strsplit(dat_wingspec$bird_id, "_")[[i]][2], sep = "_")
    dat_wingspec$BirdID[i]  = paste(strsplit(dat_wingspec$bird_id, "_")[[i]][3],strsplit(dat_wingspec$bird_id, "_")[[i]][4], sep = "_")
  }
}

# Merge the wing specific measurements with the full bird measurements
dat_bird <- merge(dat_bird,dat_wingspec, by = c("species","BirdID"))

# Correct all the units in the data frame to SI
for (i in 1:length(colnames(dat_bird))){
  # Adjust all masses to be in kg
  if (grepl("mass_g",colnames(dat_bird)[i],fixed=TRUE)){
    dat_bird[,i] = dat_bird[,i]/1000
  }
  # Adjust all lengths to be in m
  if (grepl("_cm",colnames(dat_bird)[i],fixed=TRUE)){
    dat_bird[,i] = dat_bird[,i]/100
  }
  # Adjust all lengths to be in m
  if (grepl("_mm",colnames(dat_bird)[i],fixed=TRUE)){
    dat_bird[,i] = dat_bird[,i]/1000
  }
}

#Rename columns as necessary
names(dat_bird)[names(dat_bird) == "humerus_muscles_mass_g"] <- "brachial_muscle_mass"
names(dat_bird)[names(dat_bird) == "raduln_muscles_mass_g"]  <- "antebrachial_muscle_mass"
names(dat_bird)[names(dat_bird) == "cmc_muscles_mass_g"]     <- "manus_muscle_mass"
names(dat_bird)[names(dat_bird) == "whole_body_mass_g"]      <- "total_bird_mass"
names(dat_bird)[names(dat_bird) == "tert_feathers_mass_g"]   <- "tertiary_mass"
names(dat_bird)[names(dat_bird) == "skin_coverts_mass_g"]    <- "all_skin_coverts_mass"

names(dat_bird)[names(dat_bird) == "head_length_cm"]             <- "head_length"
names(dat_bird)[names(dat_bird) == "head_mass_g"]                <- "head_mass"
names(dat_bird)[names(dat_bird) == "head_height_max_cm"]         <- "head_height"
names(dat_bird)[names(dat_bird) == "neck_mass_g"]                <- "neck_mass"
names(dat_bird)[names(dat_bird) == "neck_width_cm"]              <- "neck_width"
names(dat_bird)[names(dat_bird) == "neck_length_cm"]             <- "neck_length"
names(dat_bird)[names(dat_bird) == "torsotail_length_cm"]        <- "torsotail_length"
names(dat_bird)[names(dat_bird) == "torsotail_mass_g"]           <- "torsotail_mass"
names(dat_bird)[names(dat_bird) == "tail_length_cm"]             <- "tail_length"
names(dat_bird)[names(dat_bird) == "tail_width_cm"]              <- "tail_width"
names(dat_bird)[names(dat_bird) == "right_leg_mass_g"]           <- "right_leg_mass"
names(dat_bird)[names(dat_bird) == "left_leg_mass_"]             <- "left_leg_mass"
names(dat_bird)[names(dat_bird) == "body_width_max_cm"]          <- "body_width_max"
names(dat_bird)[names(dat_bird) == "width_at_leg_insert_cm"]     <- "body_width_at_leg_insert"
names(dat_bird)[names(dat_bird) == "x_loc_of_body_max_cm"]       <- "x_loc_of_body_max"
names(dat_bird)[names(dat_bird) == "tail_width_cm"]              <- "tail_width"
names(dat_bird)[names(dat_bird) == "x_loc_leg_insertion_cm"]     <- "x_loc_leg_insertion"
names(dat_bird)[names(dat_bird) == "x_loc_TorsotailCoG_cm"]      <- "x_loc_TorsotailCoG"
names(dat_bird)[names(dat_bird) == "z_loc_TorsotailCoG_cm"]      <- "z_loc_TorsotailCoG"
names(dat_bird)[names(dat_bird) == "x_loc_of_humeral_insert_cm"] <- "x_loc_humeral_insert"
names(dat_bird)[names(dat_bird) == "y_loc_of_humeral_insert_cm"] <- "y_loc_humeral_insert"
names(dat_bird)[names(dat_bird) == "z_loc_of_humeral_insert_cm"] <- "z_loc_humeral_insert"
names(dat_bird)[names(dat_bird) == "whole_body_mass_g"]          <- "total_bird_mass"
#Correct the sign of the x measurement
dat_bird$x_loc_TorsotailCoG = - dat_bird$x_loc_TorsotailCoG

# save a singular wing mass - left if that is the one measured right otherwise
dat_bird$wing_mass <-dat_bird$left_wing_mass_g
for (i in 1:nrow(dat_bird)){
  if (dat_bird$right_or_left[i] == "Right"){
    dat_bird$wing_mass[i] = dat_bird$right_wing_mass_g[i]
  }
}

## ----- Adjust the feather input properties ----------

# Update to include the barb information within the species specific info
dat_barb$barb_distance <- dat_barb$barb_distance_um*10^-6
dat_barb$barb_radius   <- dat_barb$barb_radius_um*10^-6
dat_bird <- merge(dat_bird,dat_barb, by=c("species"))
# need the calamus angle to be positive (i.e. from start of vane to tip of calamus)
for (k in 1:nrow(dat_feat[which(dat_feat$Component == "calamus length"),c("Angle")])){
  if (dat_feat[which(dat_feat$Component == "calamus length"),c("Angle")][k,1] < 0){
    dat_feat[which(dat_feat$Component == "calamus length"),c("Angle")][k,1] = 180 + dat_feat[which(dat_feat$Component == "calamus length"),c("Angle")][k,1]
  }
}

# need the vane angle to be negative (i.e. from start of vane to tip of feather)
for (k in 1:nrow(dat_feat[which(dat_feat$Component == "vane length"),c("Angle")])){
  if (dat_feat[which(dat_feat$Component == "vane length"),c("Angle")][k,1] > 0){
    dat_feat[which(dat_feat$Component == "vane length"),c("Angle")][k,1] = dat_feat[which(dat_feat$Component == "vane length"),c("Angle")][k,1] - 180
  }
}


# define all assumed material properties
dat_mat = list()
dat_mat$material  =	c("Muscle", "Bone", "Skin", "Cortex", "Medullary")
dat_mat$density   = c(1100, 2060, 1060, 660, 37)
# CHECK THAT THIS IS THE CORRECT PROPORTION AND LIST ALL REFERENCES

# --------------------- Initialize variables -----------------------
all_data = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
column_names = c("species","BirdID","TestID","FrameID","prop_type","component","value")
colnames(all_data) = column_names
specimens <- unique(dat_ind[,c("species","BirdID")])

# ----------- Iterate through each species ---------
for (k in 1:nrow(specimens)){

  # CAUTION: THIS JUST FOR DEBUGGING
  k = 24

  # ----------- Filter the data to the current species ---------
  species_curr  = specimens$species[k]
  birdid_curr   = specimens$BirdID[k]

  #Read in all the wing configuration data - loop through all files for the individual
  filename_wingconfigs = list.files(path = path_data_folder, pattern = paste(species_curr,birdid_curr,sep="_"))
  dat_wing_curr   = read.csv(paste(path_data_folder,filename_wingconfigs[1],sep=""))
  for (n in 2:length(filename_wingconfigs)){
    dat_wing_curr = rbind(dat_wing_curr,read.csv(paste(path_data_folder,filename_wingconfigs[n],sep="")))
  }
  # Save the side of the wing and the bird ID out of the file path
  dat_wing_curr$wing_side <- stringr::str_sub(dat_wing_curr$birdid, -1)
  dat_wing_curr$BirdID    <- stringr::str_sub(dat_wing_curr$birdid, 1, stringr::str_length(dat_wing_curr$birdid)-1)

  ## ---
  dat_bird_curr = subset(dat_bird, species == species_curr & BirdID == birdid_curr)

  # Create the bone specific data frame
  # Inputs are adjusted to be in SI units, kg and m
  dat_bone_curr = as.data.frame(matrix(nrow = 6, ncol = 5))
  names(dat_bone_curr) <- c("bone", "bone_mass", "bone_len", "bone_out_rad", "bone_in_rad")
  dat_bone_curr$bone <- c("Humerus","Ulna","Radius","Carpometacarpus","Ulnare", "Radiale")
  dat_bone_curr$bone_mass  <- c(dat_bird_curr$humerus_mass_g,dat_bird_curr$ulna_mass_g,dat_bird_curr$radius_mass_g,dat_bird_curr$cmc_mass_g,
                          dat_bird_curr$ulnare_mass_g,dat_bird_curr$radiale_mass_g)
  dat_bone_curr$bone_len  <- c(dat_bird_curr$humerus_length_mm,dat_bird_curr$ulna_length_mm,dat_bird_curr$radius_length_mm,dat_bird_curr$cmc_length_mm,
                             NA,NA)
  dat_bone_curr$bone_out_rad <- c(dat_bird_curr$humerus_diameter_mm, ## ASK ABOUT THIS
                                  dat_bird_curr$ulna_diameter_mm,
                                  dat_bird_curr$radius_diameter_mm,
                                  dat_bird_curr$cmc_diameter_mm,NA,NA)
  dat_bone_curr$bone_out_rad <- 0.5*dat_bone_curr$bone_out_rad
  dat_bone_curr$bone_in_rad  <- 0.78*dat_bone_curr$bone_out_rad # CHECK THAT THIS IS THE CORRECT PROPORTION AND LIST ALL REFERENCES

  # ------------------------ Create the feather specific data frame --------------------------
  # Inputs are adjusted to be in SI units, kg and m
  tmp = subset(dat_feat, species == species_curr & bird_id == birdid_curr)
  tmp_area   <- tidyr::spread(aggregate(x = tmp$Area, by = list(tmp$Feather, tmp$Component), "mean"),"Group.2","x")
  tmp_area   <- tmp_area[,c("Group.1","distal vane","proximal vane")]
  tmp_length <- tidyr::spread(aggregate(x = tmp$Length, by = list(tmp$Feather, tmp$Component), "mean"),"Group.2","x")
  tmp_length   <- tmp_length[,c("Group.1","calamus length","vane length","rachis width")]
  tmp_angle  <- tidyr::spread(aggregate(x = tmp$Angle, by = list(tmp$Feather, tmp$Component), "mean"),"Group.2","x")
  tmp_angle   <- tmp_angle[,c("Group.1","calamus length","vane length")] ## NEED TO ADJUST THESE TO WORK

  dat_feat_curr <- merge(tmp_area,tmp_length, by = "Group.1")

  # computes the supplement angle to the interior angle between the calamus and vane assuming that the calamus is a positive angle and the vane is negative
  dat_feat_curr$vane_angle <- 180 - (tmp_angle$"calamus length" - tmp_angle$"vane length")

  # Change all units from cm to m
  dat_feat_curr$l_vane <- dat_feat_curr$`vane length`*0.01
  dat_feat_curr$l_cal  <- dat_feat_curr$`calamus length`*0.01
  dat_feat_curr$w_cal  <- dat_feat_curr$`rachis width`*0.01
  dat_feat_curr$w_vd   <- 0.01*(dat_feat_curr$"distal vane")/(dat_feat_curr$"vane length") # compute the average width of the distal vane
  dat_feat_curr$w_vp   <- 0.01*(dat_feat_curr$"proximal vane")/(dat_feat_curr$"vane length") # compute the average width of the proximal vane
  # rename last column
  names(dat_feat_curr)[names(dat_feat_curr) == "Group.1"]         <- "feather"

  # save all feather masses
  dat_feat_curr$m_f <- NA
  for (i in 1:nrow(dat_feat_curr)){
    dat_feat_curr$m_f[i] = dat_bird_curr[1,paste(tolower(dat_feat_curr$feather[i]),"mass_g",sep = "_")]
  }
  # save the combined alula feather masses
  row_alula = nrow(dat_feat_curr)+1
  dat_feat_curr[row_alula, ]       <- NA
  dat_feat_curr$feather[row_alula] <- "alula"
  dat_feat_curr$m_f[row_alula]     <- 0 # for now have no alula masses

  dat_feat_curr$species = species_curr
  dat_feat_curr$BirdID  = birdid_curr

  # save all feather data for future use
  if (k == 1){
    dat_feat_all <- dat_feat_curr
  }else{
    dat_feat_all <- rbind(dat_feat_all,dat_feat_curr)
  }

  # --------------- Clean up wing 3D incoming data --------------------
  # CAUTION:   The read in data has the origin at the RH humeral head. Pt1X, Pt1Y, Pt1Z = 0
  #            It is critical that the data is given in the structural
  #            frame of reference with the origin at the vehicle reference
  #            point (VRP). VRP is assumed to be in the center of the body
  #            on the y axis between the two humerus bones at the clavicle point

  dat_pt = dat_wing_curr[,4:36]

  for (i in 1:length(colnames(dat_pt))){
    # x position - the sign has been verified
    if (grepl("X",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] - (dat_bird_curr$x_loc_humeral_insert)
    }
    # y position - the sign has been verified
    # CAUTION: this assumes the wings is on the right side of the bird for all calculations
    if (grepl("Y",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] + (dat_bird_curr$y_loc_humeral_insert)
    }
    # z position - the sign has been verified + positive ensures that the shoulder is "below" the clavicle origin
    if (grepl("Z",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] + (dat_bird_curr$z_loc_humeral_insert)
    }
  }

  # -------------- Iterate through the wings of this species ------------------------------

  for (ind_wing in 1:length(dat_wing_curr$frameID)){
    ind_wing = 1 # CAUTION: FOR  DEBUGGING ONLY
    # both dataframes below should only be one row of input points
    dat_id_curr = dat_wing_curr[ind_wing,c("species","BirdID","testid","frameID")]
    names(dat_id_curr)[names(dat_id_curr) == "frameID"] <- "FrameID"
    names(dat_id_curr)[names(dat_id_curr) == "testid"] <- "TestID"
    dat_pt_curr = dat_pt[ind_wing,]

    # Initialize common pts
    Pt1  = c(dat_pt_curr$pt1_X, dat_pt_curr$pt1_Y, dat_pt_curr$pt1_Z) # Shoulder
    Pt2  = c(dat_pt_curr$pt2_X, dat_pt_curr$pt2_Y, dat_pt_curr$pt2_Z) # Elbow
    Pt3  = c(dat_pt_curr$pt3_X, dat_pt_curr$pt3_Y, dat_pt_curr$pt3_Z) # Wrist
    Pt4  = c(dat_pt_curr$pt4_X, dat_pt_curr$pt4_Y, dat_pt_curr$pt4_Z) # End of carpometacarpus

    Pt8  = c(dat_pt_curr$pt8_X, dat_pt_curr$pt8_Y, dat_pt_curr$pt8_Z)    # Tip of most distal primary
    Pt9  = c(dat_pt_curr$pt9_X, dat_pt_curr$pt9_Y, dat_pt_curr$pt9_Z)    # Tip of last primary to model as if on the end of the carpometacarpus
    Pt10 = c(dat_pt_curr$pt10_X, dat_pt_curr$pt10_Y, dat_pt_curr$pt10_Z) # S1
    Pt11 = c(dat_pt_curr$pt11_X, dat_pt_curr$pt11_Y, dat_pt_curr$pt11_Z) # Wing root trailing edge
    Pt12 = c(dat_pt_curr$pt12_X, dat_pt_curr$pt12_Y, dat_pt_curr$pt12_Z) # Wing root leading edge
    clean_pts = rbind(Pt1,Pt2,Pt3,Pt4,Pt8,Pt9,Pt10,Pt11,Pt12)

        # solve the data
    curr_wing_data      = massprop_birdwing(dat_id_curr, dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat, clean_pts)
    curr_torsotail_data = massprop_restbody(dat_id_curr, dat_bird_curr)

    # Compute the full bird results
    fullbird = list()
    fullbird$I = matrix(0, nrow = 3, ncol = 3)
    fullbird$CG = matrix(0, nrow = 3, ncol = 1)
    # --- Mass ---
    fullbird$m = sum(subset(curr_torsotail_data, object == "m")$value,2*subset(curr_wing_data, object == "m" & component == "wing")$value)
    # --- Moment of Inertia tensor ---
    fullbird$I[1,1] = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Ixx")]) + 2*curr_wing_data$value[which(curr_wing_data$object == "Ixx" & curr_wing_data$component == "wing")]
    fullbird$I[2,2] = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Iyy")]) + 2*curr_wing_data$value[which(curr_wing_data$object == "Iyy" & curr_wing_data$component == "wing")]
    fullbird$I[3,3] = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Izz")]) + 2*curr_wing_data$value[which(curr_wing_data$object == "Izz" & curr_wing_data$component == "wing")]
    fullbird$I[1,3] = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Ixz")]) + 2*curr_wing_data$value[which(curr_wing_data$object == "Ixz" & curr_wing_data$component == "wing")]
    fullbird$I[3,1] = sum(curr_torsotail_data$value[which(curr_torsotail_data$object == "Ixz")]) + 2*curr_wing_data$value[which(curr_wing_data$object == "Ixz" & curr_wing_data$component == "wing")]

    fullbird$CG[1] = (subset(curr_torsotail_data, object == "CGx" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                      subset(curr_torsotail_data, object == "CGx" & component == "neck")$value*subset(curr_torsotail_data, object == "m" & component == "neck")$value +
                      subset(curr_torsotail_data, object == "CGx" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                      subset(curr_torsotail_data, object == "CGx" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                      2*subset(curr_wing_data, object == "CGx" & component == "wing")$value*subset(curr_torsotail_data, object == "m" & component == "wing")$value)/fullbird$m
    fullbird$CG[2] = (subset(curr_torsotail_data, object == "CGy" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                        subset(curr_torsotail_data, object == "CGy" & component == "neck")$value*subset(curr_torsotail_data, object == "m" & component == "neck")$value +
                        subset(curr_torsotail_data, object == "CGy" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                        subset(curr_torsotail_data, object == "CGy" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                        2*subset(curr_wing_data, object == "CGy" & component == "wing")$value*subset(curr_torsotail_data, object == "m" & component == "wing")$value)/fullbird$m
    fullbird$CG[3] = (subset(curr_torsotail_data, object == "CGz" & component == "head")$value*subset(curr_torsotail_data, object == "m" & component == "head")$value +
                        subset(curr_torsotail_data, object == "CGz" & component == "neck")$value*subset(curr_torsotail_data, object == "m" & component == "neck")$value +
                        subset(curr_torsotail_data, object == "CGz" & component == "torso")$value*subset(curr_torsotail_data, object == "m" & component == "torso")$value +
                        subset(curr_torsotail_data, object == "CGz" & component == "tail")$value*subset(curr_torsotail_data, object == "m" & component == "tail")$value +
                        2*subset(curr_wing_data, object == "CGz" & component == "wing")$value*subset(curr_torsotail_data, object == "m" & component == "wing")$value)/fullbird$m

    # save the full data
    mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
    column_names = c("species","BirdID","TestID","FrameID","prop_type","component","value")
    colnames(mass_properties) = column_names
    curr_full_bird = store_data(dat_id_curr,fullbird,mass_properties,"full")

    all_data = rbind(all_data, curr_wing_data, curr_torsotail_data, curr_full_bird)

  }
}

