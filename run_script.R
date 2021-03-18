devtools::load_all()

## ---------------------------------------------------------
# --------------------- Read in data -----------------------
## ---------------------------------------------------------

# CAUTION: All incoming measurements must be in SI units; adjust as required
# UPDATE REQUIRED: Should probably move final run files into the bird moment folder
path_data_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/rundata/"
path_dataout_folder = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/outputdata/"

no_pts = 11 # number of points that were digitized

# identification info on each of the individual specimens
dat_ind      = read.csv(file = paste(path_data_folder,"IDfile.csv",sep= ""))

# all of the non-wing based measurements for all specimens
dat_bird     = readxl::read_xlsx(paste(path_data_folder,"bird_measurements_readytorun.xlsx",sep= ""), sheet = 'Major body parts')

# all of the wing based measurements for all specimens
dat_wingspec = readxl::read_xlsx(paste(path_data_folder,"bird_measurements_readytorun.xlsx",sep= ""), sheet = 'Within wings')

# all of the feather measurements for one specimen of each species
dat_feat     = readxl::read_xlsx(paste(path_data_folder,"feathertrackingsheet.xlsx",sep= ""), sheet = 'FeatherData')

# all of the barb measurements for each species baed on measurements from previous studies
dat_barb     = readxl::read_xlsx(paste(path_data_folder,"feathertrackingsheet.xlsx",sep= ""), sheet = 'BarbInfo')

## ----------------------------------------------------
## ------- define all assumed material properties -----
## ----------------------------------------------------

dat_mat = list()
dat_mat$material  =	c("Muscle", "Bone", "Skin", "Cortex", "Medullary")
# SI units density: kg/m^3 - for references see Simplifications document
dat_mat$density   = c(1100,2060,1060,1150,80)

## ----------------------------------------------------
# ----------------- Clean Data ------------------------
## ----------------------------------------------------

# -- Clean up the data set relating to full bird measurements of each individual specimen--
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
# -- Clean up the data set relating to the wing specific measurements for each individual specimen --
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

# Merge the wing specific measurements with the full bird measurements
dat_bird = merge(dat_bird,dat_wingspec, by = c("species","BirdID"))

# Correct all the units in the data frame to SI - and ensure that all inputs are rounded to 0.0001g and 0.1mm.
# The small mass size is necessary for the tiny bird feathers.
for (i in 1:length(colnames(dat_bird))){
  # Adjust all masses to be in kg
  if (grepl("mass_g",colnames(dat_bird)[i],fixed=TRUE)){
    dat_bird[,i] = round(dat_bird[,i]/1000,7)
  }
  # Adjust all lengths to be in m
  if (grepl("_cm",colnames(dat_bird)[i],fixed=TRUE)){
    dat_bird[,i] = round(dat_bird[,i]/100,4)
  }
  # Adjust all lengths to be in m
  if (grepl("_mm",colnames(dat_bird)[i],fixed=TRUE)){
    dat_bird[,i] = round(dat_bird[,i]/1000,4)
  }
}

# Rename columns as necessary
names(dat_bird)[names(dat_bird) == "humerus_muscles_mass_g"] = "brachial_muscle_mass"
names(dat_bird)[names(dat_bird) == "raduln_muscles_mass_g"]  = "antebrachial_muscle_mass"
names(dat_bird)[names(dat_bird) == "cmc_muscles_mass_g"]     = "manus_muscle_mass"
names(dat_bird)[names(dat_bird) == "whole_body_mass_g"]      = "total_bird_mass"
names(dat_bird)[names(dat_bird) == "tert_feathers_mass_g"]   = "tertiary_mass"
names(dat_bird)[names(dat_bird) == "skin_coverts_mass_g"]    = "all_skin_coverts_mass"

names(dat_bird)[names(dat_bird) == "head_length_cm"]             = "head_length"
names(dat_bird)[names(dat_bird) == "head_mass_g"]                = "head_mass"
names(dat_bird)[names(dat_bird) == "head_height_max_cm"]         = "head_height"
names(dat_bird)[names(dat_bird) == "neck_mass_g"]                = "neck_mass"
names(dat_bird)[names(dat_bird) == "neck_width_cm"]              = "neck_width"
names(dat_bird)[names(dat_bird) == "neck_length_cm"]             = "neck_length"
names(dat_bird)[names(dat_bird) == "torsotail_length_cm"]        = "torsotail_length"
names(dat_bird)[names(dat_bird) == "torsotail_mass_g"]           = "torsotail_mass"
names(dat_bird)[names(dat_bird) == "tail_length_cm"]             = "tail_length"
names(dat_bird)[names(dat_bird) == "tail_width_cm"]              = "tail_width"
names(dat_bird)[names(dat_bird) == "right_leg_mass_g"]           = "right_leg_mass"
names(dat_bird)[names(dat_bird) == "left_leg_mass_"]             = "left_leg_mass"
names(dat_bird)[names(dat_bird) == "body_width_max_cm"]          = "body_width_max"
names(dat_bird)[names(dat_bird) == "body_height_max_cm"]         = "body_height_max"
names(dat_bird)[names(dat_bird) == "width_at_leg_insert_cm"]     = "body_width_at_leg_insert"
names(dat_bird)[names(dat_bird) == "x_loc_of_body_max_cm"]       = "x_loc_of_body_max"
names(dat_bird)[names(dat_bird) == "tail_width_cm"]              = "tail_width"
names(dat_bird)[names(dat_bird) == "x_loc_leg_insertion_cm"]     = "x_loc_leg_insertion"
names(dat_bird)[names(dat_bird) == "x_loc_TorsotailCoG_cm"]      = "x_loc_TorsotailCoG"
names(dat_bird)[names(dat_bird) == "z_loc_TorsotailCoG_cm"]      = "z_loc_TorsotailCoG"
names(dat_bird)[names(dat_bird) == "x_loc_of_humeral_insert_cm"] = "x_loc_humeral_insert"
names(dat_bird)[names(dat_bird) == "y_loc_of_humeral_insert_cm"] = "y_loc_humeral_insert"
names(dat_bird)[names(dat_bird) == "z_loc_of_humeral_insert_cm"] = "z_loc_humeral_insert"
names(dat_bird)[names(dat_bird) == "whole_body_mass_g"]          = "total_bird_mass"

# Correct the sign of the x measurement
dat_bird$x_loc_TorsotailCoG = - dat_bird$x_loc_TorsotailCoG

# Adjust all body z measurements to be from the clavicle position (VRP) rather than the back of the bird
dat_bird$z_loc_humeral_insert = dat_bird$z_loc_humeral_insert - dat_bird$z_dist_to_veh_ref_point_cm
dat_bird$z_loc_TorsotailCoG   = dat_bird$z_loc_TorsotailCoG - dat_bird$z_dist_to_veh_ref_point_cm

# save a singular wing mass - left if that is the one measured right otherwise
dat_bird$wing_mass = dat_bird$left_wing_mass_g
for (i in 1:nrow(dat_bird)){
  if (dat_bird$right_or_left[i] == "Right"){
    dat_bird$wing_mass[i] = dat_bird$right_wing_mass_g[i]
  }
}

## ----------------------------------------------------
## ----- Adjust the feather input properties ----------
## ----------------------------------------------------

# Update to include the barb information within the species specific info
dat_barb$barb_distance = dat_barb$barb_distance_um*10^-6
dat_barb$barb_radius   = dat_barb$barb_radius_um*10^-6
dat_bird               = merge(dat_bird,dat_barb, by=c("species"))

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
## -----------------------------------------------------
## ---------- Clean up the specimen list ---------------
## -----------------------------------------------------

specimens  = unique(dat_ind[,c("species","BirdID")])
# If there is a specimen that needs to be skipped remove it here
specimens  = specimens[-which(specimens$species == "tyt_alb" & specimens$BirdID == "18_1"),] # missing body data
specimens  = specimens[-which(specimens$species == "fal_per" & specimens$BirdID == "18_1"),] # missing body data
specimens  = specimens[-which(specimens$species == "ard_her" & specimens$BirdID == "20_317"),] # missing body data

specimens  = specimens[order(specimens$species),]
no_species = unique(specimens[,c("species")])

# --------------------- Initialize variables -----------------------
iter = 0

## ----------------------------------------------------
## --------- Iterate through each species -------------
## ----------------------------------------------------
#1:length(no_species)
for (m in 1:length(no_species)){

  # ----------- Filter the data to the current species ---------
  species_curr   = no_species[m]
  specimens_curr = subset(specimens, species == species_curr)

  # ---------- Acquire the appropriate feather data for the species --------------
  # Inputs are adjusted to be in SI units, kg and m
  # CAUTION: I only limit to the specific species - not specific birdID
  tmp = subset(dat_feat, species == species_curr)

  tmp_area   = tidyr::spread(aggregate(x = tmp$Area, by = list(tmp$Feather, tmp$Component), "mean"),"Group.2","x")
  tmp_area   = tmp_area[,c("Group.1","distal vane","proximal vane")]
  tmp_length = tidyr::spread(aggregate(x = tmp$Length, by = list(tmp$Feather, tmp$Component), "mean"),"Group.2","x")
  tmp_length = tmp_length[,c("Group.1","calamus length","vane length","rachis width")]
  tmp_angle  = tidyr::spread(aggregate(x = tmp$Angle, by = list(tmp$Feather, tmp$Component), "mean"),"Group.2","x")
  tmp_angle  = tmp_angle[,c("Group.1","calamus length","vane length")]

  dat_feat_curr = merge(tmp_area,tmp_length, by = "Group.1")

  # computes the supplement angle to the interior angle between the calamus and vane
  # assuming that the calamus is a positive angle and the vane is negative
  # negative value indicates that the feather tip points towards the body
  dat_feat_curr$vane_angle = 180 - (tmp_angle$"calamus length" - tmp_angle$"vane length")
  # Need to account for the "flipped" view of feathers from the left wing
  if(tmp$right_or_left[1] == "Left"){
    dat_feat_curr$vane_angle = -dat_feat_curr$vane_angle
  }

  # Change all units from cm to m
  dat_feat_curr$l_vane = dat_feat_curr$`vane length`*0.01
  dat_feat_curr$l_cal  = dat_feat_curr$`calamus length`*0.01
  dat_feat_curr$w_cal  = dat_feat_curr$`rachis width`*0.01
  dat_feat_curr$w_vd   = 0.01*(dat_feat_curr$"distal vane")/(dat_feat_curr$"vane length")   # compute the average width of the distal vane
  dat_feat_curr$w_vp   = 0.01*(dat_feat_curr$"proximal vane")/(dat_feat_curr$"vane length") # compute the average width of the proximal vane
  # rename last column
  names(dat_feat_curr)[names(dat_feat_curr) == "Group.1"] = "feather"
  # create the alula row
  row_alula = nrow(dat_feat_curr)+1
  dat_feat_curr[row_alula, ]       = NA
  dat_feat_curr$feather[row_alula] = "alula"
  # ---- Save the bird id to scale if these are not from the same individual ----
  dat_feat_ID = tmp$bird_id[1]
  dat_feat_species = dat_feat_curr # save this

  ## --------------------------------------------------------------------
  ## ----------- Iterate through each specimen within a species ---------
  ## --------------------------------------------------------------------

  for (k in 1:nrow(specimens_curr)){

    # --------------------- Initialize variables -----------------------
    all_data           = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
    column_names       = c("species","BirdID","TestID","FrameID","prop_type","component","value")
    colnames(all_data) = column_names
    iter = iter + 1

    # ----------- Filter the data to the current individual ---------
    birdid_curr   = specimens_curr$BirdID[k]
    # skip body if required
    if(dat_bird$use_body[which(dat_bird$BirdID == birdid_curr & dat_bird$species == species_curr)] == "N"){next}
    dat_feat_curr = dat_feat_species # ensures that this doesn't get overwritten by looping
    dat_bird_curr = subset(dat_bird, species == species_curr & BirdID == birdid_curr)

    # ----------- Read in all the wing configuration data ------------
    filename_wingconfigs = list.files(path = path_data_folder, pattern = paste(species_curr,birdid_curr,sep="_"))
    dat_wing_curr        = read.csv(paste(path_data_folder,filename_wingconfigs,sep=""))

    if(names(dat_wing_curr)[1] == "X"){
      dat_wing_curr <- dat_wing_curr[,-1]
    }

    ## --------------------------------------------------
    ## ------------------ Feather Info ------------------
    ## --------------------------------------------------

    ## ----- Update to the correct feather masses ---------
    dat_feat_curr$m_f = NA

    # save the combined alula feather masses
    dat_feat_curr$m_f[row_alula]     = dat_bird_curr[1,"alular_mass_g"]

    # save each individual feather mass
    for (i in 1:(nrow(dat_feat_curr)-1)){
      name_mass = paste(tolower(dat_feat_curr$feather[i]),"mass_g",sep = "_")
      dat_feat_curr$m_f[i] = dat_bird_curr[1,name_mass]

      # --- Scale the feather size as appropriate for this individual
      if(dat_feat_ID != specimens_curr$BirdID[k]){
        scaling_feat            = (dat_bird_curr[1,name_mass]/subset(dat_bird, species == species_curr & BirdID == dat_feat_ID)[1,name_mass])^(1/3)
        dat_feat_curr$l_vane[i] = scaling_feat*dat_feat_curr$l_vane[i]
        dat_feat_curr$l_cal[i]  = scaling_feat*dat_feat_curr$l_cal[i]
        dat_feat_curr$w_cal[i]  = scaling_feat*dat_feat_curr$w_cal[i]
        dat_feat_curr$w_vp[i]   = scaling_feat^2*dat_feat_curr$w_vp[i]
        dat_feat_curr$w_vd[i]   = scaling_feat^2*dat_feat_curr$w_vd[i]
      }
    }

    dat_feat_curr$species = species_curr

    # - Compute the I and CG of each feather - I origin is about feather CG and CG origin is start of feather, both in the feather FOR
    feather_inertia <- compute_feat_inertia(dat_mat, dat_feat_curr, dat_bird_curr)

    ## --------------------------------------------------
    ## -------------------- Bone Info -------------------
    ## --------------------------------------------------
    dat_bone_curr              = as.data.frame(matrix(nrow = 6, ncol = 5))
    names(dat_bone_curr)       = c("bone", "bone_mass", "bone_len", "bone_out_rad", "bone_in_rad")
    dat_bone_curr$bone         = c("Humerus","Ulna","Radius","Carpometacarpus","Ulnare", "Radiale")
    dat_bone_curr$bone_mass    = c(dat_bird_curr$humerus_mass_g,dat_bird_curr$ulna_mass_g,dat_bird_curr$radius_mass_g,dat_bird_curr$cmc_mass_g,
                                  dat_bird_curr$ulnare_mass_g,dat_bird_curr$radiale_mass_g)
    dat_bone_curr$bone_len     = c(dat_bird_curr$humerus_length_mm,dat_bird_curr$ulna_length_mm,dat_bird_curr$radius_length_mm,dat_bird_curr$cmc_length_mm,
                                   NA,NA)
    dat_bone_curr$bone_out_rad = c(dat_bird_curr$humerus_diameter_mm, ## ASK ABOUT THIS
                                    dat_bird_curr$ulna_diameter_mm,
                                    dat_bird_curr$radius_diameter_mm,
                                    dat_bird_curr$cmc_diameter_mm,NA,NA)
    dat_bone_curr$bone_out_rad = 0.5*dat_bone_curr$bone_out_rad
    dat_bone_curr$bone_in_rad  = 0.78*dat_bone_curr$bone_out_rad # Ref. De Margerie (2005)

    # --------------- Clean up wing 3D incoming data --------------------
    # CAUTION:   The wing data has the origin at the RH humeral head Pt1
    #            (for both left and right wings due to the aligning process).
    #            It is critical that the data is given in the structural
    #            frame of reference with the origin at the vehicle reference
    #            point (VRP). VRP is assumed to be in the center of the body
    #            on the y axis between the two humerus bones at the clavicle point

    for (i in 12:44){
      # x position - the sign has been verified - negative ensures that the shoulder is "behind" the clavicle origin
      if (grepl("X",colnames(dat_wing_curr)[i],fixed=TRUE)){
        dat_wing_curr[,i] = dat_wing_curr[,i] - (dat_bird_curr$x_loc_humeral_insert)
      }
      # y position - the sign has been verified
      # CAUTION: this assumes the wings is on the right side of the bird for all calculations
      if (grepl("Y",colnames(dat_wing_curr)[i],fixed=TRUE)){
        dat_wing_curr[,i] = dat_wing_curr[,i] + (dat_bird_curr$y_loc_humeral_insert)
      }
      # z position - the sign has been verified + positive ensures that the shoulder is "below" the clavicle origin
      if (grepl("Z",colnames(dat_wing_curr)[i],fixed=TRUE)){
        dat_wing_curr[,i] = dat_wing_curr[,i] + (dat_bird_curr$z_loc_humeral_insert)
      }
    }

    # Subset to only physically reasonable ranges
                          # Elbow can't go past the shoulder
    row_keep      = which(dat_wing_curr$pt2_Y > dat_wing_curr$pt1_Y &
                          #Pt 3, 4, 6, 7 can't go past the edge of body Pt 12.
                          dat_wing_curr$pt3_Y > dat_wing_curr$pt12_Y & dat_wing_curr$pt4_Y > dat_wing_curr$pt12_Y &
                          dat_wing_curr$pt6_Y > dat_wing_curr$pt12_Y & dat_wing_curr$pt7_Y > dat_wing_curr$pt12_Y &
                          # Pt 8 and 9 can't cross body centerline and body edge can't be more interior than shoulder
                          dat_wing_curr$pt8_Y > 0 & dat_wing_curr$pt9_Y > 0 & dat_wing_curr$pt12_Y > dat_wing_curr$pt1_Y)

    dat_wing_curr = dat_wing_curr[row_keep,]

    # Need a blank data frame to save each iteration
    mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
    column_names = c("species","BirdID","TestID","FrameID","prop_type","component","value")
    colnames(mass_properties) = column_names

    # -------- Save the final input info used for this specimen -------
    dat_id_curr = dat_wing_curr[1,c("species","BirdID","TestID","FrameID")]

    # Compute the CG and I for the body without the wings
    curr_torsotail_data = massprop_restbody(dat_id_curr, dat_bird_curr)


    ## --------------------------------------------------------------------------------
    ## ------------------ Iterate through all wing configurations ---------------------
    ## --------------------------------------------------------------------------------

    for (ind_wing in 1:nrow(dat_wing_curr)){

      # both data frames below should only be one row of input points
      dat_pt_curr = dat_wing_curr[ind_wing,]
      dat_id_curr = dat_pt_curr[,c("species","BirdID","TestID","FrameID")]
      dat_id_curr$TestID = paste(dat_pt_curr$BirdID_FrameSpec, dat_id_curr$TestID, sep = "_")

      # Initialize common pts
      Pt1  = c(dat_pt_curr$pt1_X, dat_pt_curr$pt1_Y, dat_pt_curr$pt1_Z) # Shoulder
      Pt2  = c(dat_pt_curr$pt2_X, dat_pt_curr$pt2_Y, dat_pt_curr$pt2_Z) # Elbow
      Pt3  = c(dat_pt_curr$pt3_X, dat_pt_curr$pt3_Y, dat_pt_curr$pt3_Z) # Wrist
      Pt4  = c(dat_pt_curr$pt4_X, dat_pt_curr$pt4_Y, dat_pt_curr$pt4_Z) # End of carpometacarpus

      Pt8  = c(dat_pt_curr$pt8_X, dat_pt_curr$pt8_Y, dat_pt_curr$pt8_Z)    # Tip of most distal primary
      Pt9  = c(dat_pt_curr$pt9_X, dat_pt_curr$pt9_Y, dat_pt_curr$pt9_Z)    # Tip of last primary to model as if on the end of the carpometacarpus

      # if the first secondary goes into where body would be rotate the last secondary back outwards - use the shoulder to approximate the width of the body at this point
      if(dat_pt_curr$pt10_Y < dat_pt_curr$pt1_Y){
        dat_pt_curr$pt10_X = dat_pt_curr$pt2_X - sqrt((dat_pt_curr$pt10_Y-dat_pt_curr$pt2_Y)^2+
                                                        (dat_pt_curr$pt10_X-dat_pt_curr$pt2_X)^2-
                                                        ((dat_pt_curr$pt1_Y + 0.005)-dat_pt_curr$pt2_Y)^2)
        dat_pt_curr$pt10_Y = dat_pt_curr$pt1_Y + 0.005 # leave 5mm gap to properly distribute the feathers in this range
        dat_wing_curr$pt10_Y[ind_wing] = dat_pt_curr$pt10_Y
        dat_wing_curr$pt10_X[ind_wing] = dat_pt_curr$pt10_X
      }

      Pt10 = c(dat_pt_curr$pt10_X, dat_pt_curr$pt10_Y, dat_pt_curr$pt10_Z) # S1

      # if the last secondary goes into where body would be rotate the last secondary back outwards - use the shoulder to approximate the width of the body at this point
      if(dat_pt_curr$pt11_Y < dat_pt_curr$pt1_Y){
        dat_pt_curr$pt11_X = dat_pt_curr$pt2_X - sqrt((dat_pt_curr$pt11_Y-dat_pt_curr$pt2_Y)^2+
                                                        (dat_pt_curr$pt11_X-dat_pt_curr$pt2_X)^2-
                                                        (dat_pt_curr$pt1_Y-dat_pt_curr$pt2_Y)^2)
        dat_pt_curr$pt11_Y = dat_pt_curr$pt1_Y
        dat_wing_curr$pt11_Y[ind_wing] = dat_pt_curr$pt11_Y
        dat_wing_curr$pt11_X[ind_wing] = dat_pt_curr$pt11_X
      }
      Pt11 = c(dat_pt_curr$pt11_X, dat_pt_curr$pt11_Y, dat_pt_curr$pt11_Z) # Wing root trailing edge
      Pt12 = c(dat_pt_curr$pt12_X, dat_pt_curr$pt12_Y, dat_pt_curr$pt12_Z) # Wing root leading edge
      clean_pts = rbind(Pt1,Pt2,Pt3,Pt4,Pt8,Pt9,Pt10,Pt11,Pt12)

      # Compute the CG and I for the wing configuration
      curr_wing_data      = massprop_birdwing(dat_id_curr, dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat, clean_pts, feather_inertia)
      # Combine the torso and wing outputs
      curr_full_bird      = combine_inertialprop(curr_torsotail_data,curr_wing_data,curr_wing_data, symmetric=TRUE)

      new_row1 = data.frame(species = dat_id_curr$species[1], BirdID = dat_id_curr$BirdID[1],TestID = dat_id_curr$TestID[1], FrameID = dat_id_curr$FrameID[1],
                            component = "wing", object = "S_proj", value = dat_pt_curr$S_proj)

      all_data = rbind(all_data, curr_wing_data, curr_full_bird, new_row1)
    } # end of the individual wing configuration loop
    # for the sake of memory need to recast from long to wide format to save
    all_data = reshape2::dcast(all_data, species + BirdID + TestID + FrameID ~ component + object, value.var="value")

    filename_output = paste(path_dataout_folder,format(Sys.Date(), "%Y_%m_%d"),"_",species_curr,"_",birdid_curr,"_results.csv",sep="")
    write.csv(all_data,filename_output)

    # ---- save all the individual info that does not change with wing position -------
    if (iter == 1){
      dat_wing_all = dat_wing_curr        # wing specimen info
      dat_feat_all = dat_feat_curr        # feather specimen info
      dat_bird_all = dat_bird_curr        # entire body specimen info
      dat_body_all = curr_torsotail_data  # CG and MOI from the torso, tail, head and neck for this individual
    }else{
      dat_wing_all = rbind(dat_wing_all,dat_wing_curr)
      dat_feat_all = rbind(dat_feat_all,dat_feat_curr)
      dat_bird_all = rbind(dat_bird_all,dat_bird_curr)
      dat_body_all = rbind(dat_body_all,curr_torsotail_data)
    }

    remove(all_data)
   } # end of the specimen loop
} # end of the species loop

# --------------------------------------------------
# --------------- Save combined data ---------------
# --------------------------------------------------
dat_wing_all$TestID = paste(dat_wing_all$BirdID_FrameSpec, dat_wing_all$TestID, sep = "_")
filename = paste(path_dataout_folder,format(Sys.Date(), "%Y_%m_%d"),"_allspecimen_winginfo.csv",sep="")
write.csv(dat_wing_all,filename)

filename = paste(path_dataout_folder,format(Sys.Date(), "%Y_%m_%d"),"_allspecimen_birdinfo.csv",sep="")
write.csv(dat_bird_all,filename)

filename = paste(path_dataout_folder,format(Sys.Date(), "%Y_%m_%d"),"_allspecimen_featherinfo.csv",sep="")
write.csv(dat_feat_all,filename)

filename = paste(path_dataout_folder,format(Sys.Date(), "%Y_%m_%d"),"_bodyCGandMOI_correctneck.csv",sep="")
write.csv(dat_body_all,filename)
