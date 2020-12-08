# This script is written to create the test cases

devtools::load_all()

# --------------------- Read in data -----------------------
alldat   = read.csv('./data/2020_05_25_OrientedWings.csv')
dat_bird = readxl::read_xlsx('./data/Master_AnatomicalData.xlsx', sheet = 'FullBird')
dat_bone = readxl::read_xlsx('./data/Master_AnatomicalData.xlsx', sheet = 'Bones')
dat_feat = readxl::read_xlsx('./data/Master_AnatomicalData.xlsx', sheet = 'Feathers')
dat_mat  = readxl::read_xlsx('./data/Master_AnatomicalData.xlsx', sheet = 'MaterialSpecs')

# adjust to SI units
dat_feat$l_cal  = 0.01*dat_feat$l_cal
dat_feat$l_vane = 0.01*dat_feat$l_vane
dat_feat$w_cal  = 0.01*dat_feat$w_cal
dat_feat$w_vp   = 0.01*dat_feat$w_vp
dat_feat$w_vd   = 0.01*dat_feat$w_vd
dat_feat$m_f    = 0.001*dat_feat$m_f

# --------------------- Initialize variables -----------------------
all_data = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")
colnames(all_data) = column_names

# ----------- Iterate through each species ---------
no_species = unique(dat_bird$species)

for (species in 1:1){ # eventually replace with length(dat_bird$species)

  # ----------- Filter the data to the current species ---------
  species_curr  = dat_bird$species[species]
  alldat_curr   = subset(alldat, species == species_curr & sweep == 0 & dihedral == 0)
  dat_bird_curr = subset(dat_bird, species == species_curr)
  dat_bone_curr = subset(dat_bone, species == species_curr)
  dat_feat_curr = subset(dat_feat, species == species_curr)

  dat_pt = alldat_curr[,8:43]

  # --------------- Clean up incoming data --------------------
  # CAUTION:   The read in data has the origin at the RH humeral head. Pt1X, Pt1Y, Pt1Z = 0
  #            It is critical that the data is given in the structural
  #            frame of reference with the origin at the vehicle reference
  #            point (VRP). VRP is assumed to be in the center of the body
  #            on the y axis between the two humerus bones at the clavicle point

  for (i in 1:length(colnames(dat_pt))){
    # x position
    if (grepl("X",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] + dat_bird_curr$x_loc_of_humeral_insert_m
    }
    # y position - CAUTION THE WING IS ASSUMED TO BE THE RIGHT SIDE WING HENCE THE POSITIVE ADDITION
    if (grepl("Y",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] + dat_bird_curr$y_loc_of_humeral_insert_m
    }
    # z position
    if (grepl("Z",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] + dat_bird_curr$z_loc_of_humeral_insert_m
    }
  }

  # -------------- Iterate through the wings of this species ------------------------------

  for (ind_wing in 1:length(alldat_curr$frameID)){

    dat_pt_curr = dat_pt[ind_wing,] # should only be one row of input points

    # Initialize common pts
    Pt1  = c(dat_pt_curr$Pt1X, dat_pt_curr$Pt1Y, dat_pt_curr$Pt1Z) # Shoulder
    Pt2  = c(dat_pt_curr$Pt2X, dat_pt_curr$Pt2Y, dat_pt_curr$Pt2Z) # Elbow
    Pt3  = c(dat_pt_curr$Pt3X, dat_pt_curr$Pt3Y, dat_pt_curr$Pt3Z) # Wrist
    Pt4  = c(dat_pt_curr$Pt4X, dat_pt_curr$Pt4Y, dat_pt_curr$Pt4Z) # End of carpometacarpus

    Pt8  = c(dat_pt_curr$Pt8X, dat_pt_curr$Pt8Y, dat_pt_curr$Pt8Z)    # Tip of most distal primary
    Pt9  = c(dat_pt_curr$Pt9X, dat_pt_curr$Pt9Y, dat_pt_curr$Pt9Z)    # Tip of last primary to model as if on the end of the carpometacarpus
    Pt10 = c(dat_pt_curr$Pt10X, dat_pt_curr$Pt10Y, dat_pt_curr$Pt10Z) # S1
    Pt11 = c(dat_pt_curr$Pt11X, dat_pt_curr$Pt11Y, dat_pt_curr$Pt11Z) # Wing root trailing edge
    Pt12 = c(dat_pt_curr$Pt12X, dat_pt_curr$Pt12Y, dat_pt_curr$Pt12Z) # Wing root leading edge
    clean_pts = rbind(Pt1,Pt2,Pt3,Pt4,Pt8,Pt9,Pt10,Pt11,Pt12)

    # solve the data
    curr_wing_data      = massprop_birdwing(alldat_curr[ind_wing,], dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat, clean_pts)
    curr_torsotail_data = massprop_restbody(alldat_curr[ind_wing,], dat_bird_curr)

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

    fullbird$CG[1] = ((curr_torsotail_data$value[which(curr_torsotail_data$object == "CGx" & curr_torsotail_data$component == "head")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "head")]) +
                    (curr_torsotail_data$value[which(curr_torsotail_data$object == "CGx" & curr_torsotail_data$component == "neck")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "neck")]) +
                    (curr_torsotail_data$value[which(curr_torsotail_data$object == "CGx" & curr_torsotail_data$component == "torso")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "torso")]) +
                    2*(curr_wing_data$value[which(curr_wing_data$object == "CGx" & curr_wing_data$component == "wing")]*curr_wing_data$value[which(curr_wing_data$object == "m" & curr_wing_data$component == "wing")]))/fullbird$m
    fullbird$CG[2] = ((curr_torsotail_data$value[which(curr_torsotail_data$object == "CGy" & curr_torsotail_data$component == "head")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "head")]) +
                        (curr_torsotail_data$value[which(curr_torsotail_data$object == "CGy" & curr_torsotail_data$component == "neck")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "neck")]) +
                        (curr_torsotail_data$value[which(curr_torsotail_data$object == "CGy" & curr_torsotail_data$component == "torso")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "torso")]) +
                        (curr_wing_data$value[which(curr_wing_data$object == "CGy" & curr_wing_data$component == "wing")]*curr_wing_data$value[which(curr_wing_data$object == "m" & curr_wing_data$component == "wing")]) -
                        (curr_wing_data$value[which(curr_wing_data$object == "CGy" & curr_wing_data$component == "wing")]*curr_wing_data$value[which(curr_wing_data$object == "m" & curr_wing_data$component == "wing")]))/fullbird$m
    fullbird$CG[3] = ((curr_torsotail_data$value[which(curr_torsotail_data$object == "CGz" & curr_torsotail_data$component == "head")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "head")]) +
                        (curr_torsotail_data$value[which(curr_torsotail_data$object == "CGz" & curr_torsotail_data$component == "neck")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "neck")]) +
                        (curr_torsotail_data$value[which(curr_torsotail_data$object == "CGz" & curr_torsotail_data$component == "torso")]*curr_torsotail_data$value[which(curr_torsotail_data$object == "m" & curr_torsotail_data$component == "torso")]) +
                        2*(curr_wing_data$value[which(curr_wing_data$object == "CGz" & curr_wing_data$component == "wing")]*curr_wing_data$value[which(curr_wing_data$object == "m" & curr_wing_data$component == "wing")]))/fullbird$m

    # save the full data
    mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
    column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")
    colnames(mass_properties) = column_names
    curr_full_bird = store_data(alldat_curr[ind_wing,],fullbird,mass_properties,"full")

    all_data = rbind(all_data, curr_wing_data, curr_torsotail_data, curr_full_bird)

  }
}

