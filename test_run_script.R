# This script is written to create the test cases

# --------------------- Read in data -----------------------
alldat   = read.csv('2020_05_25_OrientedWings.csv')
dat_bird = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'FullBird')
dat_bone = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'Bones')
dat_feat = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'Feathers')
dat_mat  = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'MaterialSpecs')

no_species = unique(dat_bird$species)

for (species in 1:no_species){
  # ----------- Filter the data to the current species ---------
  species_curr  = dat_bird$species[species]
  alldat_curr   = subset(alldat, species == species_curr)
  dat_bird_curr = subset(dat_bird, species == species_curr)
  dat_bone_curr = subset(dat_bone, species == species_curr)
  dat_feat_curr = subset(dat_feat, species == species_curr)

  w_sk = 0.045*(dat_bird_curr$total_bird_mass^0.37); # predicted distance between humerus heads (Nudds & Rayner)

  # --------------- Clean up incoming data --------------------
  # CAUTION:   The read in data has the origin at the RH humeral head.
  #            It is critical that the adj_data is given in the structural
  #            frame of reference with the origin at the vehicle reference
  #            point (VRP). VRP is assumed to be in the center of the body
  #            on the y axis between the two humerus bones

  dat_pt = alldat_curr[,8:43]

  Pt1_org = c(alldat_curr$Pt1X,alldat_curr$Pt1Y,alldat_curr$Pt1Z)

  for (i in 1:length(colnames(dat_pt))){
    if (grepl("X",colnames(alldat_curr)[i],fixed=TRUE)){
      dat_pt[,i] - Pt1_org[1]
    }
    if (grepl("Y",colnames(alldat_curr)[i],fixed=TRUE)){
      dat_pt[,i] - Pt1_org[2] + w_sk/2
    }
    if (grepl("Z",colnames(alldat_curr)[i],fixed=TRUE)){
      dat_pt[,i] - Pt1_org[3]
    }
  }

  # -------------- Begin iteration through the wings of this species ------------------------------

  for (ind_wing in 1:length(rownames(dat_pt))){
  # --------------- Bone Data ------------------------
    dat_bone_hum = subset(dat_bone_curr, bone == "Humerus")
    dat_bone_uln = subset(dat_bone_curr, bone == "Ulna")
    dat_bone_rad = subset(dat_bone_curr, bone == "Radius")
    dat_bone_car = subset(dat_bone_curr, bone == "Carpometacarpus")
    # Humerus
    massprop_bones(dat_bone_hum$bone_mass,dat_bone_hum$bone_len,dat_bone_hum$bone_out_rad,dat_bone_hum$bone_in_rad,
                   c(dat_pt$Pt1X[ind_wing], dat_pt$Pt1Y[ind_wing], dat_pt$Pt1Z[ind_wing]),
                   c(dat_pt$Pt2X[ind_wing], dat_pt$Pt2Y[ind_wing], dat_pt$Pt2Z[ind_wing]))
  }





}


