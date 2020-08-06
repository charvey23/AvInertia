# This script is written to create the test cases

# --------------------- Read in data -----------------------
alldat   = read.csv('2020_05_25_OrientedWings.csv')
dat_bird = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'FullBird')
dat_bone = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'Bones')
dat_feat = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'Feathers')
dat_mat  = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'MaterialSpecs')

# adjust to SI units
dat_feat$l_cal  = 0.01*dat_feat$l_cal
dat_feat$l_vane = 0.01*dat_feat$l_vane
dat_feat$w_r    = 0.01*dat_feat$w_r
dat_feat$w_vp   = 0.01*dat_feat$w_vp
dat_feat$w_vd   = 0.01*dat_feat$w_vd
dat_feat$m_f    = 0.001*dat_feat$m_f

no_species = unique(dat_bird$species)

for (species in 1:length(alldat$species)){
  # ----------- Filter the data to the current species ---------
  species_curr  = dat_bird$species[species]
  alldat_curr   = subset(alldat, species == species_curr)
  dat_bird_curr = subset(dat_bird, species == species_curr)
  dat_bone_curr = subset(dat_bone, species == species_curr)
  dat_feat_curr = subset(dat_feat, species == species_curr)

  dat_pt = alldat_curr[,8:43]

  # --------------- Clean up incoming data --------------------
  # CAUTION:   The read in data has the origin at the RH humeral head. Pt1X, Pt1Y, Pt1Z = 0
  #            It is critical that the data is given in the structural
  #            frame of reference with the origin at the vehicle reference
  #            point (VRP). VRP is assumed to be in the center of the body
  #            on the y axis between the two humerus bones

  for (i in 1:length(colnames(dat_pt))){
    if (grepl("Y",colnames(dat_pt)[i],fixed=TRUE)){
      dat_pt[,i] = dat_pt[,i] + w_sk/2
    }
  }

  # -------------- Begin iteration through the wings of this species ------------------------------

  for (ind_wing in 1:10){

    dat_pt_curr = dat_pt[ind_wing,] # should only be one row of input points

    # Initialize common pts
    Pt1  = c(dat_pt_curr$Pt1X, dat_pt_curr$Pt1Y, dat_pt_curr$Pt1Z) # Shoulder
    Pt2  = c(dat_pt_curr$Pt2X, dat_pt_curr$Pt2Y, dat_pt_curr$Pt2Z) # Elbow
    Pt3  = c(dat_pt_curr$Pt3X, dat_pt_curr$Pt3Y, dat_pt_curr$Pt3Z) # Wrist
    Pt4  = c(dat_pt_curr$Pt4X, dat_pt_curr$Pt4Y, dat_pt_curr$Pt4Z) # End of carpometacarpus

    Pt8  = c(dat_pt_curr$Pt8X, dat_pt_curr$Pt8Y, dat_pt_curr$Pt8Z) # Tip of most distal primary
    Pt9  = c(dat_pt_curr$Pt9X, dat_pt_curr$Pt9Y, dat_pt_curr$Pt9Z) # Tip of last primary to model as if on the end of the carpometacarpus
    Pt10 = c(dat_pt_curr$Pt10X, dat_pt_curr$Pt10Y, dat_pt_curr$Pt10Z) # S1
    Pt11 = c(dat_pt_curr$Pt11X, dat_pt_curr$Pt11Y, dat_pt_curr$Pt11Z) # Wing root trailing edge

    clean_pts = rbind(Pt1,Pt2,Pt3,Pt4,Pt8,Pt9,Pt10,Pt11)

    feather_info = list()

    test = massprop_birdwing(dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat, dat_pt_curr,clean_pts)


  }
}





