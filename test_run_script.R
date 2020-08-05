# This script is written to create the test cases

# --------------------- Read in data -----------------------
alldat   = read.csv('2020_05_25_OrientedWings.csv')
dat_bird = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'FullBird')
dat_bone = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'Bones')
dat_feat = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'Feathers')
dat_mat  = readxl::read_xlsx('Master_AnatomicalData.xlsx', sheet = 'MaterialSpecs')

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

    test = massprop_birdwing(dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat, dat_pt_curr)


  }
}





