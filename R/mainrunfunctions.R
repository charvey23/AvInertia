

store_data <- function(species_curr,alldat_row,temp,mass_properties,count,name){
  prop_type_list = c("Ixx","Iyy","Izz","Ixy","Iyz","Ixz","CGx","CGy","CGz")
  prop_type_ind  = cbind(c(1,2,3,1,2,1),c(1,2,3,2,3,3))

  species = species_curr
  wingID  = as.character(alldat_row$WingID)
  testID  = as.character(alldat_row$TestID)
  frameID = as.character(alldat_row$frameID)

  # MOI
  for (i in 1:6){
    new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = prop_type_list[i], value = temp$I[prop_type_ind[i,1],prop_type_ind[i,2]]) # saves the name and valueof the tensor component
    mass_properties = rbind(mass_properties,new_row)
  }

  #CG
  for (i in 1:3){
    new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = prop_type_list[6+i], value = temp$CG[i]) # saves the name and value of the CG component
    mass_properties = rbind(mass_properties,new_row)
  }

  return(mass_properties)
}


massprop_birdwing <- function(dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat_curr, dat_pt_curr){


  column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")
  # --------------------- Initialize variables -----------------------
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  colnames(mass_properties) = column_names
  mass_properties_bone      = mass_properties                    # specific bone data
  mass_properties_muscle    = mass_properties                    # specific muscle data
  mass_properties_feathers  = mass_properties                    # individual feather mass property data

  w_sk = 0.045*(dat_bird_curr$total_bird_mass^0.37); # predicted distance between humerus heads (Nudds & Rayner)

  # --------------------------------------------------
  # --------------- Bone Data ------------------------
  # --------------------------------------------------

  rho_bone     = dat_mat$density[which(dat_mat$material == "Bone")]
  dat_bone_hum = subset(dat_bone_curr, bone == "Humerus")
  dat_bone_uln = subset(dat_bone_curr, bone == "Ulna")
  dat_bone_rad = subset(dat_bone_curr, bone == "Radius")
  dat_bone_car = subset(dat_bone_curr, bone == "Carpometacarpus")

  # --- Humerus ---
  hum = massprop_bones(dat_bone_hum$bone_mass,dat_bone_hum$bone_len,dat_bone_hum$bone_out_rad,dat_bone_hum$bone_in_rad,rho_bone,
                       c(dat_pt_curr$Pt1X, dat_pt_curr$Pt1Y, dat_pt_curr$Pt1Z),
                       c(dat_pt_curr$Pt2X, dat_pt_curr$Pt2Y, dat_pt_curr$Pt2Z))
  # --- Ulna ---
  ulna = massprop_bones(dat_bone_uln$bone_mass,dat_bone_uln$bone_len,dat_bone_uln$bone_out_rad,dat_bone_uln$bone_in_rad,rho_bone,
                        c(dat_pt_curr$Pt2X, dat_pt_curr$Pt2Y, dat_pt_curr$Pt2Z),
                        c(dat_pt_curr$Pt3X, dat_pt_curr$Pt3Y, dat_pt_curr$Pt3Z))
  # --- Radius ---
  radius = massprop_bones(dat_bone_rad$bone_mass,dat_bone_rad$bone_len,dat_bone_rad$bone_out_rad,dat_bone_rad$bone_in_rad,rho_bone,
                          c(dat_pt_curr$Pt2X, dat_pt_curr$Pt2Y, dat_pt_curr$Pt2Z),
                          c(dat_pt_curr$Pt3X, dat_pt_curr$Pt3Y, dat_pt_curr$Pt3Z))
  # --- Carpometacarpus ---
  car = massprop_bones(dat_bone_car$bone_mass,dat_bone_car$bone_len,dat_bone_car$bone_out_rad,dat_bone_car$bone_in_rad,rho_bone,
                       c(dat_pt_curr$Pt3X, dat_pt_curr$Pt3Y, dat_pt_curr$Pt3Z),
                       c(dat_pt_curr$Pt4X, dat_pt_curr$Pt4Y, dat_pt_curr$Pt4Z))
  # --- Wrist Bones - Ulnare/Radiale---
  wristbone = massprop_pm((subset(dat_bone_curr, bone == "Ulnare")$bone_mass + subset(dat_bone_curr, bone == "Radiale")$bone_mass),
                          c(dat_pt_curr$Pt3X, dat_pt_curr$Pt3Y, dat_pt_curr$Pt3Z))
  # --- All Bones ---
  prop_bone = list()
  # simply addition as long as about the same origin in the same frame of reference (Frame of reference: VRP | Origin: VRP)
  prop_bone$I  = hum$I + ulna$I + radius$I + car$I + wristbone$I
  # weighted average of the individual center of mass (Frame of reference: VRP | Origin: VRP)
  prop_bone$CG = (dat_bone_hum$bone_mass*hum$CG + dat_bone_uln$bone_mass*ulna$CG + dat_bone_rad$bone_mass*radius$CG + dat_bone_car$bone_mass*car$CG + (subset(dat_bone_curr, bone == "Ulnare")$bone_mass + subset(dat_bone_curr, bone == "Radiale")$bone_mass)*wristbone$CG)/sum(dat_bone_curr$bone_mass)

  # --- Save Data ---
  # add data to bone specific data frame
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],hum,mass_properties_bone,count,"humerus")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],ulna,mass_properties_bone,count,"ulna")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],radius,mass_properties_bone,count,"radius")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],car,mass_properties_bone,count,"carpo")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],ulnare,mass_properties_bone,count,"wristbones")
  # save all bone data to the master list
  mass_properties      = store_data(species_curr,alldat_curr[ind_wing,],prop_bone,mass_properties,count,"bones")

  # ----------------------------------------------------
  # --------------- Muscle Data ------------------------
  # ----------------------------------------------------





  return(mass_properties)

}
