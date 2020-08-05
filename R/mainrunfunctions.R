
massprop_birdwing <- function(dat_bird, ){


  no_species = unique(dat_bird$species)
  column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")
  # --------------------- Initialize variables -----------------------
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  mass_properties_bone = mass_properties                         # specific bone data
  mass_properties_muscle = mass_properties                       # specific muscle data
  mass_properties_feathers = mass_properties                     # individual feather mass property data

  prop_type_list = c("Ixx","Iyy","Izz","Ixy","Iyz","Ixz","CGx","CGy","CGz")
  prop_type_ind  = cbind(c(1,2,3,1,2,1),c(1,2,3,2,3,3))

  for (species in 1:no_species){
    # ----------- Filter the data to the current species ---------
    species_curr  = dat_bird$species[species]
    alldat_curr   = subset(alldat, species == species_curr)
    dat_bird_curr = subset(dat_bird, species == species_curr)
    dat_bone_curr = subset(dat_bone, species == species_curr)
    dat_feat_curr = subset(dat_feat, species == species_curr)

    w_sk = 0.045*(dat_bird_curr$total_bird_mass^0.37); # predicted distance between humerus heads (Nudds & Rayner)

    # --------------- Clean up incoming data --------------------
    # CAUTION:   The read in data has the origin at the RH humeral head. Pt1X, Pt1Y, Pt1Z = 0
    #            It is critical that the adj_data is given in the structural
    #            frame of reference with the origin at the vehicle reference
    #            point (VRP). VRP is assumed to be in the center of the body
    #            on the y axis between the two humerus bones

    dat_pt = alldat_curr[,8:43]

    for (i in 1:length(colnames(dat_pt))){
      if (grepl("Y",colnames(dat_pt)[i],fixed=TRUE)){
        dat_pt[,i] = dat_pt[,i] + w_sk/2
      }
    }

    # -------------- Begin iteration through the wings of this species ------------------------------

    for (ind_wing in 1:10){
      count_start = count

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
                           c(dat_pt$Pt1X[ind_wing], dat_pt$Pt1Y[ind_wing], dat_pt$Pt1Z[ind_wing]),
                           c(dat_pt$Pt2X[ind_wing], dat_pt$Pt2Y[ind_wing], dat_pt$Pt2Z[ind_wing]))
      # add data to bone specific data frame
      mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],hum,mass_properties_bone,count,"humerus",prop_type_list,prop_type_ind)

      # --- Ulna ---
      ulna = massprop_bones(dat_bone_uln$bone_mass,dat_bone_uln$bone_len,dat_bone_uln$bone_out_rad,dat_bone_uln$bone_in_rad,rho_bone,
                            c(dat_pt$Pt2X[ind_wing], dat_pt$Pt2Y[ind_wing], dat_pt$Pt2Z[ind_wing]),
                            c(dat_pt$Pt3X[ind_wing], dat_pt$Pt3Y[ind_wing], dat_pt$Pt3Z[ind_wing]))
      # add data to bone specific data frame
      mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],ulna,mass_properties_bone,count,"ulna",prop_type_list,prop_type_ind)

      # --- Radius ---
      radius = massprop_bones(dat_bone_rad$bone_mass,dat_bone_rad$bone_len,dat_bone_rad$bone_out_rad,dat_bone_rad$bone_in_rad,rho_bone,
                              c(dat_pt$Pt2X[ind_wing], dat_pt$Pt2Y[ind_wing], dat_pt$Pt2Z[ind_wing]),
                              c(dat_pt$Pt3X[ind_wing], dat_pt$Pt3Y[ind_wing], dat_pt$Pt3Z[ind_wing]))
      # add data to bone specific data frame
      mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],radius,mass_properties_bone,count,"radius",prop_type_list,prop_type_ind)

      # --- Carpometacarpus ---
      car = massprop_bones(dat_bone_car$bone_mass,dat_bone_car$bone_len,dat_bone_car$bone_out_rad,dat_bone_car$bone_in_rad,rho_bone,
                           c(dat_pt$Pt3X[ind_wing], dat_pt$Pt3Y[ind_wing], dat_pt$Pt3Z[ind_wing]),
                           c(dat_pt$Pt4X[ind_wing], dat_pt$Pt4Y[ind_wing], dat_pt$Pt4Z[ind_wing]))
      # add data to bone specific data frame
      mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],car,mass_properties_bone,count,"carpo",prop_type_list,prop_type_ind)

      # --- Wrist Bones - Ulnare/Radiale---
      wristbone = massprop_pm((subset(dat_bone_curr, bone == "Ulnare")$bone_mass + subset(dat_bone_curr, bone == "Radiale")$bone_mass),
                              c(dat_pt$Pt3X[ind_wing], dat_pt$Pt3Y[ind_wing], dat_pt$Pt3Z[ind_wing]))
      # add data to bone specific data frame
      mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],ulnare,mass_properties_bone,count,"wristbones",prop_type_list,prop_type_ind)

      prop_bone = list()
      # simply addition as long as about the same origin in the same frame of reference (Frame of reference: VRP | Origin: VRP)
      prop_bone$I  = hum$I + ulna$I + radius$I + car$I + wristbone$I
      # weighted average of the individual center of mass (Frame of reference: VRP | Origin: VRP)
      prop_bone$CG = (dat_bone_hum$bone_mass*hum$CG + dat_bone_uln$bone_mass*ulna$CG + dat_bone_rad$bone_mass*radius$CG + dat_bone_car$bone_mass*car$CG + (subset(dat_bone_curr, bone == "Ulnare")$bone_mass + subset(dat_bone_curr, bone == "Radiale")$bone_mass)*wristbone$CG)/sum(dat_bone_curr$bone_mass)
      # save all bone data to the master list
      mass_properties = store_data(species_curr,alldat_curr[ind_wing,],prop_bone,mass_properties,count,"bones",prop_type_list,prop_type_ind)
    }


}
