# -------------------- Store Data  -------------------------------
#' Function to store moment of inertia tensor and center of gravity vector components in long format
#'
#' @param species_curr current species code
#' @param alldat_row current dataframe containing the necessary identifying information for the wing shape
#' @param dat_mass new MOI and CG data to add to mass_properties as new rows
#' @param mass_properties current master data fram that saves the data
#' @param name name of the component
#' @return
#' @export
#'
#' @examples
store_data <- function(species_curr,alldat_row,dat_mass,mass_properties,name){
  prop_type_list = c("Ixx","Iyy","Izz","Ixy","Iyz","Ixz","CGx","CGy","CGz")
  prop_type_ind  = cbind(c(1,2,3,1,2,1),c(1,2,3,2,3,3))

  species = species_curr
  wingID  = as.character(alldat_row$WingID)
  testID  = as.character(alldat_row$TestID)
  frameID = as.character(alldat_row$frameID)

  # MOI
  for (i in 1:6){
    new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = prop_type_list[i], value = dat_mass$I[prop_type_ind[i,1],prop_type_ind[i,2]]) # saves the name and valueof the tensor component
    mass_properties = rbind(mass_properties,new_row)
  }

  #CG
  for (i in 1:3){
    new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = prop_type_list[6+i], value = dat_mass$CG[i]) # saves the name and value of the CG component
    mass_properties = rbind(mass_properties,new_row)
  }

  return(mass_properties)
}

# -------------------- Mass Properties - Halfspan bird wing  -------------------------------
#' Function that reads in anatomical data and returns the moment of inertia tensor and center
#' of gravity of a wing one side of the bird
#'
#' @param dat_bird_curr
#' @param dat_bone_curr
#' @param dat_feat_curr
#' @param dat_mat_curr
#' @param clean_pts data frame of the key positions of the bird as follows:
#' pt1x, pt1y, pt1z - Point that defines the shoulder joint
#' pt2x, pt1y, pt2z - Point that defines the elbow joint
#' pt3x, pt3y, pt3z - Point that defines the wrist joint
#' pt4x, pt4y, pt4z - Point that defines the end of carpometacarpus
#'
#' @return
#' @export
#'
#' @examples
massprop_birdwing <- function(dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat_curr, clean_pts){


  column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")
  # --------------------- Initialize variables -----------------------
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  colnames(mass_properties) = column_names
  mass_properties_bone      = mass_properties                    # specific bone data
  mass_properties_muscle    = mass_properties                    # specific muscle data
  mass_properties_feathers  = mass_properties                    # individual feather mass property data

  w_sk = 0.045*(dat_bird_curr$total_bird_mass^0.37); # predicted distance between humerus heads (Nudds & Rayner)

  # define incoming points
  Pt1 = clean_pts[1,] # shoulder
  Pt2 = clean_pts[2,] # elbow
  Pt3 = clean_pts[3,] # wrist
  Pt4 = clean_pts[4,] # end of carpometacarpus

  Pt8  = clean_pts[5,] # tip of most distal primary
  Pt9  = clean_pts[6,] # tip of last primary to model as if on the end of the carpometacarpus
  Pt10 = clean_pts[7,] # tip of S1
  Pt11 = clean_pts[8,] # tip of last secondary feather at wing root

  # --------------------------------------------------
  # --------------- Bone Data ------------------------
  # --------------------------------------------------

  rho_bone     = dat_mat$density[which(dat_mat$material == "Bone")]
  dat_bone_hum = subset(dat_bone_curr, bone == "Humerus")
  dat_bone_uln = subset(dat_bone_curr, bone == "Ulna")
  dat_bone_rad = subset(dat_bone_curr, bone == "Radius")
  dat_bone_car = subset(dat_bone_curr, bone == "Carpometacarpus")

  hum       = massprop_bones(dat_bone_hum$bone_mass,dat_bone_hum$bone_len,dat_bone_hum$bone_out_rad,dat_bone_hum$bone_in_rad,rho_bone, Pt1, Pt2)
  ulna      = massprop_bones(dat_bone_uln$bone_mass,dat_bone_uln$bone_len,dat_bone_uln$bone_out_rad,dat_bone_uln$bone_in_rad,rho_bone, Pt2, Pt3)
  radius    = massprop_bones(dat_bone_rad$bone_mass,dat_bone_rad$bone_len,dat_bone_rad$bone_out_rad,dat_bone_rad$bone_in_rad,rho_bone, Pt2, Pt3)
  car       = massprop_bones(dat_bone_car$bone_mass,dat_bone_car$bone_len,dat_bone_car$bone_out_rad,dat_bone_car$bone_in_rad,rho_bone, Pt3, Pt4)
  wristbone = massprop_pm((subset(dat_bone_curr, bone == "Ulnare")$bone_mass + subset(dat_bone_curr, bone == "Radiale")$bone_mass), Pt3)

  # --- All Bones ---
  prop_bone = list()
  # simply addition as long as about the same origin in the same frame of reference (Frame of reference: VRP | Origin: VRP)
  prop_bone$I  = hum$I + ulna$I + radius$I + car$I + wristbone$I
  # weighted average of the individual center of mass (Frame of reference: VRP | Origin: VRP)
  prop_bone$CG = (dat_bone_hum$bone_mass*hum$CG + dat_bone_uln$bone_mass*ulna$CG + dat_bone_rad$bone_mass*radius$CG + dat_bone_car$bone_mass*car$CG + (subset(dat_bone_curr, bone == "Ulnare")$bone_mass + subset(dat_bone_curr, bone == "Radiale")$bone_mass)*wristbone$CG)/sum(dat_bone_curr$bone_mass)

  # ----------------------------------------------------
  # --------------- Muscle Data ------------------------
  # ----------------------------------------------------
  rho_muscle      = dat_mat$density[which(dat_mat$material == "Muscle")]
  mass_muscles    = c()
  mass_muscles[1] = dat_bird_curr$brachial_muscle_mass
  mass_muscles[2] = dat_bird_curr$antebrachial_muscle_mass
  mass_muscles[3] = dat_bird_curr$manus_muscle_mass

  brach  = massprop_muscles(mass_muscles[1],rho_muscle,Pt1,Pt2)
  abrach = massprop_muscles(mass_muscles[2],rho_muscle,Pt2,Pt3)
  manus  = massprop_muscles(mass_muscles[3],rho_muscle,Pt3,Pt4)
  # --- All Muscles ---
  prop_muscles = list()
  # simply addition as long as about the same origin in the same frame of reference (Frame of reference: VRP | Origin: VRP)
  prop_muscles$I  = brach$I + abrach$I + manus$I
  # weighted average of the individual center of mass (Frame of reference: VRP | Origin: VRP)
  prop_muscles$CG = (mass_muscles[1]*brach$CG + mass_muscles[2]*abrach$CG + mass_muscles[3]*manus$CG)/sum(mass_muscles)

  # ---------------------------------------------------------
  # --------------- Skin/Covert Data ------------------------
  # ---------------------------------------------------------
  rho_skin     = dat_mat$density[which(dat_mat$material == "Skin")]
  mass_skin    = c()
  mass_skin[1] = dat_bird_curr$propatagiale_skin_mass
  mass_skin[2] = dat_bird_curr$manus_skin_mass
  skin_prop    = massprop_skin(mass_skin[1],rho_skin,rbind(Pt1,Pt3,Pt2))
  skin_man     = massprop_skin(mass_skin[2],rho_skin,rbind(Pt3,Pt4,Pt2))
  # --- All Skin ---
  prop_skin = list()
  # simply addition as long as about the same origin in the same frame of reference (Frame of reference: VRP | Origin: VRP)
  prop_skin$I  = skin_prop$I + skin_man$I
  # weighted average of the individual center of mass (Frame of reference: VRP | Origin: VRP)
  prop_skin$CG = (mass_skin[1]*skin_prop$CG + mass_skin[2]*skin_man$CG)/sum(mass_skin)

  # -------------------------------------------------------
  # ----------------- Feather Data ------------------------
  # -------------------------------------------------------

  #pre-define storage matrices
  res_pri    = list()
  res_pri$I  = array(dim = c(3,3,no_pri))
  res_pri$CG = array(dim = c(no_pri,3))
  res_sec    = list()
  res_sec$I  = array(dim = c(3,3,no_sec))
  res_sec$CG = array(dim = c(no_sec,3))

  # density information
  rho_cor = dat_mat$density[which(dat_mat$material == "Cortex")]
  rho_med = dat_mat$density[which(dat_mat$material == "Medullary")]

  # separate out primaries and secondaries
  primaries   = dat_feat_curr[grep("P",dat_feat_curr$Feather),]
  secondaries = dat_feat_curr[grep("S",dat_feat_curr$Feather),]
  no_sec = length(secondaries$Feather)
  no_pri = length(primaries$Feather)

  # determine the orientation and normal of each feather
  feather_info = orient_feather(primaries,secondaries,no_pri,no_sec,Pt1,Pt2,Pt3,Pt4,Pt9,Pt10,Pt11)
  # --- Primaries ---
  for (i in 1:no_pri){
    feather_name = paste("P",i,sep = "")
    info = subset(dat_feat_curr,Feather == feather_name)

    tmp = massprop_feathers(info$m_f,info$l_cal,info$l_vane, info$w_r,
                      dat_bird_curr$barb_radius, dat_bird_curr$barb_distance,
                      rho_cor,rho_med,
                      info$w_vp,info$w_vd,info$vane_angle,
                      feather_info$normal[i,],feather_info$loc_start[i,],feather_info$loc_end[i,])

    res_pri$I[,,i] = tmp$I
    res_pri$CG[i,] = tmp$CG
  }
  # --- Secondaries ---
  for (i in 1:no_sec){
    feather_name = paste("S",i,sep = "")
    info = subset(dat_feat_curr,Feather == feather_name)

    tmp = massprop_feathers(info$m_f,info$l_cal,info$l_vane, info$w_r,
                            dat_bird_curr$barb_radius, dat_bird_curr$barb_distance,
                            rho_cor,rho_med,
                            info$w_vp,info$w_vd,info$vane_angle,
                            feather_info$normal[i+no_pri,],feather_info$loc_start[i+no_pri,],feather_info$loc_end[i+no_pri,])

    res_sec$I[,,i] = tmp$I
    res_sec$CG[i,] = tmp$CG
  }

  # ----------------------------------------------------
  # ----------------- Save Data ------------------------
  # ----------------------------------------------------

  # add data to bone specific data frame
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],hum,mass_properties_bone,count,"humerus")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],ulna,mass_properties_bone,count,"ulna")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],radius,mass_properties_bone,count,"radius")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],car,mass_properties_bone,count,"carpo")
  mass_properties_bone = store_data(species_curr,alldat_curr[ind_wing,],ulnare,mass_properties_bone,count,"wristbones")

  # save all combined data groups to the master list
  mass_properties      = store_data(species_curr,alldat_curr[ind_wing,],prop_bone,mass_properties,count,"bones")
  mass_properties      = store_data(species_curr,alldat_curr[ind_wing,],prop_muscles,mass_properties,count,"muscles")


  return(mass_properties)

}
