# Script containing the overarching function that calculates the moment of inertia and CG for an entire bird
# all written by Christina Harvey
# last updated: 2020-10-13

# -------------------- Mass Properties - Head and neck -------------------------------
massprop_restbody <- function(dat_wingID_curr, dat_bird_curr){
  # --------------------- Initialize variables -----------------------
  # pre-define storage matrices
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")

  # -------------------------------------------------------------
  # ----------------- Head and neck data ------------------------
  # -------------------------------------------------------------
  neck_start = c(0,0,0)
  neck_end   = c(dat_bird_curr$neck_length_m,0,0)
  head_end   = neck_end + c(dat_bird_curr$head_length_m,0,0)
  # Calculate the effects of the head
  head       = massprop_head(dat_bird_curr$head_mass_kg,0.5*dat_bird_curr$head_height_m,dat_bird_curr$head_length_m,neck_end,head_end)
  # Calculate the effects of the neck
  neck       = massprop_neck(dat_bird_curr$neck_mass_kg,0.5*dat_bird_curr$neck_width_m,dat_bird_curr$neck_length_m,neck_start,neck_end)
  # -------------------------------------------------------------
  # ------------------- Torso and tail data -------------------------
  # -------------------------------------------------------------
  tail_start  = c(-(dat_bird_curr$torsotail_length_m-dat_bird_curr$tail_length_m),0,0)
  tail_end    = c(-dat_bird_curr$torsotail_length_m,0,0)
  m_legs    = dat_bird_curr$right_leg_mass_kg + dat_bird_curr$left_leg_mass_kg
  torso     = massprop_torso(dat_bird_curr$torsotail_mass_kg, m_legs, dat_bird_curr$body_width_max_m, dat_bird_curr$body_height_max_m,
                                  dat_bird_curr$x_loc_of_body_max_m, dat_bird_curr$body_width_at_leg_insert_m, dat_bird_curr$x_loc_leg_insertion_m,
                                  dat_bird_curr$torsotail_length_m, dat_bird_curr$x_loc_TorsotailCoG_m, dat_bird_curr$z_loc_TorsotailCoG_m,
                                  neck_start, tail_end)

  tail   = massprop_tail(dat_bird_curr$tail_mass_kg, dat_bird_curr$tail_length_m,dat_bird_curr$tail_width_m,dat_bird_curr$torsotail_length_m,tail_start,tail_end)
  # ----------------------------------------------------
  # ----------------- Save Data ------------------------
  # ----------------------------------------------------
  # ---- Head ----
  mass_properties = store_data(dat_wingID_curr,head,mass_properties,"head")
  # ---- Neck ----
  mass_properties = store_data(dat_wingID_curr,neck,mass_properties,"neck")
  # ---- Torso/Legs ----
  mass_properties = store_data(dat_wingID_curr,torso,mass_properties,"torso")
  # ---- Tail ----
  mass_properties = store_data(dat_wingID_curr,tail,mass_properties,"tail")

  return(mass_properties)
}


# -------------------- Mass Properties - Halfspan bird wing  -------------------------------
#' Halfspan bird wing mass properties
#'
#' Function that reads in anatomical data and returns the moment of inertia tensor and center
#' of gravity of a wing one side of the bird
#'
#' @param dat_wingID_curr Dataframe related to the current bird wing ID info that must include the following columns:
#' \itemize{
#' \item{species}{Species ID code as a string}
#' \item{WingID}{Wing ID code as a string}
#' \item{TestID}{Test ID code as a string}
#' \item{frameID}{Video frame ID code as a string}
#' }
#'
#' @param dat_bird_curr Dataframe related to the current bird wing that must include the following columns:
#' \itemize{
#' \item{total_bird_mass}{Mass of full bird for the current wing (kg)}
#' \item{wing_mass}{Mass of one wing, should be the current wing (kg)}
#' \item{barb_radius}{Radius of feather barb  for current species (m)}
#' \item{barb_distance}{Distance between feather barbs for current species (m)}
#' \item{brachial_muscle_mass}{Mass of all muscles in the brachial region of the wing (kg)}
#' \item{antebrachial_muscle_mass}{Mass of all muscles in the antebrachial region of the wing (kg)}
#' \item{manus_muscle_mass}{Mass of all muscles in the manus region of the wing (kg)}
#' }
#'
#' @param dat_bone_curr Dataframe related to the current bird wing bones that must include the following columns:
#' \itemize{
#' \item{bone}{Bone ID code. Must include: "Humerus","Ulna","Radius","Carpometacarpus,"Ulnare" and "Radiale".}
#' \item{bone_mass}{Mass of bone in the same row as the appropriate bone ID code (kg)}
#' \item{bone_len}{Length of bone in the same row as the appropriate bone ID code (m)}
#' \item{bone_out_rad}{Outer radius of bone in the same row as the appropriate bone ID code (m)}
#' \item{bone_in_rad}{Inner radius of bone in the same row as the appropriate bone ID code (m)}
#' }
#'
#' @param dat_feat_curr Dataframe related to the current bird wing feathers input as a dataframe with the following structure:
#' \itemize{
#' \item{feather}{Feather ID code. Must be in standard format i.e. 1st primary is "P1", third secondary is "S3", etc.
#' Alula feathers should be grouped and named "alula".}
#' \item{m_f}{Mass of feather in the same row as the appropriate feather ID code (kg)}
#' \item{l_cal}{Length of calamus in the same row as the appropriate feather ID code (m)}
#' \item{l_vane}{Length of rachis/vane in the same row as the appropriate feather ID code (m)}
#' \item{w_cal}{Width (diameter) of calamus in the same row as the appropriate feather ID code (m)}
#' \item{w_vp}{Width of proximal vane (average value) in the same row as the appropriate feather ID code (m)}
#' \item{w_vd}{Width of distal vane (average value)  in the same row as the appropriate feather ID code (m)}
#' }
#' NOTE: Alula feathers will be treated as point mass so only the mass of the feathers is required. Other columns can be left blank.
#'
#' @param dat_mat_curr Dataframe related to the current species input as a dataframe with the following structure:
#' \itemize{
#' \item{material}{Material information. Must include the following: "Bone","Skin","Muscle","Cortex", "Medullary"}
#' \item{density}{Density of each material (kg/m^3)}
#' }
#'
#' @param clean_pts Dataframe of the key positions of the bird as follows:
#' \itemize{
#' \item{pt1x, pt1y, pt1z}{Point on the shoulder joint}
#' \item{pt2x, pt1y, pt2z}{Point on the elbow joint}
#' \item{pt3x, pt3y, pt3z}{Point on the wrist joint}
#' \item{pt4x, pt4y, pt4z}{Point on the end of carpometacarpus}
#' \item{pt8x, pt8y, pt8z}{Point on tip of most distal primary}
#' \item{pt9x, pt9y, pt9z}{Point that defines the end of carpometacarpus}
#' \item{pt10x, pt10y, pt10z}{Point on tip of last primary to model as if on the end of the carpometacarpus}
#' \item{pt11x, pt11y, pt11z}{Point on tip of most proximal feather (wing root trailing edge)}
#' }
#'
#' @section CAUTION:
#'          All points must all have the vehicle reference point (VRP) as their origin
#'          and the vehicle major axes as their frame of reference. This is normally selected so that
#'          the VRP is in line with the body center of gravity. Ensure the axes used represent a right-handed axis system.
#'
#' @author Christina Harvey
#'
#' @return Function returns a dataframe that includes the moment of inertia and center of gravity of
#' one wing about the VRP in the VRP frame and that of each major anatomical group i.e. skin, feathers, bones, muscles.
#'
#'
#'
#'
#'
#' @export
#'
massprop_birdwing <- function(dat_wingID_curr, dat_bird_curr, dat_bone_curr, dat_feat_curr, dat_mat_curr, clean_pts){

  # --------------------- Initialize variables -----------------------
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  column_names = c("species","WingID","TestID","FrameID","prop_type","component","value")
  colnames(mass_properties) = column_names
  mass_properties_bone      = mass_properties                    # specific bone data
  mass_properties_muscle    = mass_properties                    # specific muscle data
  mass_properties_skin      = mass_properties                    # specific skin data
  mass_properties_feathers  = mass_properties                    # individual feather mass property data

  # define incoming points
  Pt1 = clean_pts[1,] # shoulder
  Pt2 = clean_pts[2,] # elbow
  Pt3 = clean_pts[3,] # wrist
  Pt4 = clean_pts[4,] # end of carpometacarpus

  Pt8  = clean_pts[5,] # tip of most distal primary
  Pt9  = clean_pts[6,] # tip of last primary to model as if on the end of the carpometacarpus
  Pt10 = clean_pts[7,] # tip of S1
  Pt11 = clean_pts[8,] # tip of last secondary feather at wing root
  Pt12 = clean_pts[9,] # leading edge where wing was removed from the body on the wing root
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
  prop_bone$m  = sum(dat_bone_curr$bone_mass)
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
  prop_muscles$m  = sum(mass_muscles)

  # ---------------------------------------------------------
  # --------------- Skin/Covert Data ------------------------
  # ---------------------------------------------------------

  rho_skin   = dat_mat$density[which(dat_mat$material == "Skin")]
  mass_skin  = dat_bird_curr$all_skin_coverts_mass - dat_bird_curr$tertiary_mass
  prop_skin  = massprop_skin(mass_skin,rho_skin,rbind(Pt1,Pt3,Pt2))

  # -------------------------------------------------------
  # ----------------- Feather Data ------------------------
  # -------------------------------------------------------

  # density information
  rho_cor = dat_mat$density[which(dat_mat$material == "Cortex")]
  rho_med = dat_mat$density[which(dat_mat$material == "Medullary")]

  # separate out primaries and secondaries
  primaries   = dat_feat_curr[grep("P",dat_feat_curr$feather),]
  secondaries = dat_feat_curr[grep("S",dat_feat_curr$feather),]
  no_sec = length(secondaries$feather)
  no_pri = length(primaries$feather)

  #pre-define storage matrices
  res_pri     = list()
  res_pri$I   = array(dim = c(3,3,no_pri))
  res_pri$CG  = array(dim = c(no_pri,3))
  res_pri$CGm = array(dim = c(no_pri,3))
  res_pri$m   = array(dim = c(no_pri))
  res_sec     = list()
  res_sec$I   = array(dim = c(3,3,no_sec))
  res_sec$CG  = array(dim = c(no_sec,3))
  res_sec$CGm = array(dim = c(no_sec,3))
  res_sec$m   = array(dim = c(no_sec))

  # determine the orientation and normal of each feather
  feather_info = orient_feather(no_pri,no_sec,Pt1,Pt2,Pt3,Pt4,Pt9,Pt10,Pt11)
  # --------------------------- Primaries -----------------------------------
  for (i in 1:no_pri){
    feather_name = paste("P",i,sep = "")
    pri_info = subset(dat_feat_curr,feather == feather_name) # subset data to be for this specific feather
    # Calculate MOI and CG
    tmp = massprop_feathers(pri_info$m_f,pri_info$l_cal,pri_info$l_vane, pri_info$w_cal,
                      dat_bird_curr$barb_radius, dat_bird_curr$barb_distance,
                      rho_cor,rho_med,
                      pri_info$w_vp,pri_info$w_vd,pri_info$vane_angle,
                      feather_info$loc_start[i,],feather_info$loc_end[i,],feather_info$normal[i,])
    # Save MOI, CG and CG*mass
    res_pri$I[,,i] = tmp$I
    res_pri$CG[i,] = tmp$CG
    res_pri$CGm[i,] = tmp$CG*pri_info$m_f
    res_pri$m[i] = tmp$m

  }
  #  -----------------------------------Secondaries -----------------------------------
  for (i in 1:no_sec){
    feather_name = paste("S",i,sep = "")
    sec_info = subset(dat_feat_curr,feather == feather_name)  # subset data to be for this specific feather
    # Calculate MOI and CG
    tmp = massprop_feathers(sec_info$m_f,sec_info$l_cal,sec_info$l_vane, sec_info$w_cal,
                            dat_bird_curr$barb_radius, dat_bird_curr$barb_distance,
                            rho_cor,rho_med,
                            sec_info$w_vp,sec_info$w_vd,sec_info$vane_angle,
                            feather_info$loc_start[i+no_pri,],feather_info$loc_end[i+no_pri,],feather_info$normal[i+no_pri,])
    # Save MOI, CG and CG*mass
    res_sec$I[,,i]  = tmp$I
    res_sec$CG[i,]  = tmp$CG
    res_sec$CGm[i,] = tmp$CG*sec_info$m_f
    res_sec$m[i]    = tmp$m
  }

  # ----------------------------------- Alula -----------------------------------
  m_alula = subset(dat_feat_curr,feather == "alula")$m_f
  alula   = massprop_pm(m_alula, Pt3)

  # ----------------------------------- Tertiaries -----------------------------------
  last_secondary_pos = feather_info$loc_end[which.min(feather_info$loc_end[,2]),]
  prop_tertiary1 = massprop_skin(0.5*dat_bird_curr$tertiary_mass,rho_cor,rbind(Pt12,Pt2,Pt11))
  prop_tertiary2 = massprop_skin(0.5*dat_bird_curr$tertiary_mass,rho_cor,rbind(Pt11,Pt2,last_secondary_pos))

  # --- All Feathers ---
  prop_feathers = list()
  I_feathers_pri = rbind(rowSums(res_pri$I[,1,]),rowSums(res_pri$I[,2,]),rowSums(res_pri$I[,3,]))
  I_feathers_sec = rbind(rowSums(res_sec$I[,1,]),rowSums(res_sec$I[,2,]),rowSums(res_sec$I[,3,]))
  prop_feathers$I  = I_feathers_pri + I_feathers_sec +alula$I + prop_tertiary1$I + prop_tertiary2$I
  prop_feathers$m  = sum(dat_feat_curr$m_f) + dat_bird_curr$tertiary_mass
  prop_feathers$CG = (colSums(res_pri$CGm) + colSums(res_sec$CGm) + alula$CG*m_alula +
                        prop_tertiary1$CG*0.5*dat_bird_curr$tertiary_mass + prop_tertiary2$CG*0.5*dat_bird_curr$tertiary_mass)/prop_feathers$m


  # ----------------------------------------------------
  # ----------------- Save Data ------------------------
  # ----------------------------------------------------

  # add data to bone specific data frame
  mass_properties_bone = store_data(dat_wingID_curr,hum,mass_properties_bone,"humerus")
  mass_properties_bone = store_data(dat_wingID_curr,ulna,mass_properties_bone,"ulna")
  mass_properties_bone = store_data(dat_wingID_curr,radius,mass_properties_bone,"radius")
  mass_properties_bone = store_data(dat_wingID_curr,car,mass_properties_bone,"carpo")
  mass_properties_bone = store_data(dat_wingID_curr,wristbone,mass_properties_bone,"wristbones")

  # add data to muscle specific data frame
  mass_properties_muscle = store_data(dat_wingID_curr,brach,mass_properties_muscle,"brach")
  mass_properties_muscle = store_data(dat_wingID_curr,abrach,mass_properties_muscle,"abrach")
  mass_properties_muscle = store_data(dat_wingID_curr,manus,mass_properties_muscle,"manus")

  # add data to muscle specific data frame
  mass_properties_skin = store_data(dat_wingID_curr,prop_skin,mass_properties_skin,"skin_prop")

  # add data to bone specific data frame
  # --- Primaries ---
  for (i in 1:no_pri){
    feather_name = paste("P",i,sep = "")
    curr_res_pri = list()
    curr_res_pri$I  = res_pri$I[,,i]
    curr_res_pri$CG = res_pri$CG[i,]
    curr_res_pri$m  = res_pri$m[i]
    mass_properties_feathers = store_data(dat_wingID_curr,curr_res_pri,mass_properties_feathers,feather_name)
  }
  # --- Secondaries ---
  for (i in 1:no_sec){
    feather_name = paste("S",i,sep = "")
    curr_res_sec = list()
    curr_res_sec$I  = res_sec$I[,,i]
    curr_res_sec$CG = res_sec$CG[i,]
    curr_res_sec$m  = res_sec$m[i]
    mass_properties_feathers = store_data(dat_wingID_curr,curr_res_sec,mass_properties_feathers,feather_name)
  }

  # save all combined data groups to the master list
  mass_properties      = store_data(dat_wingID_curr,prop_bone,mass_properties,"bones")
  mass_properties      = store_data(dat_wingID_curr,prop_muscles,mass_properties,"muscles")
  mass_properties      = store_data(dat_wingID_curr,prop_skin,mass_properties,"skin")
  mass_properties      = store_data(dat_wingID_curr,prop_feathers,mass_properties,"feathers")

  # save all wing data
  prop_wing    = list()
  prop_wing$I  = prop_bone$I + prop_muscles$I + prop_skin$I + prop_feathers$I
  prop_wing$m  = prop_bone$m + prop_muscles$m + prop_skin$m + prop_feathers$m
  prop_wing$CG = (prop_bone$CG*prop_bone$m + prop_muscles$CG*prop_muscles$m + prop_skin$CG*prop_skin$m + prop_feathers$CG*prop_feathers$m)/prop_wing$m
  mass_properties  = store_data(dat_wingID_curr,prop_wing,mass_properties,"wing")

  # Plot to verify correct outputs
  # CGplot = plot_CGloc(clean_pts,mass_properties,mass_properties_skin,mass_properties_bone,mass_properties_feathers,mass_properties_muscle)
  # plot(CGplot)

  return(mass_properties)
}


# -------------------- Store Data  -------------------------------
#' Store data
#'
#' Function to store moment of inertia tensor and center of gravity vector components in long format
#'
#' @param dat_wingID_curr Dataframe related to the current bird wing ID info that must include the following columns:
#' \itemize{
#' \item{species}{Species ID code as a string}
#' \item{WingID}{Wing ID code as a string}
#' \item{TestID}{Test ID code as a string}
#' \item{frameID}{Video frame ID code as a string}
#' }
#'
#' @param dat_mass Dataframe containing the new MOI and CG data to add to mass_properties as new rows. Must include:
#' \itemize{
#' \item{I}{Moment of inertia tensor (kg-m^2)}
#' \item{CG}{Center of gravity with three location components (m)}
#' }
#'
#' @param mass_properties Dataframe containing any previously saved data.
#' Must have the following columns: "species","WingID","TestID","FrameID","prop_type","component","value".
#'
#' @param name Name of the component for which the moment of inertia and center of gravity were computed.
#'
#' @return This function returns mass_properties as an updated dataframe with a new row corresponding to the dat_mass information
#'
#' @export
store_data <- function(dat_wingID_curr,dat_mass,mass_properties,name){
  prop_type_list = c("Ixx","Iyy","Izz","Ixy","Iyz","Ixz","CGx","CGy","CGz")
  prop_type_ind  = cbind(c(1,2,3,1,2,1),c(1,2,3,2,3,3))

  species = dat_wingID_curr$species
  wingID  = dat_wingID_curr$WingID
  testID  = dat_wingID_curr$TestID
  frameID = dat_wingID_curr$frameID

  # Moment of inertia tensor
  for (i in 1:6){
    new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = prop_type_list[i], value = dat_mass$I[prop_type_ind[i,1],prop_type_ind[i,2]]) # saves the name and valueof the tensor component
    mass_properties = rbind(mass_properties,new_row)
  }

  # Center of gravity
  for (i in 1:3){
    new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = prop_type_list[6+i], value = dat_mass$CG[i]) # saves the name and value of the CG component
    mass_properties = rbind(mass_properties,new_row)
  }

  # Mass
  new_row = data.frame(species = species,WingID = wingID,TestID = testID, FrameID = frameID, component = name,
                         object = "m", value = dat_mass$m) # saves the name and value of the mass
  mass_properties = rbind(mass_properties,new_row)



  return(mass_properties)
}


