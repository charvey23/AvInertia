# Script containing the functions to calculate the moment of inertia and CG of different anatomical components
# all written by Christina Harvey
# last updated: 2020-08-03

# ---------------------------------------------------------------------------------------
##### ------------------------ Mass properties of Bone ---------------------------- #####
# ---------------------------------------------------------------------------------------

#' Calculate the moment of inertia of a bone modelled as a hollow cylinder with two solid end caps
#'
#' @param m full mass of the bone
#' @param l full length of the bone
#' @param r_out outer radius of bone
#' @param r_in inner radius of bone
#' @param rho density of the bone (kg/m^3)
#' @param start 3D point where bone starts. Frame of reference: VRP | Origin: VRP
#' @param end 3D point where bone ends. Frame of reference: VRP | Origin: VRP
#'
#' @return
#'
#' @section Warning:
#' the parallel axis theorem only works if the moment of inertia is calculated about components centroidal axis
#'
#' @export
#'
#' @examples
massprop_bones <- function(m,l,r_out,r_in,rho,start,end){
  # ---------------------------- Define geometry ------------------------------------
  vol   = m/rho # total volume of the bone

  # calculate how thick the cap needs to be to match the vol1 = m/rho and vol2 = geometric
  cap_t = (-(0.5*l*((r_out^2)-(r_in^2)))/(r_in^2))+(vol/(2*pi*(r_in^2)))
  cap_m = pi*cap_t*cap_r^2

  cyl_m = m - (2*cap_m)
  cyl_l = l - (2*cap_t)

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = start;
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = crossprod(z_axis,temp_vec/norm(temp_vec, type = "2"))

  # doesn't matter where the x axis points as long as: 1. we know what it is 2. it's orthogonal to z
  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # ----------------------- Calculate moment of inertia -----------------------------

  I_cyl_b = calc_inertia_cylhollow(r_out, r_in, cyl_l, cyl_m) # Frame of reference: Bone | Origin: Bone CG
  I_cap_b = calc_inertia_cylsolid(r_out, cap_t, cap_m)        # Frame of reference: Bone | Origin: Bone CG

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off = VRP2object %*% start                                  # Frame of reference: Bone | Origin: VRP

  # determine the offset vector for each component with       # Frame of reference: Bone | Origin: VRP
  I_cyl_off   = c(0,0,0.5*l_bone) + off            # the hollow cylinder is displaced halfway along the bone (in z-axis)
  I_cap1_off  = c(0,0,0.5*cap_t) + off             # Cap 1 edge centered on the beginning of the bone
  I_cap2_off  = c(0,0,(l_bone - (0.5*cap_t)))+ off # Cap 2 edge centered on the end of the bone

  # need to adjust the moment of inertia tensor               # Frame of reference: Bone | Origin: VRP
  I_cyl_vrp   = parallelaxis(I_cyl,-I_cyl_off,m_cyl)
  I_cap1_vrp  = parallelaxis(I_cap,-I_cap1_off,m_cap)
  I_cap2_vrp  = parallelaxis(I_cap,-I_cap2_off,m_cap)

  I_boneaxis  = I_cyl_vrp + I_cap1_vrp + I_cap2_vrp           # Frame of reference: Bone | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_boneaxis %*% VRP2object  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(start + end)                            # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}

# ---------------------------------------------------------------------------------------
##### ------------------------ Mass properties of Muscle -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Calculate the moment of inertia of a muscle modelled as a solid cylinder distributed along the bone length
#'
#' @param m mass of muscle (kg)
#' @param l length of muscle (m)
#' @param rho density of muscle (kg/m^3)
#' @param start 3D point where bone starts. Frame of reference: VRP | Origin: VRP
#' @param end 3D point where bone ends. Frame of reference: VRP | Origin: VRP
#'
#' @return
#' @export
#'
#' @examples
massprop_muscles <- function(m,l,rho,start,end){

  r = sqrt(m/(rho*pi*l)) # determine the pseudo radius of the muscles

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = crossprod(z_axis,temp_vec/norm(temp_vec, type = "2"))
  # doesn't matter where the x axis points as long as:1. we know what it is 2. it's orthogonal to z
  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # -------------------------- Moment of inertia --------------------------------

   I_m = calc_inertia_cylsolid(r, l, m)                    # Frame of reference: Muscle | Origin: Muscle CG

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off = VRP2object %*% start                               # Frame of reference: Muscle | Origin: VRP

  # determine the offset vector for each component
  I_m_off  = c(0,0,0.5*l)+ off                             # Frame of reference: Muscle | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_m_vrp   = parallelaxis(I_m,-I_m_off,m)                 # Frame of reference: Muscle | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_m_vrp %*% VRP2object  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(start + end)                         # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}

# ---------------------------------------------------------------------------------------
##### ------------------------- Mass properties of Skin --------------------------- #####
# ---------------------------------------------------------------------------------------

#' Calculate the moment of inertia of skin modelled as a flat triangular plate
#'
#' @param m mass of skin (kg)
#' @param rho density of skin (kg/m^3)
#' @param pts three points defining the vertices of the triangle. Frame of reference: VRP | Origin: VRP
#' Must be numbered in a counterclockwise direction for positive area.
#' each point should be a different row as follows:
#' pt1x, pt1y, pt1z
#' pt2x, pt1y, pt2z
#' pt3x, pt3y, pt3z
#'
#' @author Christina Harvey
#'
#' @return
#' @export
#'
#' @examples
massprop_skin <- function(m,rho,pts){
  # ------------------ Determine the geometry of the skin ---------------------------
  temp_cross = crossprod((pts[2,]-pts[1,]),(pts[2,]-pts[3,]))

  A = 0.5*norm(temp_cross, type = "2");
  v = m/rho
  t = v/A;

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = temp_cross
  temp_vec = pts[2,]-pts[1,] # vector points towards the second input point along the bone edge
  x_axis = crossprod(z_axis,temp_vec/norm(temp_vec, type = "2"))

  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # Adjust the input pts frame to skin axes
  adj_pts = pts # define the new matrix
  for (i in 1:4){
    adj_pts[i,] = VRP2object*pts[i,]                      # Frame of reference: Skin | Origin: VRP
  }

  # ---------------------------- Moment of inertia ------------------------------

  I_s = calc_inertia_platetri(adj_pts, A, rho, t, "I")    # Frame of reference: Skin | Origin: Skin CG
  CG_s = calc_inertia_platetri(adj_pts, A, rho, t, "CG")  # Frame of reference: Skin | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_s_vrp   = parallelaxis(I_s,-CG_s,m)                   # Frame of reference: Skin | Origin: VRP

  mass_prop = list() # pre-define

  # Adjust frames to VRP
  mass_prop$I  = t(VRP2object) %*% I_s_vrp %*% VRP2object # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_s                   # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### --------------------- Mass properties of a point mass ----------------------- #####
# ---------------------------------------------------------------------------------------

#' Calculate the moment of inertia of of any point mass
#'
#' @param m mass of skin (kg)
#' @param pt density of skin (kg/m^3)
#' @param pt three points defining the location of the point mass. Frame of reference: VRP | Origin: VRP
#' each point should be a different row as follows:
#' pt1x, pt1y, pt1z
#'
#' @author Christina Harvey
#'
#' @return
#' @export
#'
#' @examples
massprop_pm <- function(m,pt){
  # Adjust frames to VRP
  mass_prop$I  = parallelaxis(0,-pt,m)  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = pt                     # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### ---------------------- Mass properties of Feathers -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Calculate the moment of inertia of skin modelled as a flat triangular plate
#'
#' @param m mass of skin (kg)
#' @param start 3D point where feather starts. Frame of reference: VRP | Origin: VRP
#' @param end 3D point where feather ends Frame of reference: VRP | Origin: VRP
#' @param n_pts the pts that define the plane that the feather lies flat on.
#' Order counterclockwise so that the normal will pointing downwards from the ventral surface of the feather.
#' Frame of reference: VRP | Origin: VRP
#' @param l_v length of the feather vane - tip to end of calamus (m)
#' @param l_c length of the calamus - start of vane to end of calamus(m)
#' @param l_r length of rachis - tip to start of vane (m)
#' @param r_cor radius of the cortex part of the calamus (m)
#' @param r_med radius of the medullary part of the calamus (m)
#' @param r_b barb radius (m)
#' @param d_b barb distance (m)
#' @param m_f mass of the entire feather
#' @param rho_cor density of the cortex (kg/m^3)
#' @param rho_med density of the medullary (kg/m^3)
#' @param w_vp width of proximal (closest to body) vane (m)
#' @param w_vd width of distal (closest to wing tip) vane (m)
#' @param angle angle between l_c and l_r
#'
#' @author Christina Harvey
#'
#' @return
#' @export
#'
#' @examples
massprop_feathers <- function(m,rho,start,end,n_pts,l_v,l_c,l_r,r_cor,r_med,r_b,d_b,m_f,rho_cor,rho_med, w_vp, w_vd,angle){

  # ------------------ Determine the geometry of feathers ---------------------------

  # mass of each component of the feather
  m_vp = rho_cor*(l_v/d_b)*w_vp*pi*r_b^2 # mass of the proximal vane
  m_vd = rho_cor*(l_v/d_b)*w_vd*pi*r_b^2 # mass of the distal vane
  m_r  = m_f - m_vp - m_vd               # mass of the rachis and calamus

  # determine the width of the medullary material inside the rachis
  c_0 = (rho_cor*(r_cor^2)*((pi*l_c) + (4*l_v/3))) - m_r
  c_2 = (pi*l_c*(rho_med-rho_cor))
  c_3 = (4*l_v/(3*r_cor))*(rho_med - rho_cor)
  # solve the roots
  cubic_roots = polyroot(c(c_0,0,c_2,c_3))

  #save the root that is real and positive
  for (i in 1:length(test)){
    if (Imag(test[i]) == 0 & Re(test[i]) > 0 & Re(test[i]) < r_cor){
      r_med = as.numeric(Re(test[i])) # radius of the medullary component within the calamus
    }
  }

  # calculate the length of the medullary material within the rachis (height of that pyramid)
  l_r_med    = l_v*(r_med/r_cor)

  # calculate the mass of each component
  m_c_cor = rho_cor*(pi*l_c)*(r_cor^2-r_med^2)                              # mass of the cortex part of the calamus
  m_c_med = rho_med*(pi*l_c)*(r_med^2)                                      # mass of the medullary part of the calamus
  m_r_cor = rho_cor*((4/3)*l_v*r_cor^2) - rho_cor*((4/3)*l_r_med)*(r_med^2) # mass of the cortex part of the rachis
  m_r_med = rho_med*((4/3)*l_r_med)*(r_med^2)                               # mass of the medullary part of the rachis

  # ------------------------------- Adjust axis -------------------------------------

  z_axis = end - start
  temp_vec = crossprod((pt[2,]-pt[1,]),(pt[3,]-pt[1,])) # normal vector of the feather based on supplied data points
  x_axis = crossprod(z_axis,temp_vec/norm(temp_vec, type = "2"))
  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # Adjust the input pts frame to skin axes
  adj_pts = pts # define the new matrix
  for (i in 1:4){
    adj_pts[i,] = VRP2object*pts[i,]                      # Frame of reference: Skin | Origin: VRP
  }

  # --------------------------- Moment of inertia -----------------------------------

  I_s = calc_inertia_platetri(adj_pts, A, rho, t, "I")    # Frame of reference: Skin | Origin: Skin CG
  CG_s = calc_inertia_platetri(adj_pts, A, rho, t, "CG")  # Frame of reference: Skin | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_s_vrp   = parallelaxis(I_s,-CG_s,m)                   # Frame of reference: Skin | Origin: VRP

  mass_prop = list() # pre-define

  # Adjust frames to VRP
  mass_prop$I  = t(VRP2object) %*% I_s_vrp %*% VRP2object # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_s                   # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}
