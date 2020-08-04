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
  r_cap = r_out
  t_cap = (-(0.5*l*((r_out^2)-(r_in^2)))/(r_in^2))+(vol/(2*pi*(r_in^2)))
  m_cap = pi*t_cap*r_cap^2

  # cylinder geometry
  m_cyl = m - (2*m_cap)
  l_cyl = l - (2*t_cap)

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start;
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))

  # doesn't matter where the x axis points as long as: 1. we know what it is 2. it's orthogonal to z
  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # ----------------------- Calculate moment of inertia -----------------------------

  I_cyl = calc_inertia_cylhollow(r_out, r_in, l_cyl, m_cyl) # Frame of reference: Bone | Origin: Bone CG
  I_cap = calc_inertia_cylsolid(r_cap, t_cap, m_cap)        # Frame of reference: Bone | Origin: Bone CG

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off = VRP2object%*%start                                  # Frame of reference: Bone | Origin: VRP

  # determine the offset vector for each component with       # Frame of reference: Bone | Origin: VRP
  I_cyl_off   = c(0,0,0.5*l) + off                 # the hollow cylinder is displaced halfway along the bone (in z-axis)
  I_cap1_off  = c(0,0,0.5*t_cap) + off             # Cap 1 edge centered on the beginning of the bone
  I_cap2_off  = c(0,0,(l - (0.5*t_cap)))+ off      # Cap 2 edge centered on the end of the bone

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
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))
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
  temp_cross = pracma::cross((pts[2,]-pts[1,]),(pts[2,]-pts[3,]))

  A = 0.5*norm(temp_cross, type = "2");
  v = m/rho
  t = v/A;

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = temp_cross
  x_axis = pts[2,]-pts[1,] # vector points towards the second input point along the bone edge

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
#' @param l length of the feather - tip to end of calamus (m)
#' @param l_c length of the calamus - start of vane to end of calamus(m)
#' @param l_r_cor length of rachis/vane - tip to start of vane (m)
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
#' @section Warning:
#' Parallel axis theorem is only valid if relocating an object between it's CG and an arbitrary pt. I_arbitrary = I_CG + md^2
#' Cannot be used between two arbitrary points.
#'
#' @return
#' @export
#'
#' @examples
massprop_feathers <- function(m,rho,start,end,n_pts,l,l_c,l_r_cor,r_cor,r_med,r_b,d_b,m_f,rho_cor,rho_med, w_vp, w_vd,angle){

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
  x_axis = pracma::cross((pt[2,]-pt[1,]),(pt[3,]-pt[1,]), type = "2") # normal vector of the feather based on supplied data points

  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)
  off = VRP2object*start

  # --------------------------- Moment of inertia -----------------------------------
  # ------- Calamus -------
  # Moment of inertia tensors in the calamus
  I_c_cor = calc_inertia_cylhollow(r_cor, r_med, l_c, m_c_cor)    # Frame of reference: Feather Calamus | Origin: Calamus CG
  I_c_med = calc_inertia_cylsolid(r_med, l_c, m_c_med)            # Frame of reference: Feather Calamus | Origin: Calamus CG

  # ------- Rachis -------
  # Moment of inertia tensor in the rachis - inner medullary component
  I_r_in_med  = calc_inertia_pyrasolid(r_med, l_r_med, m_c_med)       # Frame of reference: Feather Vane | Origin: Inner Rachis Medullary hollow pyramid CG

  # Moment of inertia tensor in the rachis - as if solid inner and outer cortex components
  I_r_out_cor = calc_inertia_pyrasolid(r_cor, l_r_cor, m_r_cor)       # Frame of reference: Feather Vane | Origin: Start of the vane (center)
  I_r_in_cor  = calc_inertia_pyrasolid(r_cor, l_r_med, m_c_med)       # Frame of reference: Feather Vane | Origin: Start of the vane (center)
  I_r_cor = I_r_out_cor - I_r_in_cor # hollow square pyramid          # Frame of reference: Feather Vane | Origin: Start of the vane (center)

  # mass of the solid inner and outer cortex components
  mass_inner     = (4*rho_cor*r_med^2*l_r_med/3) # mass of solid inner cortex pyramid
  mass_outer     = (4*rho_cor*r_cor^2*l_r_cor/3) # mass of solid outer cortex pyramid

  # Calculate the CG of hollow cortex pyramid and solid medullary pyramid
  CG_pyrasquare = c(0,0,0.25) # multiple by length
  CG_r_out_cor  = CG_pyrasquare*l_r_cor # Frame of reference: Feather Vane | Origin: Start of the vane (center)
  CG_r_in_cor   = CG_pyrasquare*l_r_med # Frame of reference: Feather Vane | Origin: Start of the vane (center)
  CG_r_cor1 = (mass_outer*CG_r_out_cor - mass_inner*CG_r_in_cor)/m_r # Frame of reference: Feather Vane | Origin: Start of the vane (center)
  CG_r_med1 = CG_pyrasquare*l_r_med # Frame of reference: Feather Vane | Origin: Start of the vane

  # To properly use parallel axis need to relocate origin to the CG of the hollow pyramid
  I_r_cor = parallelaxis(I_r_cor,CG_r_cor,m_r_cor)        # Frame of reference: Feather Vane | Origin: Outer Rachis Cortex hollow pyramid CG

  # ------- Vane -------
  I_vd = calc_inertia_platerect(w_vd, l_r_cor, m_vd)    # Frame of reference: Feather Distal Vane | Origin: Distal Vane CG
  CG_vd_prerot1 = c(0,-0.5*w_vd,0.5*l_r_cor)            # Frame of reference: Feather Vane | Origin: Start of the vane (edge)

  I_vp = calc_inertia_platerect(w_vp, l_r_cor, m_vp)    # Frame of reference: Feather Proximal Vane | Origin: Proximal Vane CG
  CG_vp_prerot1 = c(0,0.5*w_vp,0.5*l_r_cor)             # Frame of reference: Feather Vane | Origin: Start of the vane (edge)


  # ------- Adjust the Vanes and Rachis to be Feather Calamus frame of reference -------
  # postive angle is a ccw rotation about x (normal to the feather plane) ** all angles should be negative for bird feathers
  rot_r  = rotx(angle)
  rot_vd = rotx(angle - atand(r_cor/l_r_cor))
  rot_vp = rotx(angle + atand(r_cor/l_r_cor))

  # rotate all rachis and vane MOI to be in the same frame of reference as the calamus
  I_r_cor = rot_r %*% I_r_cor %*% t(rot_r)    # Frame of reference: Feather Calamus | Origin: Outer Rachis Cortex hollow pyramid CG
  I_r_med = rot_r %*% I_r_med %*% t(rot_r)    # Frame of reference: Feather Calamus | Origin: Inner Rachis Medullary hollow pyramid CG
  I_vd    = rot_vd %*% I_vd %*% t(rot_vd)     # Frame of reference: Feather Calamus | Origin: Distal Vane CG
  I_vp    = rot_vp %*% I_vp %*% t(rot_vp)     # Frame of reference: Feather Calamus | Origin: Proximal Vane CG

  # rotate all rachis CG to be in the same frame of reference as the calamus
  CG_r_cor_rot   = rot_r*CG_r_cor1              # Frame of reference: Feather Calamus | Origin: Start of the vane (center of pyramid base)
  CG_r_med_rot   = rot_r*CG_r_med1              # Frame of reference: Feather Calamus | Origin: Start of the vane (center of pyramid base)

  # rotate all vane CG about the lower corner close to the rachis to be in the same frame of reference as the calamus need to also relocate
  CG_vd_prerot2  =  rotx(-atand(r_c/l_v))*CG_vd_prerot1 # Frame of reference: Feather Vane | Origin: Start of the vane (edge)
  CG_vp_prerot2  =  rotx(atand(r_c/l_v))*CG_vp_prerot1  # Frame of reference: Feather Vane | Origin: Start of the vane (edge)

  CG_vd_prerot3  =  CG_vd_prerot2 - c(0,r_c,0)          # Frame of reference: Feather Vane | Origin: Start of the vane (center)
  CG_vp_prerot3  =  CG_vp_prerot2 + c(0,r_c,0)          # Frame of reference: Feather Vane | Origin: Start of the vane (center)

  CG_vd_postrot =  rot_vd*CG_vd_prerot3                 # Frame of reference: Feather Calamus | Origin: Start of the vane (center)
  CG_vp_postrot =  rot_vp*CG_vp_prerot3                 # Frame of reference: Feather Calamus | Origin: Start of the vane (center)

  # ---- Determine the offset between the start of the feather and the current I origin ----
  CG_c     = c(0,0,0.5*l_c) + off                       # Frame of reference: Feather Calamus | Origin: VRP
  CG_r_cor = c(0,0,l_c) + CG_r_cor_rot + off            # Frame of reference: Feather Calamus | Origin: VRP
  CG_r_med = c(0,0,l_c) + CG_r_med_rot + off            # Frame of reference: Feather Calamus | Origin: VRP
  CG_vd    = c(0,0,l_c) + CG_vd_postrot + off           # Frame of reference: Feather Calamus | Origin: VRP
  CG_vp    = c(0,0,l_c) + CG_vp_postrot + off           # Frame of reference: Feather Calamus | Origin: VRP

  CG_calaxis = (m_c_cor*CG_c + m_c_med*CG_c + m_r_cor*CG_r_cor + m_r_med*CG_r_med + m_vd*CG_vd + m_vp*CG_vp)/m_f # Frame of reference: Feather Calamus | Origin: VRP

  # ---- Relocate the MOI from object CG to VRP ----
  I_c_cor = parallelaxis(I_c_cor,-CG_c,m_c_cor)         # Frame of reference: Feather Calamus | Origin: VRP
  I_c_med = parallelaxis(I_c_med,-CG_c,m_c_med)         # Frame of reference: Feather Calamus | Origin: VRP
  I_r_cor = parallelaxis(I_r_cor,-CG_r_cor,m_r_cor)     # Frame of reference: Feather Calamus | Origin: VRP
  I_r_med = parallelaxis(I_r_med,-CG_r_med,m_r_med)     # Frame of reference: Feather Calamus | Origin: VRP
  I_vd    = parallelaxis(I_vd,-CG_vd,m_vd)              # Frame of reference: Feather Calamus | Origin: VRP
  I_vp    = parallelaxis(I_vp,-CG_vp,m_vp)              # Frame of reference: Feather Calamus | Origin: VRP

  I_calaxis = I_c_cor + I_c_med + I_r_cor + I_r_med + I_vd + I_vp # Frame of reference: Feather Calamus | Origin: VRP

  # ---- Adjust frames to VRP ----
  mass_prop = list() # pre-define
  mass_prop$I  = t(VRP2object) %*% I_calaxis %*% VRP2object # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_calaxis               # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}
