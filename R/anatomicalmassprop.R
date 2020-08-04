# Script containing the functions to calculate the moment of inertia and CG of different anatomical components
# all written by Christina Harvey
# last updated: 2020-08-03

##### ---------------------- Mass properties of Bone -------------------------- #####
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
  # ---------------------------------------------------------------------------------
  # ---------------------------- Define geometry ------------------------------------
  # ---------------------------------------------------------------------------------
  vol   = m/rho # total volume of the bone

  # calculate how thick the cap needs to be to match the vol1 = m/rho and vol2 = geometric
  cap_t = (-(0.5*l*((r_out^2)-(r_in^2)))/(r_in^2))+(vol/(2*pi*(r_in^2)))
  cap_m = pi*cap_t*cap_r^2

  cyl_m = m - (2*cap_m)
  cyl_l = l - (2*cap_t)

  # ---------------------------------------------------------------------------------
  # ------------------------------- Adjust axis -------------------------------------
  # ---------------------------------------------------------------------------------
  z_axis = start;
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = pracma::cross(z_axis,temp_vec/pracma::Norm(temp_vec));

  # doesn't matter where the x axis points as long as:
  # 1. we know what it is 2. it's orthogonal to z
  VRP2object = calc_rot(z_axis,x_axis);

  # ---------------------------------------------------------------------------------
  # ----------------------- Calculate moment of inertia -----------------------------
  # ---------------------------------------------------------------------------------

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


##### ---------------------- Mass properties of Muscle -------------------------- #####
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

  # ---------------------------------------------------------------------------------
  # ------------------------------- Adjust axis -------------------------------------
  # ---------------------------------------------------------------------------------
  z_axis = end-start
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = pracma::cross(z_axis,temp_vec/pracma::Norm(temp_vec))

  # doesn't matter where the x axis points as long as:
  # 1. we know what it is 2. it's orthogonal to z
  VRP2object = calc_rot(z_axis,x_axis)

  # -----------------------------------------------------------------------------
  # -------------------------- Moment of inertia --------------------------------
  # -----------------------------------------------------------------------------

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


##### ---------------------- Mass properties of Skin -------------------------- #####
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
  # ---------------------------------------------------------------------------------
  # ------------------ Determine the geometry of the skin ---------------------------
  # ---------------------------------------------------------------------------------
  temp_cross = pracma::cross((pts[2,]-pts[1,]),(pts[2,]-pts[3,]))

  A = 0.5*pracma::Norm(temp_cross);
  v = m/rho
  t = v/A;
  # ---------------------------------------------------------------------------------
  # ------------------------------- Adjust axis -------------------------------------
  # ---------------------------------------------------------------------------------
  z_axis = temp_cross
  temp_vec = pts[2,]-pts[1,] # vector points towards the second input point along the bone edge
  x_axis = pracma::cross(z_axis,temp_vec/pracma::Norm(temp_vec))

  # doesn't matter where the x axis points as long as:
  # 1. we know what it is 2. it's orthogonal to z
  VRP2object = calc_rot(z_axis,x_axis)

  # Adjust the input pts frame to skin axes
  adj_pts = pts # define the new matrix
  for (i in 1:4){
    adj_pts[i,] = VRP2object*pts[i,]                      # Frame of reference: Skin | Origin: VRP
  }

  # -----------------------------------------------------------------------------
  # ---------------------------- Moment of inertia ------------------------------
  # -----------------------------------------------------------------------------

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


##### ---------------------- Mass properties of Feathers -------------------------- #####
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
  # ---------------------------------------------------------------------------------
  # ------------------ Determine the geometry of muscles ----------------------------
  # ---------------------------------------------------------------------------------
  temp_cross = pracma::cross((pts[2,]-pts[1,]),(pts[2,]-pts[3,]))

  A = 0.5*pracma::Norm(temp_cross);
  v = m/rho
  t = v/A;
  # ---------------------------------------------------------------------------------
  # ------------------------------- Adjust axis -------------------------------------
  # ---------------------------------------------------------------------------------
  z_axis = temp_cross
  temp_vec = pts[2,]-pts[1,] # vector points towards the second input point along the bone edge
  x_axis = pracma::cross(z_axis,temp_vec/pracma::Norm(temp_vec))

  # doesn't matter where the x axis points as long as:
  # 1. we know what it is 2. it's orthogonal to z
  VRP2object = calc_rot(z_axis,x_axis)

  # Adjust the input pts frame to skin axes
  adj_pts = pts # define the new matrix
  for (i in 1:4){
    adj_pts[i,] = VRP2object*pts[i,]                      # Frame of reference: Skin | Origin: VRP
  }

  # -----------------------------------------------------------------------------
  # ----------------------- Moment of inertia - muscle axis----------------------
  # -----------------------------------------------------------------------------

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
