# Script containing the functions to calculate the moment of inertia and CG of different anatomical components
# all written by Christina Harvey
# last updated: 2020-09-10

# ---------------------------------------------------------------------------------------
##### ------------------------ Mass properties of Bone ---------------------------- #####
# ---------------------------------------------------------------------------------------

#' Bone mass properties
#'
#' Calculate the moment of inertia of a bone modeled as a hollow cylinder with two solid end caps
#'
#' @param m Mass of bone (kg)
#' @param l Length of bone (kg)
#' @param r_out Outer radius of bone (m)
#' @param r_in Inner radius of bone (m)
#' @param rho Density of the bone (kg/m^3)
#' @param start a 1x3 vector (x,y,z) representing the 3D point where bone starts. Points from the VRP to the bone start. Frame of reference: VRP | Origin: VRP
#' @param end a 1x3 vector (x,y,z) representing the 3D point where bone ends. Points from the VRP to the bone start. Frame of reference: VRP | Origin: VRP
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a bone modeled as a hollow cylinder with two solid end caps}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a bone modeled as a hollow cylinder with two solid end caps}
#' \item{m}{a double that returns the input bone mass}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export
#'
#' @examples
#'
#' # This example takes the humerus data contained within the package to predict
#' # the humerus moment of inertia and center of gravity with the origin at the VRP and in the VRP frame of reference.
#'
#' library(birdmoment)
#'
#' massprop_bones(m,l,r_out,r_in,rho,start,end)
#'
#'
massprop_bones <- function(m,l,r_out,r_in,rho,start,end){
  # ---------------------------- Define geometry ------------------------------------
  vol   = m/rho # total volume of the bone material

  # calculate how thick the cap needs to be to match the vol_true = m/rho and vol2 = geometric
  r_cap = r_out
  t_cap = 0.5*((-(l*(r_out^2-r_in^2))/(r_in^2))+(vol/(pi*r_in^2)))
  m_cap = pi*t_cap*r_cap^2

  # cylinder geometry
  m_cyl = m - (2*m_cap)
  l_cyl = l - (2*t_cap)

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start;                   # Frame of reference: VRP | Origin: VRP
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))
  # doesn't matter where the x axis points as long as: 1. we know what it is 2. it's orthogonal to z
  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)  # Frame of reference: VRP | Origin: VRP

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
  I_cyl_vrp   = parallelaxis(I_cyl,-I_cyl_off,m_cyl,"CG")
  I_cap1_vrp  = parallelaxis(I_cap,-I_cap1_off,m_cap,"CG")
  I_cap2_vrp  = parallelaxis(I_cap,-I_cap2_off,m_cap,"CG")

  I_boneaxis  = I_cyl_vrp + I_cap1_vrp + I_cap2_vrp           # Frame of reference: Bone | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_boneaxis %*% VRP2object  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(start + end)                            # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m

  return(mass_prop)
}



# ---------------------------------------------------------------------------------------
##### ------------------------ Mass properties of Muscle -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Muscle mass properties
#'
#' Calculate the moment of inertia of a muscle modeled as a solid cylinder distributed along the bone length
#'
#' @param m Mass of muscle (kg)
#' @param rho Density of muscle (kg/m^3)
#' @param start a 1x3 vector (x,y,z) representing the 3D point where bone starts. Frame of reference: VRP | Origin: VRP
#' @param end a 1x3 vector (x,y,z) representing the 3D point where bone ends. Frame of reference: VRP | Origin: VRP
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a muscle modeled as a solid cylinder distributed along the bone length}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a muscle modeled as a solid cylinder distributed along the bone length}
#' \item{m}{a double that returns the input muscle mass}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export

massprop_muscles <- function(m,rho,start,end){

  l = norm((end - start), type = "2")
  r = sqrt(m/(rho*pi*l)) # determine the pseudo radius of the muscles

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start
  temp_vec = c(1,1,1) # arbitrary vector as long as it's not the z-axis
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))
  # doesn't matter where the x axis points as long as:1. we know what it is 2. it's orthogonal to z
  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # -------------------------- Moment of inertia --------------------------------

  I_m = calc_inertia_cylsolid(r, l, m)                     # Frame of reference: Muscle | Origin: Muscle CG

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off = VRP2object %*% start                               # Frame of reference: Muscle | Origin: VRP

  # determine the offset vector for each component
  I_m_off  = c(0,0,0.5*l)+ off                             # Frame of reference: Muscle | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_m_vrp   = parallelaxis(I_m,-I_m_off,m,"CG")            # Frame of reference: Muscle | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_m_vrp %*% VRP2object  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(start + end)                         # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m
  return(mass_prop)
}

# ---------------------------------------------------------------------------------------
##### ------------------------- Mass properties of Skin --------------------------- #####
# ---------------------------------------------------------------------------------------

#' Skin mass properties
#'
#' Calculate the moment of inertia of skin/ tertiaries modeled as a flat triangular plate
#'
#' @param m Mass of skin (kg)
#' @param rho Density of skin (kg/m^3)
#' @param pts a 3x3 matrix that represent three points that define the vertices of the triangle. Frame of reference: VRP | Origin: VRP
#' Must be numbered in a counterclockwise direction for positive area, otherwise signs will be reversed.
#' each point should be a different of the matrix as follows:
#' \itemize{
#' \item{pt1x, pt1y, pt1z}
#' \item{pt2x, pt1y, pt2z}
#' \item{pt3x, pt3y, pt3z}
#' }
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' point mass
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of skin modeled as a flat triangular plate}
#' \item{CG}{a 1x3 vector representing the center of gravity position of skin modeled as a flat triangular plate}
#' \item{m}{a double that returns the input skin mass}
#' }
#'
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' Caution: The skin frame of reference assumes that the z axis is normal to the incoming points.
#'
#' @export

massprop_skin <- function(m,rho,pts){
  # ------------------ Determine the geometry of the skin ---------------------------
  temp_cross = pracma::cross((pts[3,]-pts[1,]),(pts[3,]-pts[2,]))

  A = 0.5*norm(temp_cross, type = "2");
  v = m/rho
  t = v/A;

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = temp_cross      # normal to the incoming points
  x_axis = pts[3,]-pts[1,] # vector points towards the second input point along the bone edge

  # calculate the rotation matrix between VRP frame of reference and the object
  VRP2object = calc_rot(z_axis,x_axis)

  # Adjust the input pts frame to skin axes
  adj_pts = pts # define the new matrix
  for (i in 1:3){
    adj_pts[i,] = VRP2object%*%pts[i,]                    # Frame of reference: Skin | Origin: VRP
  }

  # ---------------------------- Moment of inertia ------------------------------

  I_s  = calc_inertia_platetri(adj_pts, A, rho, t, "I")   # Frame of reference: Skin | Origin: Skin CG
  CG_s = calc_inertia_platetri(adj_pts, A, rho, t, "CG")  # Frame of reference: Skin | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_s_vrp   = parallelaxis(I_s,-CG_s,m,"CG")              # Frame of reference: Skin | Origin: VRP

  mass_prop = list() # pre-define

  # Adjust frames to VRP
  mass_prop$I  = t(VRP2object) %*% I_s_vrp %*% VRP2object # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_s                   # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m
  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### --------------------- Mass properties of a point mass ----------------------- #####
# ---------------------------------------------------------------------------------------

#' Point-mass mass properties
#'
#' Calculate the moment of inertia of any point mass
#'
#' @param m Mass of skin (kg)
#' @param pt Density of skin (kg/m^3)
#' @param pt a 1x3 vector (x,y,z) representing the  location of the point mass. Frame of reference: VRP | Origin: VRP
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a point mass}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a point mass}
#' \item{m}{a double that returns the input mass}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export

massprop_pm <- function(m,pt){

  emtpy_I = matrix(0, nrow = 3, ncol = 3) # 0 tensor for point mass about it's own CG

  mass_prop = list() # pre-define

  # Adjust frames to VRP
  mass_prop$I  = parallelaxis(emtpy_I,-pt,m,"CG")  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = pt                                # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m

  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### ---------------------- Mass properties of Feathers -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Feather mass properties
#'
#' Calculate the moment of inertia of skin modeled as a flat triangular plate
#'
#' @param m_f Mass of the entire feather (kg)
#' @param l_c Length of the calamus; start of vane to end of calamus(m)
#' @param l_r_cor Length of rachis; tip to start of vane (m)
#' @param w_cal Width (diameter) of the cortex part of the calamus (m)
#' @param r_b Radius of feather barbs (m)
#' @param d_b Distance between barbs (m)
#' @param rho_cor Density of the cortex (kg/m^3)
#' @param rho_med Density of the medullary (kg/m^3)
#' @param w_vp Width of proximal (closest to body) vane (m)
#' @param w_vd Width of distal (closest to wing tip) vane (m)
#' @param angle Angle between calamus and the vane taken the supplement angle to the interior angle. Negative indicates the feather tip is rotated proximally relative to the start of the feather vane.
#' @param start a 1x3 vector (x,y,z) representing the 3D point where feather starts. (Frame of reference: VRP | Origin: VRP)
#' @param end a 1x3 vector (x,y,z) representing the 3D point where feather tip ends (Frame of reference: VRP | Origin: VRP)
#' @param normal Vector that defines the normal to each feather plane. (Frame of reference: VRP | Origin: VRP)
#' @author Christina Harvey
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' CAUTION: While computing the variable components of the feather the x axis is the normal of the feather.
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a simplified feather}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a simplified feather}
#' \item{m}{a double that returns the feather mass}
#' }
#'
#' @export

massprop_feathers <- function(m_f,l_c,l_r_cor,w_cal,r_b,d_b,rho_cor,rho_med,w_vp,w_vd,angle){
  # ------------------ Determine the geometry of feathers ---------------------------
  r_cor = 0.5*w_cal # radius of the cortex part of the calamus

  # mass of each component of the feather
  m_vp = rho_cor*(l_r_cor/d_b)*w_vp*pi*r_b^2 # mass of the proximal vane
  m_vd = rho_cor*(l_r_cor/d_b)*w_vd*pi*r_b^2 # mass of the distal vane
  m_rc = m_f - m_vp - m_vd                   # mass of the rachis and calamus

  #assumes that the interior medullary pyramid has the same height as the cortex exterior
  r_med   = sqrt((m_rc -  rho_cor*r_cor^2*((pi*l_c)+(4/3)*l_r_cor))/((4/3)*l_r_cor*(rho_med-rho_cor)-l_c*pi*rho_cor))
  l_r_med = l_r_cor
  # calculate the mass of each component
  m_c        = rho_cor*(pi*l_c)*(r_cor^2-r_med^2)                # mass of the cortex part of the calamus
  m_r_cor    = (4/3)*rho_cor*(l_r_cor*r_cor^2 - l_r_med*r_med^2) # mass of the cortex part of the rachis
  m_r_med    = (4/3)*rho_med*(r_med^2)*l_r_med                   # mass of the medullary part of the rachis
  mass_outer = (4/3)*rho_cor*r_cor^2*l_r_cor                     # mass as if entire rachis was solid cortex
  mass_inner = (4/3)*rho_cor*r_med^2*l_r_med                     # mass as if hollow part of rachis was solid cortex

  if (r_med > r_cor){
    warning("Incorrect feather medullary radius. Larger than the exterior radius.")
  }

  # --------------------------- Moment of inertia -----------------------------------

  # 1. For each component compute their inertia about their center of mass in the centroidal axes
  # 2. Rotate the axes to be in the feather rachis axis (Vanes only)
  # 3. Parallel axis to the start of feather vane (Vanes only)
  # 4. Sum the vanes and rachis components
  # 5. Rotate the axes to be in the feather calamus axis (Rachis and Vanes only)
  # 6. Parallel axis to the start of feather
  # 7. Sum all feather components
  # 8. Rotate axes so that the feather tip will fall on the z-axis
  # 9. Parallel axes to VRP axes
  # 10. Rotate axes to the VRP
  # ------- Calamus -------
  # 1. Moment of inertia tensors in the calamus
  I_cCG    = calc_inertia_cylhollow(r_cor, r_med, l_c, m_c) # Frame of reference: Feather Calamus | Origin: Calamus CG
  # 6. Adjust so that the origin is at the start of the feather - READY TO BE SUMMED WITH OTHER COMPONENTS
  CG_c1 = c(0,0,0.5*l_c)                                        # Frame of reference: Feather Calamus | Origin: Start of Feather
  I_c1  = parallelaxis(I_cCG,-CG_c1,m_c, "CG")                  # Frame of reference: Feather Calamus | Origin: Start of Feather

  # ------- Rachis -------
  # - Medullary - solid square pyramid
  # Moment of inertia tensor - inner rachis medullary component
  # CAUTION: this is about the center of base
  I_r_med_base  = calc_inertia_pyrasolid(2*r_med, l_r_med, m_r_med)        # Frame of reference: Feather Rachis | Origin: Start of the vane (center)
  # - Cortex - hollow square pyramid
  # Moment of inertia tensor - as if solid inner and outer cortex components
  # CAUTION: this is about the center of base
  I_r_out_cor = calc_inertia_pyrasolid(2*r_cor, l_r_cor, mass_outer)         # Frame of reference: Feather Rachis | Origin: Start of the vane
  I_r_in_cor  = calc_inertia_pyrasolid(2*r_med, l_r_med, mass_inner)         # Frame of reference: Feather Rachis | Origin: Start of the vane
  I_r_cor_base = I_r_out_cor - I_r_in_cor # hollow square pyramid            # Frame of reference: Feather Rachis | Origin: Start of the vane

  CG_pyrasquare = c(0,0,0.25) # multiply by length
  CG_r_out_cor  = CG_pyrasquare*l_r_cor                                      # Frame of reference: Feather Rachis | Origin: Start of the vane
  CG_r_in_cor   = CG_pyrasquare*l_r_med                                      # Frame of reference: Feather Rachis | Origin: Start of the vane
  # Calculate the CG of hollow cortex pyramid
  CG_r_cor1     = (mass_outer*CG_r_out_cor - mass_inner*CG_r_in_cor)/m_r_cor # Frame of reference: Feather Rachis | Origin: Start of the vane
  # Calculate the CG of solid medullary pyramid
  CG_r_med1     = CG_pyrasquare*l_r_med                                      # Frame of reference: Feather Rachis | Origin: Start of the vane
  # Calculate the CG of the entire rachis
  CG_r1         = ((CG_r_cor1*m_r_cor)+(CG_r_med1*m_r_med))/(m_r_cor+m_r_med)# Frame of reference: Feather Rachis | Origin: Start of the vane
  # -------------------

  # ------- Vane -------
  # 1. Moment of inertia for each vane about their centroidal axes
  I_vd1  = calc_inertia_platerect(w_vd, l_r_cor, m_vd)    # Frame of reference: Feather Distal Vane | Origin: Distal Vane CG
  CG_vd1 = c(0,-0.5*w_vd,0.5*l_r_cor)                     # Frame of reference: Feather Distal Vane | Origin: Start of the vane (distal edge)

  I_vp1  = calc_inertia_platerect(w_vp, l_r_cor, m_vp)    # Frame of reference: Feather Proximal Vane | Origin: Proximal Vane CG
  CG_vp1 = c(0,0.5*w_vp,0.5*l_r_cor)                      # Frame of reference: Feather Proximal Vane | Origin: Start of the vane (proximal edge)

  # 2. Rotate positive angle is a ccw rotation about x (normal to the feather plane) ** all angles should be negative for bird feathers
  rot_vd = rotx(-pracma::atand(r_cor/l_r_cor))
  rot_vp = rotx(pracma::atand(r_cor/l_r_cor))
  I_vd2  = rot_vd %*% I_vd1 %*% t(rot_vd)                 # Frame of reference: Feather Rachis | Origin: Distal Vane CG
  I_vp2  = rot_vp %*% I_vp1 %*% t(rot_vp)                 # Frame of reference: Feather Rachis | Origin: Proximal Vane CG
  CG_vd2 = rot_vd %*% CG_vd1                              # Frame of reference: Feather Rachis | Origin: Start of the vane (distal edge)
  CG_vp2 = rot_vp %*% CG_vp1                              # Frame of reference: Feather Rachis | Origin: Start of the vane (proximal edge)

  # 3. adjust I and CG to the start of the vane
  CG_vd3 =  CG_vd2 - c(0,r_cor,0)                         # Frame of reference: Feather Rachis | Origin: Start of the vane
  CG_vp3 =  CG_vp2 + c(0,r_cor,0)                         # Frame of reference: Feather Rachis | Origin: Start of the vane
  I_vd3  = parallelaxis(I_vd2,-CG_vd3,m_vd, "CG")         # Frame of reference: Feather Rachis | Origin: Start of the vane
  I_vp3  = parallelaxis(I_vp2,-CG_vp3,m_vp, "CG")         # Frame of reference: Feather Rachis | Origin: Start of the vane

  # ----------- Rotate all Rachis and Vane components relative to the calamus -----------------
  # 4. sum the rachis and vane components that are in Frame of reference: Feather Rachis | Origin: Start of the vane
  m_r    = m_r_cor+m_r_med+m_vd+m_vp
  I_vr1  = I_vp3 + I_vd3 + I_r_med_base + I_r_cor_base
  CG_vr1 = ((CG_r_cor1*m_r_cor)+(CG_r_med1*m_r_med)+(CG_vd3*m_vd)+(CG_vp3*m_vp))/(m_r) # Frame of reference: Feather Rachis | Origin: Start of the vane

  # 5. Rotate the frame of reference from the rachis to the calamus
  I_vr2  = rotx(angle) %*% I_vr1 %*% t(rotx(angle))              # Frame of reference: Feather Calamus | Origin: Start of the vane
  CG_vr2 = rotx(angle) %*% CG_vr1                                # Frame of reference: Feather Calamus | Origin: Start of the vane

  # 6. Adjust the origin to first the CG of the entire rachis and then to the start of the feather
  I_vrCG  = parallelaxis(I_vr2,CG_vr2,m_r,"A")                   # Frame of reference: Feather Calamus | Origin: Rachis & Vane CG
  CG_vr3  = CG_vr2 + c(0,0,l_c)                                  # Frame of reference: Feather Calamus | Origin: Start of Feather
  I_vr3   = parallelaxis(I_vrCG,-CG_vr3,m_r,"CG")                # Frame of reference: Feather Calamus | Origin: Start of Feather

  # 7. Sum all feather components - must be in Frame of reference: Feather Calamus | Origin: Start of Feather -----------------
  I_1  = I_c1 + I_vr3
  CG_1 = ((CG_c1*m_c)+(CG_vr3*m_r))/m_f

  # 8. Rotate the axes about the x-axis with the origin at the start of the feather so that the length of the calamus is no longer directly along the z-axis.
  full_rot = rotx(pracma::atand(l_r_cor*abs(pracma::sind(angle))/(l_c + l_r_cor*abs(pracma::cosd(angle)))))
  I_2     = full_rot %*% I_1 %*% t(full_rot)                  # Frame of reference: Feather start to tip | Origin: Start of Feather
  CG_2    = full_rot %*% CG_1                                 # Frame of reference: Feather start to tip | Origin: Start of Feather

  I_fCG   = parallelaxis(I_2,CG_2,m_f,"A")                    # Frame of reference: Feather start to tip | Origin: Feather CG

  mass_prop    = list() # pre-define
  mass_prop$I  = I_fCG
  mass_prop$CG = CG_2
  mass_prop$m  = m_f

  return(mass_prop)
}

structural2VRP_feat <- function(m_f, I_fCG, CG_start, start, end, normal){
  # ------------------------------- Adjust axis -------------------------------------
  # first find the frame where z points towards the tip then rotate to frame where z axis points straight along the calamus
  z_axis = end - start                                   # Frame of reference: VRP | Origin: VRP
  x_axis = normal                                        # Frame of reference: VRP | Origin: VRP

  # calculate the rotation matrix between VRP frame of reference and the feather start to tip
  VRP2object = calc_rot(z_axis,x_axis)

  # To properly use parallel axis first, return the origin to the center of gravity of the full feather
  # 9. Need to return origin to VRP before rotating to the final axis system
  off     = VRP2object %*% start                         # Frame of reference: Feather start to tip | Origin: VRP
  CG_3    = CG_start + off                               # Frame of reference: Feather start to tip | Origin: VRP
  I_3     = parallelaxis(I_fCG,-CG_3,m_f,"CG")           # Frame of reference: Feather start to tip | Origin: VRP

  # 10. Rotate the FOR to the VRP axes
  mass_prop    = list() # pre-define
  mass_prop$I  = t(VRP2object) %*% I_3 %*% VRP2object    # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_3                  # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m_f

  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### ------------------------ Mass properties of Neck -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Neck mass properties
#'
#' Calculate the moment of inertia of a neck modeled as a solid cylinder
#'
#' @param m Mass of muscle (kg)
#' @param rho Density of muscle (kg/m^3)
#' @param start a 1x3 vector (x,y,z) representing the 3D point where neck starts. Frame of reference: VRP | Origin: VRP
#' @param end a 1x3 vector (x,y,z) representing the 3D point where neck ends. Frame of reference: VRP | Origin: VRP
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a neck modeled as a solid cylinder}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a neck modeled as a solid cylinder}
#'\item{m}{a double that returns the neck mass}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export

massprop_neck <- function(m,r,l,start,end){

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start
  temp_vec = c(0,-1,0) # In this case this is not arbitrary and defines the y axis
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))
  VRP2object = calc_rot(z_axis,x_axis)

  # -------------------------- Moment of inertia --------------------------------

  I_n = calc_inertia_cylsolid(r, l, m)                    # Frame of reference: Neck | Origin: Neck CG

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off = VRP2object %*% start                               # Frame of reference: Neck | Origin: VRP

  # determine the offset vector for each component
  I_n_off  = c(0,0,0.5*l)+ off                             # Frame of reference: Neck | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_n_vrp   = parallelaxis(I_n,-I_n_off,m,"CG")            # Frame of reference: Neck | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_n_vrp %*% VRP2object  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(start + end)                         # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m
  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### ------------------------ Mass properties of Head -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Head mass properties
#'
#' Calculate the moment of inertia of a head modeled as a solid cone
#'
#' @param m Mass of head (kg)
#' @param r Maximum head radius (m)
#' @param l Maximum head length (m)
#' @param start a 1x3 vector (x,y,z) representing the 3D point where head starts. Frame of reference: VRP | Origin: VRP
#' @param end a 1x3 vector (x,y,z) representing the 3D point where head ends. Frame of reference: VRP | Origin: VRP
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a head modeled as a solid cone}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a head modeled as a solid cone}
#'\item{m}{a double that returns the head mass}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export

massprop_head <- function(m,r,l,start,end){

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start
  temp_vec = c(0,-1,0) # In this case this is not arbitrary and defines the y axis
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))
  VRP2object = calc_rot(z_axis,x_axis)

  # -------------------------- Moment of inertia --------------------------------

  I_h = calc_inertia_conesolid(r, l, m)                    # Frame of reference: Head | Origin: Head CG

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off = VRP2object %*% start                               # Frame of reference: Head | Origin: VRP

  # determine the offset vector for each component
  I_h_off  = c(0,0,0.25*l)+ off                             # Frame of reference: Head | Origin: VRP

  # need to adjust the moment of inertia tensor
  I_h_vrp   = parallelaxis(I_h,-I_h_off,m,"CG")            # Frame of reference: Head | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_h_vrp %*% VRP2object  # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(start + end)                         # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m

  return(mass_prop)
}

# ---------------------------------------------------------------------------------------
##### ------------------------- Mass properties of Tail -------------------------- #####
# ---------------------------------------------------------------------------------------


#' Head mass properties
#'
#' Calculate the moment of inertia of a tail modeled as a solid cone
#'
#' @param m Mass of muscle (kg)
#' @param l Length of the tail (m)
#' @param w Width of the tail (m)
#' @param start a 1x3 vector (x,y,z) representing the 3D point where head starts. Frame of reference: VRP | Origin: VRP
#' @param end a 1x3 vector (x,y,z) representing the 3D point where head ends. Frame of reference: VRP | Origin: VRP
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of a head modeled as a solid cone}
#' \item{CG}{a 1x3 vector representing the center of gravity position of a head modeled as a solid cone}
#'\item{m}{a double that returns the head mass}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export

massprop_tail <- function(m,l_tail,w_tail,l_torso,start,end){

  # ------------------------------- Adjust axis -------------------------------------
  z_axis     = end-start
  temp_vec   = c(0,-1,0) # must be this value to ensure that the x axis in the tail FOR is the VRP z axis
  x_axis     = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2"))
  VRP2object = calc_rot(z_axis,x_axis)

  # -------------------------- Moment of inertia --------------------------------

  I_t1 = calc_inertia_platerect(w_tail, l_tail, m)         # Frame of reference: Tail | Origin: Tail CG
  CG_t = c(0,0,l_torso + 0.5*l_tail)                       # Frame of reference: Tail | Origin: VRP
  I_t2 = parallelaxis(I_t1,-CG_t,m,"CG")                   # Frame of reference: Tail | Origin: VRP

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_t2 %*% VRP2object     # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_t                    # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m

  return(mass_prop)
}

# ---------------------------------------------------------------------------------------
##### ------------------------- Mass properties of Torso -------------------------- #####
# ---------------------------------------------------------------------------------------

#' Torso and leg mass properties
#'
#' Calculate the moment of inertia of a head modeled as a solid cone
#'
#' @param m_true Mass of the torso and legs - no tail (kg)
#' @param m_legs Mass of the legs only (kg)
#' @param w_max Maximum width of the body (m)
#' @param h_max Maximum height of the body (m)
#' @param l_bmax x location of the maximum width of the body (m)
#' @param w_leg width of the body at leg insertion (m)
#' @param w_leg height of the body at leg insertion (m)
#' @param l_leg x location of the leg insertion point (m)
#' @param l_tot length of body from clavicle to end of the tail (m)
#' @param CG_true_x x location of the CG for the torso and legs, origin is at the VRP  (m)
#' @param CG_true_z z location of the CG for the torso and legs, origin is at the VRP (m)
#' @param start a 1x3 vector (x,y,z) representing the 3D point where torso starts. Frame of reference: VRP | Origin: VRP
#' @param end a 1x3 vector (x,y,z) representing the 3D point where tail ends. Frame of reference: VRP | Origin: VRP
#'
#' @author Christina Harvey
#'
#' @return This function returns a list that includes:
#' \itemize{
#' \item{I}{a 3x3 matrix representing the moment of inertia tensor of the torso and leg composite body}
#' \item{CG}{a 1x3 vector representing the center of gravity position of the torso and leg composite body}
#'  \item{m}{a double that returns the mass of the torso and leg composite body}
#' }
#'
#' @section Warning:
#' Parallel axis theorem does not apply between two arbitrary points. One point must be the object's center of gravity.
#'
#' @export

massprop_torso <- function(m_true, m_legs, w_max, h_max, l_bmax, w_leg, h_leg, l_leg, l_tot, CG_true_x, CG_true_z, start, end){

  # ------------------------------- Adjust axis -------------------------------------
  z_axis = end-start
  temp_vec = c(0,-1,0) # In this case this is not arbitrary and defines the y axis
  x_axis = pracma::cross(z_axis,temp_vec/norm(temp_vec, type = "2")) # z points along the VRP positive z axis
  VRP2object = calc_rot(z_axis,x_axis)

  # -------------------------- Moment of inertia --------------------------------
  # pre-define info about the partial elliptic cone
  h_leg  = w_leg*h_max/w_max
  l_par  = l_leg-l_bmax
  l_full = l_par*w_max/(w_max-w_leg)
  l_end  = l_tot - l_leg

  # --------------- Legs - point mass -------------------------
  # placed at the bottom of the bird body
  CG_leg_right = c(-0.5*h_leg, 0.5*w_leg, l_leg)        # Frame of reference: Torso | Origin: VRP
  CG_leg_left  = c(-0.5*h_leg,-0.5*w_leg, l_leg)        # Frame of reference: Torso | Origin: VRP
  leg_right = massprop_pm(0.5*m_legs, CG_leg_right)     # Frame of reference: Torso | Origin: VRP
  leg_left  = massprop_pm(0.5*m_legs, CG_leg_left)      # Frame of reference: Torso | Origin: VRP

  # -------------- Calculate the body volume to have an estimate for the density --------------
  # - 1. volume and CGx of the hemiellipsoid -
  v_ell  = (1/6)*(pi*w_max*h_max*l_bmax)
  CG_ell = (5/8)*l_bmax                                   # Frame of reference: Torso | Origin: VRP

  # - 2. volume of the partial elliptical cone -
  # volume as if the interior cone went to full length
  v_full  = (1/12)*pi*w_max*h_max*l_full
  #volume of the ghost part of the cone
  v_cut   = (1/12)*pi*w_leg*h_leg*(l_full-l_par)
  v_par   = v_full-v_cut
  CG_full = l_bmax + 0.25*l_full                          # Frame of reference: Torso | Origin: VRP
  CG_cut  = l_leg + 0.25*(l_full-l_par)                   # Frame of reference: Torso | Origin: VRP
  CG_par  = (1/v_par)*(v_full*CG_full - v_cut*CG_cut)     # Frame of reference: Torso | Origin: VRP

  # - 3. volume and CGx of the full elliptical cone -
  # volume of the full end cone
  v_end  = (1/12)*(pi*w_leg*h_leg*l_end)
  CG_end  = l_leg+0.25*l_end                              # Frame of reference: Torso | Origin: VRP

  v_body  = v_ell + v_par + v_end
  m_body  = (m_true-m_legs)
  rho_avg = m_body/v_body

  #adjust CG of torso to remove the effects of the legs
  CG_body_z = (CG_true_z*m_true - CG_leg_left[1]*0.5*m_legs -CG_leg_right[1]*0.5*m_legs)/m_body
  CG_body_x = (CG_true_x*m_true - CG_leg_left[3]*0.5*m_legs -CG_leg_right[3]*0.5*m_legs)/m_body # Yes this is supposed to in index 3 - in the torso FOR

  # ----------- Estimate the density in each section of the body -----------------
  A    = matrix(c(v_par,v_par*CG_par,v_end,v_end*CG_end),2,2)
  test = pracma::lsqnonlin(density_optimizer, rho_avg, options=list(tolx=1e-12, tolg=1e-12), A, m_body, v_ell, CG_body_x, CG_ell, rho_avg)
  rho_ell  = test$x
  b        = c((m_body-(rho_ell*v_ell)),((m_body*CG_body_x)-(rho_ell*v_ell*CG_ell)))
  output   = solve(A,b)
  rho_par  = output[1]
  rho_back = output[2]

  # -------------- Hemiellipsoid -------------------
  # mass of the front
  m_ell   = rho_ell*v_ell

  # MOI about the base
  I_ell1  = calc_inertia_ellipse(h_max/2, w_max/2, l_bmax, m_ell)   # Frame of reference: Torso | Origin: Hemiellipsoid base
  # Define center of gravity wrt to I origin
  CG_ell1 = c(0,0,-(3/8)*l_bmax)                                    # Frame of reference: Torso | Origin: Hemiellipsoid base

  # Adjust the moment of inertia to be about the center of gravity of the hemiellipsoid
  I_ellCG = parallelaxis(I_ell1,CG_ell1,m_ell,"A")                  # Frame of reference: Torso | Origin: Hemiellipsoid CG
  # Define center of gravity wrt to VRP origin
  # NOTE: we need to offset the body from the x axis to ensure that the VRP z-location of the CG is correct in the Torso frame of ref this is the positive x direction
  CG_ell2 = c(CG_body_z,0,(5/8)*l_bmax)                             # Frame of reference: Torso | Origin: VRP
  # Adjust the MOI to be about the VRP
  I_ell_vrp = parallelaxis(I_ellCG,-CG_ell2,m_ell,"CG")             # Frame of reference: Torso | Origin: VRP

  # -------------- Partial elliptic cone -------------------
  # mass of the partial cone
  m_par  = rho_par*v_par

  # CG as if the interior cone went to full length
  CG_full = c(0,0,0.25*l_full)                                      # Frame of reference: Torso | Origin: Hemiellipsoid base
  # CG of the ghost part of the cone
  CG_cut  = c(0,0,l_par+0.25*(l_full-l_par))                        # Frame of reference: Torso | Origin: Hemiellipsoid base
  # CG of the partial cone
  CG_par1 = (1/v_par)*(v_full*CG_full - v_cut*CG_cut)               # Frame of reference: Torso | Origin: Hemiellipsoid base

  # MOI of the partial cone about the wide base
  I_par1  = calc_inertia_ellcone(w_max, h_max, l_par, l_full, rho_par)  # Frame of reference: Torso | Origin: Hemiellipsoid base
  # Adjust the moment of inertia to be about the center of gravity of the partial cone
  I_parCG = parallelaxis(I_par1,CG_par1,m_par,"A")                      # Frame of reference: Torso | Origin: Partial Cone CG
  # Define CG wrt to VRP origin
  # NOTE: we need to offset the body from the x axis to ensure that the VRP z-location of the CG is correct in the Torso frame of ref this is the positive x direction
  CG_par2 = CG_par1 + c(CG_body_z,0,l_bmax)                             # Frame of reference: Torso | Origin: VRP
  # Adjust the MOI to be about the VRP
  I_par_vrp = parallelaxis(I_parCG,-CG_par2,m_par,"CG")                 # Frame of reference: Torso | Origin: VRP

  # -------------- Full elliptic cone -------------------
  # mass of the end cone
  m_end  = rho_back*v_end

  # CG of the end cone
  CG_end1 = c(0,0,0.25*l_end)                                              # Frame of reference: Torso | Origin: Base of the end cone
  # MOI of the end cone about the wide base
  I_end1  = calc_inertia_ellcone(w_leg, h_leg, l_end, l_end, rho_back)     # Frame of reference: Torso | Origin: Base of the end cone

  # Adjust the moment of inertia to be about the center of gravity of the partial cone
  I_endCG = parallelaxis(I_end1,CG_end1,m_end,"A")                         # Frame of reference: Torso | Origin: End Cone CG
  # Define CG wrt to VRP origin
  # NOTE: we need to offset the body from the x axis to ensure that the VRP z-location of the CG is correct in the Torso frame of ref this is the positive x direction
  CG_end2 = CG_end1 + c(CG_body_z,0,(l_bmax + l_par))                      # Frame of reference: Torso | Origin: VRP
  # Adjust the MOI to be about the VRP
  I_end_vrp = parallelaxis(I_endCG,-CG_end2,m_end,"CG")                    # Frame of reference: Torso | Origin: VRP

  # Check that the mass matches the expected value
  if(abs(m_end+m_par+m_ell-m_body)>0.001){
    warning("The mass of the predicted body shape does not match the expected value")
  }

  # -------------- Sum data and adjust axes -------------------
  I_torso_vrp  = I_ell_vrp + I_par_vrp + I_end_vrp + leg_right$I + leg_left$I
  CG_torso_vrp = (1/m_true)*(m_ell*CG_ell2 + m_par*CG_par2 + m_end*CG_end2 + 0.5*m_legs*CG_leg_left + 0.5*m_legs*CG_leg_right)

  # we need to offset the body from the x axis to ensure that the VRP z-location of the CG is correct in the Torso frame of ref this is the positive x direction

  mass_prop = list() # pre-define
  # Adjust frame to VRP axes
  mass_prop$I  = t(VRP2object) %*% I_torso_vrp %*% VRP2object               # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = t(VRP2object) %*% CG_torso_vrp                             # Frame of reference: VRP | Origin: VRP
  mass_prop$m  = m_true

  # Check that the CG matches the expected value
  if(abs(abs(mass_prop$CG[1])-CG_true_x)>0.01){
    warning("The center of gravity the predicted body shape does not match the expected value.")
  }
  return(mass_prop)
}


# ---------------------------------------------------------------------------------------
##### -------------------- Body component density optimizer ----------------------- #####
# ---------------------------------------------------------------------------------------

density_optimizer <- function(x,A,m_body,v_ell,CG_body_x,CG_ell,rho_avg){
  rho = x
  b      = c((m_body-(rho*v_ell)),
             ((m_body*CG_body_x)-(rho*v_ell*CG_ell)))
  output = solve(A,b)

  err = ((rho-rho_avg)^2 + (output[1]-rho_avg)^2 + (output[2]-rho_avg)^2 + (rho-output[1])^2 + (rho-output[2])^2)

  if(output[1] < 0 | output[2] < 0 | rho < 0){
    err = err + rho_avg
  }

  return(err)
}
