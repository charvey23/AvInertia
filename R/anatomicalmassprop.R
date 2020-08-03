# Script containing the functions to calculate the moment of inertia and CG of different anatomical components
# all written by Christina Harvey
# last updated: 2020-08-03


#' Calculate the moment of inertia of a bone modelled as a hollow cylinder with two solid end caps
#'
#' @param m full mass of the bone
#' @param l full length of the bone
#' @param r_out outer radius of bone
#' @param r_in inner radius of bone
#' @param rho density of the bone (kg/m^3)
#' @param start 3D point where bone starts. This point must be defined relative to the vehicle reference point (VRP) in the VRP frame of reference
#' @param end 3D point where bone ends. This point must be defined relative to the vehicle reference point (VRP) in the VRP frame of reference
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
  z_axis = bone_end-bone_start;
  temp_vec = c(1,1,1) # arbitrary
  x_axis = cross(z_axis_b,temp_vec/Norm(temp_vec));

  # doesn't matter where the x axis points as long as:
  # 1. we know what it is 2. it's orthogonal to z
  VRP2object = calc_rot(z_axis,x_axis);

  # ---------------------------------------------------------------------------------
  # ----------------------- Calculate moment of inertia -----------------------------
  # ---------------------------------------------------------------------------------

  # --- Calculate the moment of inertia about the CG of the individual components --
  # Frame of reference: Bone | Origin: Bone CG
  I_cyl_b = calc_inertia_cylhollow(r_out, r_in, cyl_l, cyl_m)
  I_cap_b = calc_inertia_cylsolid(r_out, cap_t, cap_m)

  # want to move origin from CG to VRP but need to know where the bone is relative to the VRP origin
  off_bone = VRP2object %*% bone_start  # Frame of reference: Bone | Origin: VRP

  # determine the offset vector for each component with Frame of reference: Bone | Origin: VRP
  I_cyl_off   = c(0,0,0.5*l_bone) + off_bone            # the hollow cylinder is displaced halfway along the bone (in z-axis)
  I_cap1_off  = c(0,0,0.5*cap_t) + off_bone             # Cap 1 edge centered on the beginning of the bone
  I_cap2_off  = c(0,0,(l_bone - (0.5*cap_t)))+ off_bone # Cap 2 edge centered on the end of the bone

  # need to adjust the moment of inertia tensor - Frame of reference: Bone | Origin: VRP
  I_cyl_vrp   = parallelaxis(I_cyl,-I_cyl_off,m_cyl)
  I_cap1_vrp  = parallelaxis(I_cap,-I_cap1_off,m_cap)
  I_cap2_vrp  = parallelaxis(I_cap,-I_cap2_off,m_cap)

  I_boneaxis  = I_cyl_vrp + I_cap1_vrp + I_cap2_vrp      # Frame of reference: Bone | Origin: VRP

  mass_prop = list()
  mass_prop$I  = t(VRP2object) %*% bone_start %*% VRP2object # Frame of reference: VRP | Origin: VRP
  mass_prop$CG = 0.5*(bone_start + bone_end) # Frame of reference: VRP | Origin: VRP

  return(mass_prop)
}
