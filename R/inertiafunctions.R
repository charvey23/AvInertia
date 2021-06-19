

# ------------------------------------------------------------------------------
# ------------------------- Moment of inertia tensors Properties  --------------
# ------------------------------------------------------------------------------

# ---------------------- Mass Properties - Solid cylinder ----------------------
#' Moment of inertia tensor of a solid cylinder
#'
#' Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#'
#' @param r radius of the cylinder (m)
#' @param h height of the cylinder (m)
#' @param m mass of the cylinder (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of
#' a solid cylinder about its center of gravity with z oriented
#' through it's major axis
#'
#' @export
#'
calc_inertia_cylsolid <- function(r, h, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (1/12)*(3*r^2+h^2) # Ixx
  I[2,2] = (1/12)*(3*r^2+h^2) # Iyy
  I[3,3] = (0.5)*(r^2)        # Izz
  I = m*I

  return(I)
}
# ------------------------------------------------------------------------------
# ---------------------- Mass Properties - Hollow cylinder ---------------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a hollow cylinder
#'
#' Reference:https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#'
#' @param r_out outer radius of the cylinder (m)
#' @param r_in inner radius of the cylinder (m)
#' @param h height of the cylinder (m)
#' @param m mass of the cylinder (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a hollow
#' cylinder about its center of gravity with z oriented through it's major axis
#'
#' @export
#'
calc_inertia_cylhollow <- function(r_out, r_in, h, m){

  temp_r = r_out^2 + r_in^2

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (1/12)*(3*temp_r+h^2) # Ixx
  I[2,2] = (1/12)*(3*temp_r+h^2) # Iyy
  I[3,3] = (0.5)*(temp_r)        # Izz
  I = m*I

  return(I)
}

# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Solid Ellipse   -----------------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of solid ellipse CG or a half ellipse centered on
#' the base
#'
#' Reference:https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#'
#' @param a half the height along the x direction (m)
#' @param b half the width along the y direction  (m)
#' @param c half the length along the z direction (m)
#' @param m mass of the ellipse (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid
#' ellipse about its center of gravity with the major axes aligned.
#'
#' @section
#' CAUTION: Origin is at the center of gravity for a full ellipse or at the
#' center of the base if modelling a half ellipse.
#'
#' @export
#'
calc_inertia_ellipse <- function(a, b, c, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (b^2 + c^2) # Ixx
  I[2,2] = (a^2 + c^2) # Iyy
  I[3,3] = (a^2 + b^2) # Izz
  I = (1/5)*m*I

  return(I)
}
# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Solid Circular Cone   -----------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a solid circular cone pyramid
#'
#' All outputs are based on an origin at the centered point on the base
#' Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#'
#' @param r radius of the cone base (m)
#' @param h height of the cone (m)
#' @param m mass of the cone (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid
#' circular cone about
#' its center of gravity with z oriented through it's major axis.
#'
#' @section
#' CAUTION: Origin of the output tensor is NOT at the center of gravity but
#' at the center of the base.
#'
#' @export
#'
calc_inertia_conesolid <- function(r, h, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (1.5*r^2+h^2) # Ixx
  I[2,2] = (1.5*r^2+h^2) # Iyy
  I[3,3] = 3*(r^2)       # Izz
  I = (1/10)*m*I

  return(I)
}
# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Solid Elliptical Cone  ----------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a solid elliptical cone - end of purple notebook
#' derivation verified in green
#'
#' @param A   half height of the wide base (m)
#' @param B   half width of the wide base (m)
#' @param l   length of the cone (m)
#' @param m   mass of the cone (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid
#' elliptical cone about the center of the wider base with z oriented
#' towards the end.
#'
#' @section
#' CAUTION: Origin of the output tensor is NOT at the center of gravity
#' but at the center of the base.
#' @export
#'
calc_inertia_ellcone <- function(A, B, l, m){

  I = matrix(0, nrow = 3, ncol = 3)
  tmpx2 = (3/20)*A^2
  tmpy2 = (3/20)*B^2
  tmpz2  = (1/10)*l^2
  I[1,1] = tmpy2 + tmpz2 # Ixx
  I[2,2] = tmpx2 + tmpz2 # Iyy
  I[3,3] = tmpx2 + tmpy2 # Izz
  I = m*I

  return(I)
}
# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Solid Elliptical Cylinder  ------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a solid elliptical cylinder
#' Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#'
#' @param a   half height of the base - oriented along the x axis in torso FOR
#' (z axis in full bird FOR) (m)
#' @param b   half width of the base - oriented along the y axis in torso FOR
#' (y axis in full bird FOR)  (m)
#' @param l   length of the cylinder - oriented along the z axis in torso FOR
#' (x axis in full bird FOR)  (m)
#' @param m mass of the cylinder (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid
#' elliptical cylinder about it's center of gravity
#'
#' @export
#'
calc_inertia_ellcyl <- function(a, b, l, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (1/4)*b^2 + (1/12)*l^2 # Ixx
  I[2,2] = (1/4)*a^2 + (1/12)*l^2 # Iyy
  I[3,3] = (1/4)*(a^2+b^2)  # Izz
  I = m*I

  return(I)
}

# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Solid Square Pyramid  -----------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a solid square pyramid
#' Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#' All outputs are based on an origin at the centered point on the base
#'
#' @param w entire width of one side of the pyramid base (m)
#' @param h entire height of the pyramid (m)
#' @param m mass of the pyramid (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid
#' square pyramid about its center of gravity with z oriented through it's
#' major axis.
#'
#' @section
#' CAUTION: Origin is NOT at the center of gravity but at the center of the base.
#'
#' @export
#'
calc_inertia_pyrasolid <- function(w, h, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (w^2+2*h^2) # Ixx
  I[2,2] = (w^2+2*h^2) # Iyy
  I[3,3] = 2*(w^2)     # Izz
  I = (1/20)*m*I

  return(I)
}

# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Flat Rectangular Plate  ---------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a flat rectangular plate - assumes thickness is
#' approximately zero
#' Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
#' @param w full width of one side of the plate - short edge (m)
#' @param h height of the plate - long edge (m)
#' @param m mass of the plate (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a flat
#' plate about its center of gravity with z oriented parallel with it's long
#' edge (h) and y along its short edge (w)
#'
#' @export
#'
calc_inertia_platerect <- function(w, h, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (w^2+h^2) # Ixx
  I[2,2] = (h^2)     # Iyy
  I[3,3] = (w^2)     # Izz
  I = (1/12)*m*I

  return(I)
}
# ------------------------------------------------------------------------------
# -------------------- Mass Properties - Flat Triangular Plate  ----------------
# ------------------------------------------------------------------------------

#' Moment of inertia tensor of a flat triangular plate
#'
#' Reference: https://apps.dtic.mil/dtic/tr/fulltext/u2/a183444.pdf
#' page 4 eqns 2.16-2.20
#'
#' @param pts a matrix of the three 3D points that define a point on
#' the triangular plate.
#' Frame of reference: Muscle | Origin: VRP
#' each point should be a different row as follows:
#' pt1x, pt1y, pt1z
#' pt2x, pt1y, pt2z
#' pt3x, pt3y, pt3z
#'
#' @param a  area of the triangular plate (m)
#' @param rho density of the material (kg/m^3)
#' @param t  thickness of the plate (m)
#' @param desired_prop a string containing either "I" or "CG" depending on the
#' desired output
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a flat
#' triangular plate about
#' its center of gravity. Z axis defined as the normal to the input points.
#'
#' @section Warning:
#' The input points should be defined in a counterclockwise direction around the
#' plate in the triangular plate frame of reference
#' i.e. all z components should be equal
#'
#' @export
#'
calc_inertia_platetri <- function(pts, A, rho, t, desired_prop){
  # need to add the first point to allow circular calculation
  pts = rbind(pts,pts[1,])

  # ------------ Center of Gravity Calculation -----------------
  # Returns the centroid times the area. Ref: page 2
  x1  = 0
  y1  = 0
  for (i in 1:3){
    x1   = x1 + (1/(6*A))*(((pts[i,1]*pts[i+1,2])-
                              (pts[i+1,1]*pts[i,2]))*(pts[i,1]+
                                                        pts[i+1,1])); # eqn 2.16
    y1   = y1 + (1/(6*A))*(((pts[i,1]*pts[i+1,2])-
                              (pts[i+1,1]*pts[i,2]))*(pts[i,2]+
                                                        pts[i+1,2])); # eqn 2.17
  }
  # Frame of reference: Incoming points | Origin: Incoming points
  CG = c(x1, y1, pts[1,3]);

  # ------------ Moment of Inertia Calculation -----------------
  # Returns the second moment of area
  if(desired_prop == "I"){
    # adjust the incoming points to be so that the origin of the pts is
    # on the plate's center of gravity
    adj_pts = pts # define the new matrix
    for (i in 1:4){
      adj_pts[i,] = pts[i,] - CG
    }

    # calculate the moment of inertia of the plate
    x2  = 0;
    y2  = 0;
    Ixy = 0;
    for (i in 1:3){
      temp = ((adj_pts[i,1]*adj_pts[i+1,2])-(adj_pts[i+1,1]*adj_pts[i,2]))
      x2   = x2  + (temp*(adj_pts[i,1]^2+
                            (adj_pts[i,1]*adj_pts[i+1,1])+
                            adj_pts[i+1,1]^2)); # eqn 2.18
      y2   = y2  + (temp*(adj_pts[i,2]^2+
                            (adj_pts[i,2]*adj_pts[i+1,2])+
                            adj_pts[i+1,2]^2)); # eqn 2.20
      Ixy  = Ixy + (0.5)*(temp*((2*adj_pts[i,1]*adj_pts[i,2]) +
                                  adj_pts[i,1]*adj_pts[i+1,2] +
                                  adj_pts[i,2]*adj_pts[i+1,1] +
                                  (2*adj_pts[i+1,1]*adj_pts[i+1,2]))); # eqn 2.19
    }


    I = matrix(0, nrow = 3, ncol = 3) # define matrix
    I[1,1] = y2       # Ixx - (y^2)dA
    I[2,2] = x2       # Iyy - (x^2)dA
    I[3,3] = (x2+y2)  # Izz - (x^2 + y^2)dA
    I[1,2] = -Ixy     # Ixy - (xy)dA
    I[2,1] = -Ixy     # Iyx - (xy)dA

    I = (1/12)*rho*t*I

    return(I)
  }


  if(desired_prop == "CG"){
    return(CG)
  }
}

# ------------------------------------------------------------------------------
# ------------------- Mass Properties - Flight Feather -------------------------
# ------------------------------------------------------------------------------

#' Compute the inertia of the individual feathers
#'
#' @param dat_mat Dataframe related to the current species input as a
#' dataframe with the following structure:
#' \itemize{
#' \item{material}{Material information. Must include the following:
#' "Cortex", "Medullary"}
#' \item{density}{Density of each material (kg/m^3)}
#' }
#' @param dat_feat_curr Dataframe related to the current bird wing feathers
#' input as a dataframe with the following structure:
#' \itemize{
#' \item{feather}{Feather ID code. Must be in standard format i.e.
#' 1st primary is "P1", third secondary is "S3", etc.
#' Alula feathers should be grouped and named "alula".}
#' \item{m_f}{Mass of feather in the same row as the
#' appropriate feather ID code (kg)}
#' \item{l_cal}{Length of calamus in the same row as the
#' appropriate feather ID code (m)}
#' \item{l_vane}{Length of rachis/vane in the same row as the
#' appropriate feather ID code (m)}
#' \item{w_cal}{Width (diameter) of calamus in the same row as the
#' appropriate feather ID code (m)}
#' \item{w_vp}{Width of proximal vane (average value) in the same row as the
#' appropriate feather ID code (m)}
#' \item{w_vd}{Width of distal vane (average value)  in the same row as the
#' appropriate feather ID code (m)}
#' \item{vane_angle}{Interior angle between the rachis and calamus  (degrees)}
#' }
#' NOTE: Alula feathers will be treated as point mass so only the mass of the
#' feathers is required. Other columns can be left blank.
#' @param dat_bird_curr Dataframe related to the current bird wing that must
#' include the following columns:
#' \itemize{
#' \item{barb_radius}{Radius of feather barb  for current species (m)}
#' \item{barb_distance}{Distance between feather barbs for current species (m)}
#' }
#'
#' @return
#' @export
#'

compute_feat_inertia <- function(dat_mat, dat_feat_curr, dat_bird_curr){

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
  res_pri$I_pri   = array(dim = c(3,3,no_pri))
  res_pri$CG_pri  = array(dim = c(no_pri,3))
  res_pri$m_pri   = array(dim = c(no_pri))
  res_sec     = list()
  res_sec$I_sec   = array(dim = c(3,3,no_sec))
  res_sec$CG_sec  = array(dim = c(no_sec,3))
  res_sec$m_sec   = array(dim = c(no_sec))

  # --------------------------- Primaries --------------------------------------
  #  P1 -> P10
  for (i in 1:no_pri){
    feather_name = paste("P",i,sep = "")
    # subset data to be for this specific feather
    pri_info = subset(dat_feat_curr,feather == feather_name)
    # Calculate MOI and CG
    tmp = massprop_feathers(pri_info$m_f,
                            pri_info$l_cal,
                            pri_info$l_vane,
                            pri_info$w_cal,
                            dat_bird_curr$barb_radius,
                            dat_bird_curr$barb_distance,
                            rho_cor,rho_med,
                            pri_info$w_vp,
                            pri_info$w_vd,
                            pri_info$vane_angle)
    # Save MOI, CG and CG*mass
    res_pri$I_pri[,,i] = tmp$I
    res_pri$CG_pri[i,] = tmp$CG
    res_pri$m_pri[i]   = tmp$m

  }
  #  -----------------------------------Secondaries ----------------------------
  #  S1 -> last secondary
  for (i in 1:no_sec){
    feather_name = paste("S",i,sep = "")
    # subset data to be for this specific feather
    sec_info = subset(dat_feat_curr,feather == feather_name)
    # Calculate MOI and CG
    tmp = massprop_feathers(sec_info$m_f,
                            sec_info$l_cal,
                            sec_info$l_vane,
                            sec_info$w_cal,
                            dat_bird_curr$barb_radius,
                            dat_bird_curr$barb_distance,
                            rho_cor,rho_med,
                            sec_info$w_vp,
                            sec_info$w_vd,
                            sec_info$vane_angle)
    # Save MOI, CG and CG*mass
    res_sec$I_sec[,,i]  = tmp$I
    res_sec$CG_sec[i,]  = tmp$CG
    res_sec$m_sec[i]    = tmp$m
  }

  return(c(res_pri,res_sec))
}
