# Script containing all basic mathematical functions
# all written by Christina Harvey
# last updated: 2020-09-10

# ---------------------- Kroneckerdelta function ---------------------------
#' Kroneckerdelta function
#'
#' @param i  a scalar index (usually row of a matrix)
#' @param j  a scalar index (usually column of a matrix)
#'
#' @author Christina Harvey
#'
#' @return a scalar value. Returns 1 if i and j are equal otherwise returns 0
#'
#' @examples
#'
#'  # should return 1
#'  kronecker_delta(7,7)
#'
#'  # should return 0
#'  kronecker_delta(5,4)
#'
#' @export
#'
kronecker_delta <- function(i,j){

  if (i == j) {
    return(1)
  }
  else {
    return(0)
  }
}

# ---------------------- Rotation about the x axis ---------------------------
#' A 3x3 rotation matrix allowing rotation about the x-axis. Constructed using a cosine rotation matrix where the rotation angle in degrees
#' is measured counterclockwise allowing positive rotation under the right hand rule.
#'
#' @param angle a scalar representing the angle to rotate (degrees)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the rotation about the x-axis by the given angle
#' @export
#'
#' @examples
#'
rotx <- function(angle){
  R = matrix(0, nrow = 3, ncol = 3)
  R[1,1] = 1
  R[2,2] = pracma::cosd(angle)
  R[3,3] = pracma::cosd(angle)
  R[3,2] = pracma::sind(angle)
  R[2,3] = -pracma::sind(angle)
  return(R)
}

# ---------------------- Calculate a unit vector --------------------------
#' Determine the unit vector of any input vector
#'
#' @param vector any vector or array with only one dimension
#'
#' @author Christina Harvey
#'
#' @return the unit vector in the size of the input vector
#'
#' @examples
#'
#' output = calc_univec(c(2,6,5))
#' length_output = pracma::Norm(output) #should equal 1 if a unit vector
#'
#' @export

calc_univec <- function(vector){
  unit_vector = vector/pracma::Norm(vector)
  return(unit_vector)
}


# ---------------------- Calculate a rotation matrices ---------------------------
#' A 3x3 rotation matrix constructed by projecting the new axes onto the original system. Likely results in rotation about all axes.
#'
#' @param z_vector a 1x3 vector representing the direction for the desired z axis of the new frame of reference. Frame of reference: VRP
#' @param x_vector a 1x3 vector representing the direction for the desired x axis of the new frame of reference. Frame of reference: VRP
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the rotation matrix that transforms between VRP frame and object frame
#'
#' @export
#'
calc_rot <- function(z_vector, x_vector){

  #CAUTION: incoming vectors must be in the structural frame at the VRP
  #         x and z vectors must already be orthogonal axes

  # ensure vectors are in unit form
  unit_z_vector = calc_univec(z_vector)
  unit_x_vector = calc_univec(x_vector)

  # cross z with x to get the right-handed axis - verified
  y_vector = pracma::cross(unit_z_vector,unit_x_vector)
  # ensure vectors are in unit form
  unit_y_vector = calc_univec(y_vector)

  # rotation matrix representing the rotated basis
  VRP2object = rbind(unit_x_vector,unit_y_vector,unit_z_vector)
  # strip row names
  rownames(VRP2object)<-NULL

  return(VRP2object)
}

# ---------------------- Calculate a matrix translation with parallel axis theorem  ---------------------------
#' Parallel axis theory
#'
#' Reads in an initial tensor and an offset to compute the transformed tensor.
#' Will be in the same frame of reference as the input tensor.
#'
#' @param I_CG Moment of intertia tensor (3x3) about the center of gravity of the object (kg-m^2)
#' @param offset_vec Distance between the objects CG and the arbitrary pt A. (m) Vector should always point from the CG to the arbitrary point A.
#' @param m Mass of the object (kg)
#' @param cg_a If input I is about the CG enter "CG" or if I is about an arbitrary axis enter "A".
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the transformed moment of inertia tensor after a solid body translation
#'
#' @export
#'

parallelaxis <- function(I, offset_vec, m, cg_a){

  # CAUTION: the parallel axis theorem only works if the I_CG is given about
  #          the object's centroidal axis

  I_new = matrix(0, nrow = 3, ncol = 3) # predefine matrix

  if(cg_a == "CG"){
    sign = 1
  }

  if(cg_a == "A"){
    sign = -1
  }

  for (i in 1:3){
    for (j in 1:3){
      I_new[i,j] = I[i,j] + sign*m*((kronecker_delta(i,j)*pracma::dot(offset_vec,offset_vec)) - (offset_vec[i]*offset_vec[j]))
    }
  }
  return(I_new)
}

# ---------------------------------------------------------------------------------------
##### -------------------------- Feather Orientation ------------------------------ #####
# ---------------------------------------------------------------------------------------

#' Determine the feather orientation
#'
#' Code that returns the orientation of each primary and secondary feather on the wing.
#'
#' @param no_pri a scalar representing the amount of primary feathers
#' @param no_sec a scalar representing the amount of secondary feathers
#' @param Pt1 a matrix of 3 values representing the point on the shoulder joint (m)
#' @param Pt2 a matrix of 3 values representing the point on the elbow joint (m)
#' @param Pt3 a matrix of 3 values representing the point on the wrist joint (m)
#' @param Pt4 a matrix of 3 values representing the point on the end of carpometacarpus (m)
#' @param Pt9 a matrix of 3 values representing the point on tip of most distal primary (m)
#' @param Pt10 a matrix of 3 values representing the point on tip of last primary to model as if on the end of the carpometacarpus (m)
#' @param Pt11 a matrix of 3 values representing the point on tip of most proximal secondary feather (m)
#'
#' @author Christina Harvey
#'
#' @return a list called "feather". This contains three matrices.
#' 1. "loc_start" a matrix defining the 3D point where each feather starts. Rows are the different feathers and columns are x, y, z coordinates respectively.
#' 2. "loc_end" a matrix defining the 3D point where each feather end. Rows are the different feathers and columns are x, y, z coordinates respectively.
#' 3. "normal" a matrix that gives the vector that defines the normal to each feather plane.  Rows are the different feathers and columns are x, y, z vector directions respectively.
#'
#' @export

orient_feather <- function(no_pri,no_sec,Pt1,Pt2,Pt3,Pt4,Pt9,Pt10,Pt11){

  # pre-define variables
  count             = 1
  feather           = list()
  feather$loc_start = matrix(0, nrow = (no_pri+no_sec), ncol = 3)
  feather$loc_end   = matrix(0, nrow = (no_pri+no_sec), ncol = 3)
  feather$normal    = matrix(0, nrow = (no_pri+no_sec), ncol = 3)
  k = 1
  # --- Primaries ---
  # Calculate the start and end of the primaries
  for(i in 1:no_pri){

    if (i < (no_pri-2)) {
      # Note: We don't want last primary to be exact same spot as S1 or P7 to be at same spot as P8
      #       this requires that we skip the first and last position

      # -- Start -- Primaries less than P7 distribute linearly along the carpometacarpus
      feather$loc_start[count,] = Pt3 + (1/8)*(i)*(Pt4-Pt3)
      # -- End -- Distributes linearly between P7 and S1
      feather$loc_end[count,]   = Pt10 + (1/8)*(i)*(Pt9-Pt10)

    } else {
      # -- Start -- Any primary past P7 starts on the end of the carpometacarpus
      feather$loc_start[count,] = Pt4;
      # -- End -- Distributes linearly
      feather$loc_end[count,]   = Pt9 + (1/(no_pri-7))*(k)*(Pt8-Pt9);
      k = k + 1
    }
    count = count + 1
  }

  # --- Secondaries ---
  # Calculate the start and end of the secondaries
  if(Pt11[2]<0){ # given the orientation it is possible that the secondary moves into a positive area
    sec_vec = Pt11-Pt10;                   # vector between S1 and the wing root trailing edge
    t       = (Pt12[2]-Pt10[2])/sec_vec[2]; # proportion along the vector where it intersects y = Pt12y
    sec_end = c((t*sec_vec[1]+Pt10[1]), Pt12[2], (t*sec_vec[3]+Pt10[3])); # Point along the vector at Pt12Y
  }else{
    sec_end = Pt11
  }
  for(i in 1:no_sec){
    # -- Start --
    #equally space the start of the secondaries along the forearm ulna/radius
    feather$loc_start[count,] = Pt3 + (1/(no_sec-1))*(i-1)*(Pt2-Pt3);
    # -- End --
    feather$loc_end[count,]   = Pt10 + (1/(no_sec-1))*(i-1)*(sec_end-Pt10);
    count = count + 1
  }

  # --- Calculate normal for the feather ----

  for(i in 1:(no_pri+no_sec)){
    if(i <= no_pri){
      feather$normal[i,] = pracma::cross((Pt3-Pt4),(feather$loc_end[i,]-Pt4)) # primaries on the carpometacarpus plane
    } else {
      feather$normal[i,] = pracma::cross((Pt2-Pt3),(feather$loc_end[i,]-Pt3)) # secondaries on the ulna/radius plane
    }
  }

  return(feather)
}
