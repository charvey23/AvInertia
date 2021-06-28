# Script containing all basic mathematical functions

# ------------------------------------------------------------------------------
# ---------------------- Kronecker-delta function ------------------------------
# ------------------------------------------------------------------------------

#' Kroneckerdelta function
#'
#' @param i  a scalar index (usually row of a matrix)
#' @param j  a scalar index (usually column of a matrix)
#'
#' @author Christina Harvey
#'
#' @return a scalar value. Returns 1 if i and j are equal otherwise returns 0
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
# ------------------------------------------------------------------------------
# ---------------------- Rotation about the x axis -----------------------------
# ------------------------------------------------------------------------------
#' A 3x3 rotation matrix allowing rotation about the x-axis.
#' Constructed using a cosine rotation matrix where the rotation angle in
#' degrees is measured counterclockwise allowing positive rotation under
#' the right hand rule.
#'
#' @param angle a scalar representing the angle to rotate (degrees)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the rotation about the x-axis by the
#' given angle
#'
#' @export
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

# ------------------------------------------------------------------------------
# ------------------------ Calculate a unit vector -----------------------------
# ------------------------------------------------------------------------------

#' Determine the unit vector of any input vector
#'
#' @param vector any vector or array with only one dimension
#'
#' @author Christina Harvey
#'
#' @return the unit vector in the size of the input vector
#'
#' @export

calc_univec <- function(vector){
  unit_vector = vector/pracma::Norm(vector)
  return(unit_vector)
}

# ------------------------------------------------------------------------------
# ---------------------- Calculate a rotation matrices -------------------------
# ------------------------------------------------------------------------------

#' A 3x3 rotation matrix constructed by projecting the new axes onto the
#' original system. Likely results in rotation about all axes.
#'
#' @param z_vector a 1x3 vector representing the direction for the desired z
#' axis of the new frame of reference. Frame of reference: VRP
#' @param x_vector a 1x3 vector representing the direction for the desired x
#' axis of the new frame of reference. Frame of reference: VRP
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the rotation matrix that transforms between
#' VRP frame and object frame
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
# ------------------------------------------------------------------------------
# ------------------------- Parallel axis theorem  -----------------------------
# ------------------------------------------------------------------------------

#' Parallel axis theory
#'
#' Reads in an initial tensor and an offset to compute the transformed tensor.
#' Will be in the same frame of reference as the input tensor.
#'
#' @param I a 3x3 matrix representing the moment of inertia tensor about the
#' center of gravity of the object (kg-m^2).
#' @param offset_vec a 1x3 vector representing the distance (x,y,z) between
#' the objects CG and the arbitrary pt A (m).
#' Vector should always point from the CG to the arbitrary point A.
#' @param m Mass of the object (kg).
#' @param cg_a If input I is about the CG enter "CG" or if I is about an
#' arbitrary axis enter "A".
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the transformed moment of inertia tensor
#' after a solid body translation defined by the offset vector.
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
      I_new[i,j] = I[i,j] +
        sign*m*((kronecker_delta(i,j)*pracma::dot(offset_vec,offset_vec)) -
                  (offset_vec[i]*offset_vec[j]))
    }
  }
  return(I_new)
}
