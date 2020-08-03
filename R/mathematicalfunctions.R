# Script containing all basic mathematical functions
# all written by Christina Harvey
# last updated: 2020-08-03

library(pracma) # where should this go

# ---------------------- Dirac delta function ---------------------------
#' Dirac delta function
#'
#' @param i  an index (usually row of a matrix)
#' @param j  an index (usually column of a matrix)
#'
#' @author Christina Harvey
#'
#' @return If i and j are equal dirac returns 1 otherwise returns 0
#'
#' @examples
#'
#'  # should return 1
#'  dirac_delta(7,7)
#'
#'  # should return 0
#'  dirac_delta(5,4)
#'
#' @export
#'
dirac_delta <- function(i,j){

  if (i == j) {
    return(1)
  }
  else {
    return(0)
  }
}

# ---------------------- Calculate a unit vector --------------------------
#' Determine the unit vector of any input vector
#'
#' @param vector
#'
#' @author Christina Harvey
#'
#' @return This function returns the unit vector of any input vector
#'
#' @examples
#'
#' library(pracma)
#'
#' output = unit_vector(c(2,6,5))
#' #should equal one if a unit vector
#' length_output = Norm(output)
#'
#' @export

calc_univec <- function(vector){
  unit_vector = vector/Norm(vector)
  return(unit_vector)
}


# ---------------------- Calculate a rotation matrices ---------------------------
#' Title
#'
#' @param z_vector
#' @param x_vector
#'
#' @author Christina Harvey
#'
#' @return Function returns the cosine rotation matrices to transform between VRP frame and object frame
#'
#' @examples
#'
#' @export
#'
calc_rot <- function(z_vector, x_vector){

  #CAUTION: incoming vectors must be in the structural frame at the VRP
  #         x and z vectors must already be orthogonal axes

  # ensure vectors are in unit form
  unit_z_vector = calc_univec(z_vector)
  unit_x_vector = calc_univec(x_vector)

  #cross z with x to get the righthanded axis
  y_vector = cross(unit_z_vector,unit_x_vector)
  # ensure vectors are in unit form
  y_vector = calc_univec(y_vector)

  rotation = list()
  rotation$VRP2object = rbind(x_vector,y_vector,z_vector)
  rotation$object2VRP = t(rotation$VRP2object)

  return(rotation)
}

# ---------------------- Calculate a matrix translation with parallel axis theorem  ---------------------------
#' Reads in an initial tensor and an offset to compute the transformed tensor.
#' Will be in the same frame of reference as the input tensor.
#'
#' @param I_CG
#' @param offset_vec
#' @param m
#'
#' @author Christina Harvey
#'
#' @return Function returns the transformed tensor after a solid body translation
#'
#'
#' @examples
#'
#' @export
#'

parallelaxis <- function(I_CG, offset_vec, m){

  # CAUTION: the parallel axis theorem only works if the I_CG is given about
  #          the object's centroidal axis

  I_off = matrix(0, nrow = 3, ncol = 3) # predefine matrix

  for (i in 1:3){
    for (j in 1:3){
      I_off[i,j] = I_CG[i,j] + m*((dirac_delta(i,j)*dot(offset_vec,offset_vec)) - (offset_vec[i]*offset_vec[j]))
    }
  }

  return(I_off)
}


# ---------------------- Mass Properties - Solid cylinder -------------------------------
