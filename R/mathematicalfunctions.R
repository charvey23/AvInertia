# Script containing all basic mathematical functions
# all written by Christina Harvey
# last updated: 2020-08-03

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
  unit_vector = vector/pracma::Norm(vector)
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
  y_vector = pracma::cross(unit_z_vector,unit_x_vector)
  # ensure vectors are in unit form
  y_vector = calc_univec(y_vector)

  VRP2object = rbind(x_vector,y_vector,z_vector)

  return()
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
      I_off[i,j] = I_CG[i,j] + m*((dirac_delta(i,j)*pracma::dot(offset_vec,offset_vec)) - (offset_vec[i]*offset_vec[j]))
    }
  }

  return(I_off)
}



# ----------------------------------------------------------------------------------------
# ------------------------- Moment of interia tensors Properties  ------------------------
# ----------------------------------------------------------------------------------------

# ---------------------- Mass Properties - Solid cylinder -------------------------------
#' Moment of inertia tensor of a solid cylinder
#'
#' @param r radius of the cylinder
#' @param h height of the cylinder
#' @param m mass of the cylinder
#'
#' @author Christina Harvey
#'
#' @return Function returns the moment of inertia tensor of a solid cylinder about
#' its center of gravity with z oriented through it's major axis
#'
#'
#' @examples
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

# ---------------------- Mass Properties - Hollow cylinder -------------------------------
#' Moment of inertia tensor of a solid cylinder
#'
#' @param r_out outer radius of the cylinder
#' @param r_in inner radius of the cylinder
#' @param h height of the cylinder
#' @param m mass of the cylinder
#'
#' @author Christina Harvey
#'
#' @return Function returns the moment of inertia tensor of a hollow cylinder about
#' its center of gravity with z oriented through it's major axis
#'
#'
#' @examples
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

# -------------------- Mass Properties - solid square pyramid  -------------------------------
#' Moment of inertia tensor of a solid square pyramid
#'
#' @param w full width of one side of the pyramid base
#' @param h height of the pyramid
#' @param m mass of the pyramid
#'
#' @author Christina Harvey
#'
#' @return Function returns the moment of inertia tensor of a solid square pyramid cylinder about
#' its center of gravity with z oriented through it's major axis
#'
#'
#' @examples
#'
#' @export
#'
calc_inertia_pyrasolid <- function(w, h, m){

  temp_r = r_out^2 + r_in^2

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (w^2+0.5*h^2) # Ixx
  I[2,2] = (w^2+0.5*h^2) # Iyy
  I[3,3] = 2*(w^2)       # Izz
  I = (1/5)*m*I

  return(I)
}

# -------------------- Mass Properties - flat rectangular plate  -------------------------------
#' Moment of inertia tensor of a solid square pyramid
#'
#' @param w full width of one side of the pyramid base
#' @param h height of the pyramid
#' @param m mass of the pyramid
#'
#' @author Christina Harvey
#'
#' @return Function returns the moment of inertia tensor of a flat plate about
#' its center of gravity with z oriented parallel with it's long edge
#'
#'
#' @examples
#'
#' @export
#'
calc_inertia_platerect <- function(w, h, m){

  temp_r = r_out^2 + r_in^2

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (w^2+h^2) # Ixx
  I[2,2] = (h^2)     # Iyy
  I[3,3] = (w^2)     # Izz
  I = (1/12)*m*I

  return(I)
}


# -------------------- Mass Properties - flat triangular plate  -------------------------------
#' Moment of inertia tensor of a solid square pyramid
#'
#' @param pts a matrix of the three 3D points that define a point on the triangular plate;
#' each point should be a different row as follows:
#' pt1x, pt1y, pt1z
#' pt2x, pt1y, pt2z
#' pt3x, pt3y, pt3z
#'
#' @param a  area of the triangular plate
#' @param rho density of the material
#' @param t  thickness of the plate
#' @param desired_prop either "I" or "CG" depending on the desired output
#'
#' @author Christina Harvey
#'
#' @return Function returns the moment of inertia tensor of a flat triangular plate about
#' its center of gravity z axis defined by the input pts
#'
#' @section Warning:
#' The input points should be defined in a counterclockwise direction around the plate
#' in the triangular plate frame of reference
#' i.e. all z components should be equal
#'
#'
#' @examples
#'
#' @export
#'
calc_inertia_platetri <- function(pts, a, rho, t, desired_prop){

  pts = rbind(pts,pts[1,]) # need to add the first point to allow circular calculation

  # calculate the center of gravity of the plate
  x1  = 0
  y1  = 0
  for (i in 1:3){
    x1   = x1 + (1/(6*a))*(((pts[i,1]*pts[i+1,2])-(pts[i+1,1]*pts[i,2]))*(pts[i,1]+pts[i+1,1]));
    y1   = y1 + (1/(6*a))*(((pts[i,1]*pts[i+1,2])-(pts[i+1,1]*pts[i,2]))*(pts[i,2]+pts[i+1,2]));
  }

  CG = c(x1, y1, pts[1,3]);

  # adjust the incoming points to be so that the origin of the pts is on the plate's center of gravity
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
    x2   = x2  + (temp*(adj_pts[i,1]^2+(adj_pts[i,1]*adj_pts[i+1,1])+adj_pts[i+1,1]^2));
    y2   = y2  + (temp*(adj_pts[i,2]^2+(adj_pts[i,2]*adj_pts[i+1,2])+adj_pts[i+1,2]^2));
    Ixy  = Ixy + (0.5)*(temp*((2*adj_pts[i,1]*adj_pts[i,2])+ adj_pts[i,1]*adj_pts[i+1,2] + adj_pts[i,2]*adj_pts[i+1,1] + (2*adj_pts[i+1,1]*adj_pts[i+1,2])));
  }


  I = matrix(0, nrow = 3, ncol = 3) # define matrix
  I[1,1] = y2      # Ixx
  I[2,2] = x2      # Iyy
  I[3,3] = (x2+y2) # Izz
  I[1,2] = Ixy     # Ixy
  I[2,1] = Ixy     # Iyx

  I = (1/12)*rho*t*I


  if(desired_prop == "I"){
    return(I)
  }
  if(desired_prop == "CG"){
    return(CG)
  }
}

