# Script containing all basic mathematical functions
# all written by Christina Harvey
# last updated: 2020-08-03

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
#' @param z_vector a 1x3 vector representing the direction for the desired z axis of the new frame of reference.
#' @param x_vector a 1x3 vector representing the direction for the desired x axis of the new frame of reference.
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

  #cross z with x to get the righthanded axis
  y_vector = pracma::cross(unit_z_vector,unit_x_vector)
  # ensure vectors are in unit form
  unit_y_vector = calc_univec(y_vector)

  # rotation matrix representing the rotated basis
  VRP2object = rbind(unit_x_vector,unit_y_vector,unit_z_vector)

  return(VRP2object)
}

# ---------------------- Calculate a matrix translation with parallel axis theorem  ---------------------------
#' Parallel axis theory
#'
#' Reads in an initial tensor and an offset to compute the transformed tensor.
#' Will be in the same frame of reference as the input tensor.
#'
#' @param I_CG Moment of intertia tensor (3x3) about the center of gravity of the object (kg-m^2)
#' @param offset_vec Distance between the objects CG and the arbitrary pt A. (m) Should always point from the CG to the arbitrary point A
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



# ----------------------------------------------------------------------------------------
# ------------------------- Moment of inertia tensors Properties  ------------------------
# ----------------------------------------------------------------------------------------

# ---------------------- Mass Properties - Solid cylinder -------------------------------
#' Moment of inertia tensor of a solid cylinder
#'
#' @param r radius of the cylinder (m)
#' @param h height of the cylinder (m)
#' @param m mass of the cylinder (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid cylinder about
#' its center of gravity with z oriented through it's major axis
#'
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
#' Moment of inertia tensor of a hollow cylinder
#'
#' @param r_out outer radius of the cylinder (m)
#' @param r_in inner radius of the cylinder (m)
#' @param h height of the cylinder (m)
#' @param m mass of the cylinder (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a hollow cylinder about
#' its center of gravity with z oriented through it's major axis
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
#' All outputs are based on an origin at the centered point on the base
#'
#' @param w_2 half width of one side of the pyramid base (m)
#' @param h height of the pyramid (m)
#' @param m mass of the pyramid (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid square pyramid cylinder about
#' its center of gravity with z oriented through it's major axis. Origin is NOT at the center of gravity but at the center of the base.
#'
#' @export
#'
calc_inertia_pyrasolid <- function(w_2, h, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = (w_2^2+0.5*h^2) # Ixx
  I[2,2] = (w_2^2+0.5*h^2) # Iyy
  I[3,3] = 2*(w_2^2)       # Izz
  I = (1/5)*m*I

  return(I)
}

# -------------------- Mass Properties - flat rectangular plate  -------------------------------
#' Moment of inertia tensor of a flat rectangular plate
#'
#' @param w full width of one side of the plate (m)
#' @param h height of the plate (m)
#' @param m mass of the plate (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a flat plate about
#' its center of gravity with z oriented parallel with it's long edge
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


# -------------------- Mass Properties - flat triangular plate  -------------------------------
#' Moment of inertia tensor of a flat triangular plate
#'
#' @param pts a matrix of the three 3D points that define a point on the triangular plate.
#' Frame of reference: Muscle | Origin: VRP
#' each point should be a different row as follows:
#' pt1x, pt1y, pt1z
#' pt2x, pt1y, pt2z
#' pt3x, pt3y, pt3z
#'
#' @param a  area of the triangular plate (m)
#' @param rho density of the material (kg/m^3)
#' @param t  thickness of the plate (m)
#' @param desired_prop a string containing either "I" or "CG" depending on the desired output
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a flat triangular plate about
#' its center of gravity z axis defined by the input pts
#'
#' @section Warning:
#' The input points should be defined in a counterclockwise direction around the plate
#' in the triangular plate frame of reference
#' i.e. all z components should be equal
#'
#' @export
#'
calc_inertia_platetri <- function(pts, A, rho, t, desired_prop){

  pts = rbind(pts,pts[1,]) # need to add the first point to allow circular calculation

  # calculate the center of gravity of the plate
  x1  = 0
  y1  = 0
  for (i in 1:3){
    x1   = x1 + (1/(6*A))*(((pts[i,1]*pts[i+1,2])-(pts[i+1,1]*pts[i,2]))*(pts[i,1]+pts[i+1,1]));
    y1   = y1 + (1/(6*A))*(((pts[i,1]*pts[i+1,2])-(pts[i+1,1]*pts[i,2]))*(pts[i,2]+pts[i+1,2]));
  }

  CG = c(x1, y1, pts[1,3]); # Frame of reference: Incoming points | Origin: Incoming points

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
#' @param Pt11 a matrix of 3 values representing the point on tip of most proximal feather (wing root trailing edge) (m)
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
  # Note: Secondaries will at most be in line with the elbow point along the
  #       vector drawn between S1 and the wing root trailing edge
  sec_vec = Pt11-Pt10;                   # vector between S1 and the wing root trailing edge
  t       = (Pt2[2]-Pt10[2])/sec_vec[2]; # proportion along the vector where it intersects y = Pt2y
  sec_end = c((t*sec_vec[1]+Pt10[1]), Pt2[2], (t*sec_vec[3]+Pt10[3])); # Point along the vector at Pt2Y

  # Calculate the start and end of the secondaries
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
