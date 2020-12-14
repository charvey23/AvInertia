

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


# -------------------- Mass Properties - solid ellipse   -------------------------------
#' Moment of inertia tensor of solid ellipse CG or a half ellipse centered on the base
#'
#' All outputs are based on an origin at the centered point on the base
#'
#' @param a half the height along the x direction (m)
#' @param b half the width along the y direction  (m)
#' @param c half the length along the z direction (m)
#' @param m mass of the ellipse (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid ellipse about
#' its center of gravity with the major axes aligned. Origin is at the center of gravity for
#' a full ellipse or at the center of the base if modelling a half ellipse.
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

# -------------------- Mass Properties - solid circular cone   -------------------------------
#' Moment of inertia tensor of a solid circular cone pyramid
#'
#' All outputs are based on an origin at the centered point on the base
#'
#' @param r radius of the cone base (m)
#' @param h height of the cone (m)
#' @param m mass of the cone (kg)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid circular cone about
#' its center of gravity with z oriented through it's major axis. Origin is NOT at the center of gravity but at the center of the base.
#'
#' @export
#'
calc_inertia_conesolid <- function(r, h, m){

  I = matrix(0, nrow = 3, ncol = 3)
  I[1,1] = ((1.5*r^2)+h^2) # Ixx
  I[2,2] = ((1.5*r^2)+h^2) # Iyy
  I[3,3] = 3*(r^2)         # Izz
  I = (1/10)*m*I

  return(I)
}

# -------------------- Mass Properties - partial solid elliptical cone  -------------------------------
#' Moment of inertia tensor of a partial solid elliptical cone
#'
#' @param w   full width of the wide base (m)
#' @param h   full width of the wide base (m)
#' @param lp  length of the partial cone (m)
#' @param lf  length of the full cone as if it had not been cut (m)
#' @param rho density of the partial cone (kg/m^3)
#'
#' @author Christina Harvey
#'
#' @return a 3x3 matrix representing the moment of inertia tensor of a partial solid elliptical cone about
#' the center of the wider base with z oriented towards the end
#'
#' @export
#'
calc_inertia_ellcone <- function(w, h, lp, lf, rho){

  I = matrix(0, nrow = 3, ncol = 3)
  tmpxy = (((lp-lf)^5/lf^4)+lf)
  tmpz  = (lp^3/(3*lf^2))*((6*lp^2)-(15*lp*lf)+(10*lf^2))
  I[1,1] = rho*w*h*pi*(tmpxy*(w^2/8) + tmpz) # Ixx
  I[2,2] = rho*w*h*pi*(tmpxy*(h^2/8) + tmpz) # Iyy
  I[3,3] = (1/8)*rho*w*h*pi*tmpxy*(w^2+h^2)  # Izz
  I = (1/40)*I

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
#' @return a 3x3 matrix representing the moment of inertia tensor of a solid square pyramid about
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
#' its center of gravity with z oriented parallel with it's long edge and y along its short edge
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

