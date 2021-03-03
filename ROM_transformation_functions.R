


## This code reads in three points in 3D and computs the angle where the second point is the central point
## Requires standard deviation on the calculation as well.

## To calculate the reflex angle this must be checked for each specimen independently due to the orientation and truncation errors

## Written by: Christina Harvey
## Date: Feb 2020

#' Title
#'
#' @param x1 1st component of the 1st point
#' @param y1 2nd component of the 1st point
#' @param z1 3rd component of the 1st point
#' @param x2 1st component of the 2nd point
#' @param y2 2nd component of the 2nd point
#' @param z2 3rd component of the 2nd point
#' @param x3 1st component of the 3rd point
#' @param y3 2nd component of the 3rd point
#' @param z3 3rd component of the 3rd point
#'
#' @return the interior angle between the input three points
#' @export
#'
#' @examples
jointangles <- function(x1,y1,z1,x2,y2,z2,x3,y3,z3){

  #pt 2 is the central point to compute the angle around
  # Calculate the vectors
  u1 = c(x1,y1,z1)-c(x2,y2,z2)
  v1 = c(x3,y3,z3)-c(x2,y2,z2)

  interior = dot(u1,v1)/(Norm(u1)*Norm(v1))

  angle = acos(interior)*180/pi

  return(angle)
}


allDup <- function (value) # found online by jholtman at gmail.com
{duplicated(value) | duplicated(value, fromLast = TRUE)}




################################ Custom functions ##############################
## Just source these and move on to the next section

## This function scales wings up or down based on a scale factor
scale_the_wing <- function(specimens, ## MUST BE $rotated component of procSym output
                           scale_factor ## numeric scaling factor
) {

  ## figure out dimensionality of array
  n <- dim(specimens)[3]
  k <- dim(specimens)[1]
  m <- dim(specimens)[2]
  x <- specimens

  resizelist <-  array(numeric(),c(k, m ,n))
  for (i in 1:n) {
    resizelist[,,i] <-
      x[,,i] * scale_factor
  }

  ## export as a 3D array
  return(resizelist)
}

## This function re-sizes the optitrack data for a specimen (specimen_to_adjust)
## based on the size of a second specimen (target_specimen)
resize_to <- function(specimen_to_adjust,
                      target_specimen,
                      tol = 5e-3){

  ## make names shorter for ease of coding later
  x1 <- specimen_to_adjust %>%
    select(-c(frame, time_sec)) %>%
    drop_na() %>%
    arrayspecs(p = 11, k = 3) %>%
    procSym(scale = FALSE)
  x2 <- target_specimen %>%
    select(-c(frame, time_sec)) %>%
    drop_na() %>%
    arrayspecs(p = 11, k = 3) %>%
    procSym(scale = FALSE)

  ## get into 2d array
  x1_two_d <- two.d.array(x1$rotated)
  x2_two_d <- two.d.array(x2$rotated)

  ## get the length of the humerus in each specimen
  x1_humerus_lengths <-
    calc_dist(x1_two_d[,7:12])
  x2_humerus_lengths <-
    calc_dist(x2_two_d[,7:12])

  ## set the scale_factor based on the means of the humerus lengths
  scale_factor <- mean(x2_humerus_lengths)/mean(x1_humerus_lengths)

  ## rescale the wing according to the scale factor
  rescaled <- scale_the_wing(x1$rotated, scale_factor = scale_factor)

  rescaled_two_d <- two.d.array(rescaled)

  ## check humerus lengths
  rescaled_hums <-
    calc_dist(rescaled_two_d[,7:12])

  hum_diffs <- mean(rescaled_hums) - mean(x2_humerus_lengths)


  if (abs(hum_diffs) > tol) {
    stop("shit happened. give up.")
  }

  ## isolate the frames and times
  spec_framtimes <-
    specimen_to_adjust %>%
    drop_na() %>%
    select(frame, time_sec)

  ## paste it back together
  export_dat <-
    data.frame(spec_framtimes,
               rescaled_two_d)
  colnames(export_dat) <- colnames(specimen_to_adjust)

  ## export
  return(export_dat)

}

################################## calc_dist ###################################
# this function calculates the 3D vector between two points
# input: a dataframe with 6 columns, first 3 columns are coordinates for one point
# last 3 columns are coordinates for another point
# output: a numeric value
calc_dist <- function(pair) {

  dx = pair[,1]-pair[,4]
  dy = pair[,2]-pair[,5]
  dz = pair[,3]-pair[,6]
  dist = sqrt(dx^2+dy^2+dz^2)

  return(dist)
}
