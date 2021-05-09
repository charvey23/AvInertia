


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
## based on the size of a second specimen (target_length)
resize_to <- function(specimen_to_adjust,char_colnames,adjust_length,target_length,tol = 5e-3){

  x1 <- specimen_to_adjust %>%
    select(-char_colnames) %>%
    drop_na() %>%
    arrayspecs(p = 11, k = 3) %>%
    procSym(scale = FALSE, CSinit=FALSE, reflect = FALSE)
  ## get into 2d array
  x1_two_d <- two.d.array(x1$rotated)

  ## set the scale_factor based on the means of the bone lengths
  # length of target bone/length of bone in wing that will be resized
  scale_factor <- target_length/adjust_length

  ## rescale the wing according to the scale factor
  rescaled <- scale_the_wing(x1$rotated, scale_factor = scale_factor)

  rescaled_two_d <- two.d.array(rescaled)

  ## isolate the frames and times
  spec_framtimes <-
    specimen_to_adjust %>%
    drop_na() %>%
    select(char_colnames)

  ## paste it back together
  export_dat           <- data.frame(spec_framtimes, rescaled_two_d)
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




subsample_frames <- function(dat_keep,dat_dup,curr_BirdID,max_sample_no){

  dat_dup_uni   <- unique(dat_dup[,c("elbow_round","manus_round")])
  # Loop through each duplicated values and save one random sample
  for (j in 1:nrow(dat_dup_uni)){
    # limit to the current unique bin for the specific individual we are looking at
    tmp = subset(dat_dup, elbow_round == dat_dup_uni$elbow_round[j] &
                   manus_round == dat_dup_uni$manus_round[j] &
                   BirdID == curr_BirdID)

    # if there are more frames than desired samples than choose at random
    if (nrow(tmp) > max_sample_no){
      tmp = subset(tmp, FrameID %in% sample(tmp$FrameID, max_sample_no))
    }
    # if there are not enough frames than desired samples than keep those but add in randomly sampled frames from other individuals
    if (nrow(tmp) < max_sample_no){
      tmp_oth = subset(dat_dup, elbow_round == dat_dup_uni$elbow_round[j] &
                         manus_round == dat_dup_uni$manus_round[j] &
                         BirdID != curr_BirdID)
      row_rand = sample(tmp_oth[,c("FrameID","BirdID")], max_sample_no-nrow(tmp))
      tmp_oth  = subset(tmp_oth, FrameID %in% sample(tmp_oth$FrameID, max_sample_no-nrow(tmp)))
      tmp      = rbind(tmp,tmp_oth)
    }

    # Save the selected rows
    if (j == 1){
      dat_subsample = rbind(dat_keep,tmp)
    }else{
      dat_subsample = rbind(dat_subsample,tmp)
    }

  } # end of the loop through joint angle bins

  return(dat_subsample)
}



# ---- This code is to be used with the LLT takes raw digitized points and prepares them for testing
#When re-running for new data double check that the transformations are working correctly.

# CAUTION: This code is now set up to orient wings as follows:
#          1. Pt3 (Wrist) will be inline with Pt1 (humeral head) along the y axis
#          2. Pt3 (Wrist) will be inline with Pt1 (humeral head) along the z axis
#          3. Pt3 (Wrist) will be inline with Pt10 (S1) along the x axis

# Created: Christina Harvey
# Last updated: 07-May-2021

reorient_wings<- function(dat_subsample){

  # ----------------- Define set variables ------------
  point_list <- c("pt1","pt2","pt3","pt4","pt6","pt7","pt8","pt9","pt10","pt11","pt12")
  dim_list   <- c("X","Y","Z")
  xaxis <- c(1, 0, 0)
  yaxis <- c(0, 1, 0)
  zaxis <- c(0, 0, 1)

  ## ------- Reorient Wings ---------
  dat_clean <- dat_subsample

  #------------------------------- Step 1 -------------------------------
  ##### error doesn't adjust all points ######
  #Make pt 1 on humerus the beginning location
  for (j in 2:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      pt1_name = paste(point_list[1],dim_list[k],sep = "_")
      dat_clean[,col_name] <- dat_clean[,col_name] - dat_clean[,pt1_name]
    }
  }
  dat_clean[,c("pt1_X","pt1_Y","pt1_Z")] <- 0
  #-------------------------------  Step 2 -------------------------------
  # Set pt 3 in line with Pt 1 along the wingspan
  #Calculate angle between axis that runs along wing length and pt 3 (rotate about z)
  #Project onto xaxis & make negative if j1 > 0:
  i1 = dat_clean$pt3_X
  j1 = dat_clean$pt3_Y

  dot_xaxis = i1
  vec_norm  = sqrt(i1^2 + j1^2)
  interior  = dot_xaxis/vec_norm
  thetaz    = acos(interior) + pi / 2
  thetaz[which(j1 > 0)] = -acos(interior[which(j1 > 0)]) + pi / 2

  # Rotate about the z axis
  tmp <- dat_clean
  for (j in 1:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      if (k == 1){
        dat_clean[,col_name] = cos(thetaz) * tmp[,col_name] - sin(thetaz) * tmp[,paste(point_list[j],dim_list[2],sep = "_")]
      } else{
        if(k == 2){
          dat_clean[,col_name] = sin(thetaz) * tmp[,paste(point_list[j],dim_list[1],sep = "_")] + cos(thetaz) * tmp[,col_name]
        }
      }
    }
  }
  # CAUTION: Verify that Pt3X=x_adj1=0 after this step - View(subset(pts,Pt.No == 3))

  #------------------------------- Step 3 -------------------------------
  # Set pt 3 in line with Pt 1 looking down the length of wing
  # Calculate angle between x axis and pt 3 (rotate about y)
  # Project onto xaxis & make negative if k2 > 0
  j2 = dat_clean$pt3_Y
  k2 = dat_clean$pt3_Z

  dot_yaxis = j2
  vec_norm  = sqrt(j2^2 + k2^2)
  interior = dot_yaxis/vec_norm

  thetax = acos(interior)
  thetax[which(k2 > 0)] = -acos(interior[which(k2 > 0)])

  # Rotate about the x axis
  tmp <- dat_clean
  for (j in 1:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      if (k == 2){
        dat_clean[,col_name] = cos(thetax) * tmp[,col_name] - sin(thetax) * tmp[,paste(point_list[j],dim_list[3],sep = "_")]
      } else{
        if(k == 3){
          dat_clean[,col_name] = sin(thetax) * tmp[,paste(point_list[j],dim_list[2],sep = "_")] + cos(thetax) * tmp[,col_name]
        }
      }
    }
  }
  # CAUTION: Verify that Pt3Z=z_adj2=0 after this step - View(subset(pts,Pt.No == 3))

  #------------------------------- Step 4 -------------------------------
  # Set pt 10 in line with Pt 3 looking from back of wing
  # Calculate angle between x axis and pt 10 (rotate about y)
  # Project onto yaxis & make negative if k3 > 0
  i3 = dat_clean$pt10_X
  k3 = dat_clean$pt10_Z

  dot_xaxis = i3
  vec_norm  = sqrt(i3^2 + k3^2)
  interior  = dot_xaxis/vec_norm

  thetay = -acos(interior)
  thetay[which(k3 > 0)] = acos(interior[which(k3 > 0)])

  # Rotate about the y axis
  tmp <- dat_clean
  for (j in 1:11){
    for (k in 1:3){
      col_name = paste(point_list[j],dim_list[k],sep = "_")
      if (k == 1){
        dat_clean[,col_name] = cos(thetay) * tmp[,col_name] + sin(thetay) * tmp[,paste(point_list[j],dim_list[3],sep = "_")]
      } else{
        if(k == 3){
          dat_clean[,col_name] = -sin(thetay) * tmp[,paste(point_list[j],dim_list[1],sep = "_")] + cos(thetay) * tmp[,col_name]
        }
      }
    }
  }
  # CAUTION: Verify that Pt10Z=z_adj3=0 after this step - View(subset(pts,Pt.No == 10))

  #------------------------------- Step 5 -------------------------------
  for (j in 1:11) {
    col_name = paste(point_list[j], dim_list[1], sep = "_")
    dat_clean[, col_name] <- -dat_clean[, col_name]

    col_name = paste(point_list[j], dim_list[3], sep = "_")
    dat_clean[which(dat_clean$wing == "R"), col_name] <- -dat_clean[which(dat_clean$wing == "R"), col_name]
  }

  return(dat_clean)

}
