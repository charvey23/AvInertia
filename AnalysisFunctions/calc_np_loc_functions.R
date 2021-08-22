## ----------------------------------------------------------------
## --------- Define functions to determine key values -------------
## ----------------------------------------------------------------

find_c4_x <- function(peri_pts,target_chord){
  distal_y_pt = which.max(peri_pts[,2])
  no_dis = 1000
  y_test = seq(from=min(peri_pts[,2]),to=max(peri_pts[,2]),length.out=no_dis)
  step_size = ((max(peri_pts[,2]))/(no_dis - 1))
  LE_x = c(1:no_dis)*0
  TE_x = c(1:no_dis)*0
  LE_z = c(1:no_dis)*0
  TE_z = c(1:no_dis)*0
  chord = c(1:no_dis)*0
  c2 = 0
  c1 = 0
  cy = 0
  cx = 0
  for(j in 1:no_dis){
    LE_x[j] = calc_LE_xz(y_test[j],distal_y_pt,peri_pts,"x")
    TE_x[j] = calc_TE_xz(y_test[j],distal_y_pt,peri_pts,"x")
    LE_z[j] = calc_LE_xz(y_test[j],distal_y_pt,peri_pts,"z")
    TE_z[j] = calc_TE_xz(y_test[j],distal_y_pt,peri_pts,"z")
    chord[j] = sqrt((LE_x[j]-TE_x[j])^2+(LE_z[j]-TE_z[j])^2)

    if(j < no_dis){
      c2 = c2 + chord[j]^2*step_size
      c1 = c1 + chord[j]*step_size
      cy = cy + chord[j]*y_test[j]*(1+0.5*step_size)*step_size
      theta = atan((TE_z[j]-LE_z[j])/(TE_x[j]-LE_x[j]))
      if(TE_x[j] > LE_x[j]){ # accounts for situations where the TE is in front of the LE and so the c/4 is in front of the LE
        cx = cx + chord[j]*(LE_x[j] + 0.25*chord[j]*cos(theta))*step_size
      } else{
        cx = cx + chord[j]*(LE_x[j] - 0.25*chord[j]*cos(theta))*step_size
      }

    }
  }

  # identify span location where chord = mean chord

  if (target_chord == "MAC"){
    target_chord = c2/c1
  }

  if (target_chord == "area centroid"){
    target_y_chord = cy/c1
    LE_x_cen = calc_LE_xz(target_y_chord,distal_y_pt,peri_pts,"x")
    TE_x_cen = calc_TE_xz(target_y_chord,distal_y_pt,peri_pts,"x")
    LE_z_cen = calc_LE_xz(target_y_chord,distal_y_pt,peri_pts,"z")
    TE_z_cen = calc_TE_xz(target_y_chord,distal_y_pt,peri_pts,"z")
    target_chord = sqrt((LE_x_cen-TE_x_cen)^2+(LE_z_cen-TE_z_cen)^2)

    chord_LE_x = LE_x_cen
    # account for the twist of the chord
    theta = atan((TE_z_cen-LE_z_cen)/(TE_x_cen-LE_x_cen))
  } else if(target_chord != "standard mean"){
    # if not centroid need to find the approximate y location
    neg_val = which(chord-target_chord < 0)
    y_mean = which(neg_val != row(as.matrix(neg_val)))[1] # pulls out the first time that the result crosses zero
    if (y_mean == 1){
      pos_val = which(chord-target_chord > 0)
      y_mean = which(pos_val != row(as.matrix(pos_val)))[1] # pulls out the first time that the result crosses zero
      if(is.na(y_mean)){
        y_mean = max(pos_val) + 1
      }
    }
    chord_LE_x = LE_x[y_mean]
    # account for the twist of the chord
    theta = atan((TE_z[y_mean]-LE_z[y_mean])/(TE_x[y_mean]-LE_x[y_mean]))
  }

  if (target_chord == "standard mean"){
    c4 = cx/c1
  } else{
    c4 = chord_LE_x - 0.25*target_chord*cos(theta)
  }

  return(c4)
}

calc_LE_xz <- function(LE_y,distal_pt,peri_pts,x_or_z){
  if(x_or_z =="x"){
    coord = 1
  }else{
    coord = 3
  }
  # Define linear portion interior to Pt6
  if (LE_y <= peri_pts[2,2] & distal_pt >= 1){
    LE = calc_lin_edge(peri_pts[1,coord],peri_pts[2,coord],
                       peri_pts[1,2],peri_pts[2,2],LE_y)}
  # Define linear portion interior to Pt7
  if (LE_y <= peri_pts[3,2] & LE_y > peri_pts[2,2] & distal_pt >= 2){
    LE = calc_lin_edge(peri_pts[2,coord],peri_pts[3,coord],
                       peri_pts[2,2],peri_pts[3,2],LE_y)}
  # Define linear portion interior to Pt8
  if (LE_y <= peri_pts[4,2] & LE_y > peri_pts[3,2] & distal_pt >= 3){
    LE = calc_lin_edge(peri_pts[3,coord],peri_pts[4,coord],
                       peri_pts[3,2],peri_pts[4,2],LE_y)}
  # Define linear portion interior to Pt9
  if (LE_y <= peri_pts[5,2] & LE_y > peri_pts[4,2] & distal_pt >= 4){
    LE = calc_lin_edge(peri_pts[4,coord],peri_pts[5,coord],
                       peri_pts[4,2],peri_pts[5,2],LE_y)}
  # Define linear portion interior to Pt10
  if (LE_y <= peri_pts[6,2] & LE_y > peri_pts[5,2] & distal_pt >= 5){
    LE = calc_lin_edge(peri_pts[5,coord],peri_pts[6,coord],
                       peri_pts[5,2],peri_pts[6,2],LE_y)}
  return(LE)
}

calc_TE_xz <- function(TE_y,distal_pt,peri_pts,x_or_z){
  if(x_or_z =="x"){
    coord = 1
  }else{
    coord = 3
  }
  # Define linear portion interior to Pt11
  if(nrow(peri_pts) == 8){
    if (TE_y <= peri_pts[7,2] & distal_pt <= 7){
      TE = calc_lin_edge(peri_pts[8,coord],peri_pts[7,coord],
                         peri_pts[8,2],peri_pts[7,2],TE_y)}
    # Define linear portion interior to Pt10
    if (TE_y <= peri_pts[6,2] & TE_y > peri_pts[7,2] & distal_pt <= 6){
      TE = calc_lin_edge(peri_pts[7,coord],peri_pts[6,coord],
                         peri_pts[7,2],peri_pts[6,2],TE_y)}
  } else{
    # Define linear portion interior to Pt10
    if (TE_y <= peri_pts[6,2] & distal_pt <= 6){
      TE = calc_lin_edge(peri_pts[7,coord],peri_pts[6,coord],
                         peri_pts[7,2],peri_pts[6,2],TE_y)}
  }

  # Define linear portion interior to Pt9
  if (TE_y <= peri_pts[5,2] & TE_y > peri_pts[6,2] & distal_pt <= 5){
    TE = calc_lin_edge(peri_pts[6,coord],peri_pts[5,coord],
                       peri_pts[6,2],peri_pts[5,2],TE_y)}
  # Define linear portion interior to Pt8
  if (TE_y <= peri_pts[4,2] & TE_y > peri_pts[5,2] & distal_pt <= 4){
    TE = calc_lin_edge(peri_pts[5,coord],peri_pts[4,coord],
                       peri_pts[5,2],peri_pts[4,2],TE_y)}
  # Define linear portion interior to Pt7
  if (TE_y <= peri_pts[3,2] & TE_y > peri_pts[4,2] & distal_pt <= 3){
    TE = calc_lin_edge(peri_pts[4,coord],peri_pts[3,coord],
                       peri_pts[4,2],peri_pts[3,2],TE_y)}
  # Define linear portion interior to Pt6
  if (TE_y <= peri_pts[2,2] & TE_y > peri_pts[3,2] & distal_pt <= 2){
    TE = calc_lin_edge(peri_pts[3,coord],peri_pts[2,coord],
                       peri_pts[3,2],peri_pts[2,2],TE_y)}
  return(TE)
}


calc_lin_edge <- function(x1,x2,y1,y2,y_curr){
  m = (x1-x2)/(y1-y2)
  b = x1 - m*y1
  x_curr = m*y_curr+b
  return(x_curr)
}
