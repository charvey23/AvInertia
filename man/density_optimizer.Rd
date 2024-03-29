% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anatomicalmassprop.R
\name{density_optimizer}
\alias{density_optimizer}
\title{Optimize torso section densities}
\usage{
density_optimizer(x, m_body, A, v_ell, CG_body_x, CG_ell, rho_avg)
}
\arguments{
\item{x}{a scalar guess at the density of the torso hemi-ellipsoid section (kg/m^3)}

\item{m_body}{a scalar representing the mass of the full torso}

\item{A}{a 2x2 matrix representing the mass and volume calculations of the
final two torso section}

\item{v_ell}{a scalar representing the volume of the hemi-ellipsoid}

\item{CG_body_x}{a scalar representing the position of the center of gravity
of the full torso along the x axis with the origin at the VRP}

\item{CG_ell}{a scalar representing the position of the center of gravity
of the hemi-ellipsoid section along the x axis with the origin at the VRP}

\item{rho_avg}{average density of the full torso}
}
\value{
the squared error between the three section densities and the
average torso density
}
\description{
Function that is used within an optimization routine to select the
appropriate torso section density
}
\examples{
# refer to the vignette
library(AvInertia)

# load data
data(dat_id_curr, package = "AvInertia")
data(dat_bird_curr, package = "AvInertia")
data(dat_feat_curr, package = "AvInertia")
data(dat_bone_curr, package = "AvInertia")
data(dat_mat, package = "AvInertia")
data(clean_pts, package = "AvInertia")

# 1. Determine the center of gravity of the bird's torso (including the legs)
dat_torsotail_out = massprop_restbody(dat_id_curr, dat_bird_curr)
# 2. Calculate the inertia of the flight feathers about the tip of the calamus
feather_inertia <- compute_feat_inertia(dat_mat, dat_feat_curr, dat_bird_curr)
# 3. Determine the center of gravity of one of the bird's wings
dat_wing_out      = massprop_birdwing(dat_id_curr, dat_bird_curr,
dat_bone_curr, dat_feat_curr, dat_mat, clean_pts,
feather_inertia, plot_var = 0)
# Visualize the center of gravity of each wing component in the x and y axis
dat_wing_out      = massprop_birdwing(dat_id_curr, dat_bird_curr,
dat_bone_curr, dat_feat_curr, dat_mat, clean_pts,
feather_inertia, plot_var = "yx")
# or the y and z axis
dat_wing_out      = massprop_birdwing(dat_id_curr, dat_bird_curr,
dat_bone_curr, dat_feat_curr, dat_mat, clean_pts,
feather_inertia, plot_var = "yz")
# 4. Combine all data and obtain the center of gravity, moment of inertia
# and principal axes of the bird
curr_full_bird      = combine_inertialprop(dat_torsotail_out,dat_wing_out,
dat_wing_out, dat_id_curr, dat_bird_curr, symmetric=TRUE)

}
