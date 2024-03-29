% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anatomicalmassprop.R
\name{structural2VRP_feat}
\alias{structural2VRP_feat}
\title{Transform feather specific inertial properties to current position}
\usage{
structural2VRP_feat(m_f, I_fCG, CG_start, start, end, normal)
}
\arguments{
\item{m_f}{a scalar representing the mass of the feather (kg)}

\item{I_fCG}{a 3x3 matrix representing the moment of inertia tensor with the
origin at the feather calamus end and within the feather frame of reference}

\item{CG_start}{a 1x3 vector representing the center of gravity of the
feather with the origin at the feather calamus end and within the feather
frame of reference}

\item{start}{a 1x3 vector representing the location of the feather calamus
end with the origin at the VRP and within the full bird frame of reference}

\item{end}{a 1x3 vector representing the location of the feather
tip with the origin at the VRP and within the full bird frame of reference}

\item{normal}{a 1x3 vector representing the normal to the plane of the
feather vanes within the full bird frame of reference}
}
\value{
a list that includes:
\itemize{
\item{I}{a 3x3 matrix representing the moment of inertia tensor of a
simplified feather with the origin at the VRP and within the full bird
frame of reference}
\item{CG}{a 1x3 vector representing the center of gravity position of a
simplified feather with the origin at the VRP and within the full bird
frame of reference}
\item{m}{a double that returns the feather mass}
}
}
\description{
Transform feather specific inertial properties to current position
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
