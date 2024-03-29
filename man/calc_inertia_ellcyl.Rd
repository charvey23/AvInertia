% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inertiafunctions.R
\name{calc_inertia_ellcyl}
\alias{calc_inertia_ellcyl}
\title{Moment of inertia tensor of a solid elliptical cylinder
Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf}
\usage{
calc_inertia_ellcyl(a, b, l, m)
}
\arguments{
\item{a}{half height of the base - oriented along the x axis in torso FOR
(z axis in full bird FOR) (m)}

\item{b}{half width of the base - oriented along the y axis in torso FOR
(y axis in full bird FOR)  (m)}

\item{l}{length of the cylinder - oriented along the z axis in torso FOR
(x axis in full bird FOR)  (m)}

\item{m}{mass of the cylinder (kg)}
}
\value{
a 3x3 matrix representing the moment of inertia tensor of a solid
elliptical cylinder about it's center of gravity
}
\description{
Moment of inertia tensor of a solid elliptical cylinder
Reference: https://apps.dtic.mil/sti/pdfs/AD0274936.pdf
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
\author{
Christina Harvey
}
