#' The identified peripheral and joint 3D positions.
#'
#' A dataset containing the identified peripheral and joint 3D positions for a
#' pigeon BirdID = 20_0300, TestID = 20_03004, FrameID = 317
#'
#' @format A matrix with 10 rows and 3 columns. Columns represent x, y, and z coordinate respectively:
#' \describe{
#' \item{pt1x, pt1y, pt1z}{Point on the shoulder joint (m).}
#' \item{pt2x, pt1y, pt2z}{Point on the elbow joint (m).}
#' \item{pt3x, pt3y, pt3z}{Point on the wrist joint (m).}
#' \item{pt4x, pt4y, pt4z}{Point on the end of carpometacarpus (m).}
#' \item{pt6x, pt6y, pt6z}{Point on the leading edge of the wing in front of the
#' wrist joint (m).}
#' \item{pt8x, pt8y, pt8z}{Point on tip of most distal primary (m).}
#' \item{pt9x, pt9y, pt9z}{Point on the tip of the last primary to model as if
#' it is on the end of the carpometacarpus (m).}
#' \item{pt10x, pt10y, pt10z}{Point on tip of last primary to model as if
#' it was distributed along the carpometacarpus (m).}
#' \item{pt11x, pt11y, pt11z}{Point on tip of most proximal feather (m).}
#' \item{pt12x, pt12y, pt12z}{Point on exterior shoulder position
#' (wing root leading edge) (m).}
#' }
#'
"clean_pts"

#' Identification variables for current configuration
#'
#' A dataset containing the identification data for a
#' pigeon BirdID = 20_0300, TestID = 20_03004, FrameID = 317
#'
#' @format A dataframe with 1 row and 4 required variables:
#' \describe{
#' \item{species}{Species ID code as a string}
#' \item{BirdID}{Bird ID code as a string}
#' \item{TestID}{Test ID code as a string}
#' \item{frameID}{Video frame ID code as a string}
#' }
#'
"dat_id_curr"

#' Bird specific morphology dataset
#'
#' A dataset containing the bird specific morphology data for a
#' pigeon BirdID = 20_0300, TestID = 20_03004, FrameID = 317
#'
#' @format A dataframe with 1 row and 125 variables.
#' Only 29 required for proper operation:
#' \describe{
#' \item{total_bird_mass}{Mass of full bird for the current wing (kg)}
#' \item{wing_mass}{Mass of one wing, should be the current wing (kg)}
#' \item{barb_radius}{Radius of feather barb  for current species (m)}
#' \item{barb_distance}{Distance between feather barbs for current species (m)}
#' \item{brachial_muscle_mass}{Mass of all muscles in the brachial region
#' of the wing (kg)}
#' \item{antebrachial_muscle_mass}{Mass of all muscles in the antebrachial region
#' of the wing (kg)}
#' \item{manus_muscle_mass}{Mass of all muscles in the manus region
#' of the wing (kg)}
#' \item{all_skin_coverts_mass}{Mass of all skin and covert feathers (kg).}
#' \item{tertiary_mass}{Mass of tertiary feathers (kg).}
#' \item{extend_neck}{Logical input defining whether the neck should be
#' modeled as extended or not}
#' \item{head_length}{Length of the head from base to tip (m)}
#' \item{head_mass}{Mass of the head (kg)}
#' \item{head_height}{Height of the head at the base (m)}
#' \item{neck_mass}{Mass of the neck (kg)}
#' \item{neck_width}{OPTIONAL: Average width of the stretched neck (m)}
#' \item{neck_length}{OPTIONAL: Length of the stretched neck (m)}
#' \item{torsotail_length}{Length from the beginning of the torso to the tip of the tail (m)}
#' \item{torsotail_mass}{Mass of the torso and tail (kg)}
#' \item{tail_length}{Length of the tail (m)}
#' \item{tail_mass}{Mass of the tail (kg)}
#' \item{tail_width}{Average width of the furled tail (m)}
#' \item{right_leg_mass}{Mass of the right leg (kg)}
#' \item{left_leg_mass}{Mass of the left leg (kg)}
#' \item{body_width_max}{Maximum width of the torso (m)}
#' \item{body_width_at_leg_insert}{Width of the body at the point
#' where the legs are inserted (m).}
#' \item{x_loc_of_body_max}{x coordinate from the VRP in the
#' full bird frame of reference of the maximum body width (m).}
#' \item{x_loc_leg_insertion}{x coordinate from the VRP in the
#' full bird frame of reference of the leg insertion location (m).}
#' \item{x_loc_TorsotailCoG}{x coordinate from the VRP in the
#' full bird frame of reference of the center of gravity of the torso and tail (m).}
#' \item{z_loc_TorsotailCoG}{x coordinate from the VRP in the
#' full bird frame of reference of the the center of gravity of the torso and tail (m).}
#' }
#'
"dat_bird_curr"

#' Wing bone specific morphology dataset
#'
#' A dataset containing the wing bone specific data for a
#' pigeon BirdID = 20_0300, TestID = 20_03004, FrameID = 317
#'
#' @format A dataframe with 6 rows and 5 variables:
#' \describe{
#' \item{bone}{Bone ID code. Must include:
#' "Humerus","Ulna","Radius","Carpometacarpus","Ulnare" and "Radiale".}
#' \item{bone_mass}{Mass of bone in the same row as the appropriate
#' bone ID code (kg)}
#' \item{bone_len}{Length of bone in the same row as the appropriate
#' bone ID code (m)}
#' \item{bone_out_rad}{Outer radius of bone in the same row as the appropriate
#' bone ID code (m)}
#' \item{bone_in_rad}{Inner radius of bone in the same row as the appropriate
#' bone ID code (m)}
#' }
#'
"dat_bone_curr"

#' Flight feather specific morphology dataset
#'
#' A dataset containing the wing bone specific data for a
#' pigeon BirdID = 20_0300, TestID = 20_03004, FrameID = 317
#'
#' @format A dataframe with 20 rows and 15 variables.
#' Only 8 variables required for proper functioning:
#' \describe{
#' \item{feather}{Feather ID code. Must be in standard format i.e.
#' 1st primary is "P1", third secondary is "S3", etc.
#' Alula feathers should be grouped and named "alula".}
#' \item{m_f}{Mass of feather in the same row as the
#' appropriate feather ID code (kg)}
#' \item{l_cal}{Length of calamus in the same row as the
#' appropriate feather ID code (m)}
#' \item{l_vane}{Length of rachis/vane in the same row as the
#' appropriate feather ID code (m)}
#' \item{w_cal}{Width (diameter) of calamus in the same row as the
#' appropriate feather ID code (m)}
#' \item{w_vp}{Width of proximal vane (average value) in the same row as the
#' appropriate feather ID code (m)}
#' \item{w_vd}{Width of distal vane (average value)  in the same row as the
#' appropriate feather ID code (m)}
#' \item{vane_angle}{Interior angle between the rachis and calamus  (degrees)}
#' }
#'
"dat_feat_curr"

#' Material properties.
#'
#' A dataset containing the material properties for a
#' pigeon BirdID = 20_0300, TestID = 20_03004, FrameID = 317.
#'
#' @format A dataframe with 20 rows and 15 variables.
#' Only 8 variables required for proper functioning:
#' \describe{
#' \item{material}{Material information. Must include the following:
#' "Bone","Skin","Muscle","Cortex", "Medullary"}
#' \item{density}{Density of each material (kg/m^3)}
#' }
#'
"dat_mat"
