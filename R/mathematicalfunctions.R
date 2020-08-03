# Basic mathematical functions

calc_rot(z_vector, x_vector){
  #CAUTION: incoming vectors must be in the structural frame at the VRP
  #         x and z vectors must already be orthogonal axes

  unit_z_vector = z_vector/norm(z_vector)
  unit_x_vector = x_vector/norm(x_vector)

  #cross z with x to get the righthanded axis
  y_vector = cross(unit_z_vector,unit_x_vector)
  y_vector = y_vector/norm(y_vector)

  VRP2object = rbind(x_vector,y_vector,z_vector)
  object2VRP = t(structural2in)

  return()
}
