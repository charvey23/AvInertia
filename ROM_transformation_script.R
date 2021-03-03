##### This scripts reads in all range of motion tracks for a specific species.
##### Resizes down to each individual, subsamples and finally reorients each wing shape

## Load required libraries
library(tidyverse)
library(Morpho)
library(geomorph)
library(stringr)

source('ROM_transformation_functions.R')


#-------- Find all file names necessary -------
working_path  = "/Users/christinaharvey/Dropbox (University of Michigan)/Bird Mass Distribution/03_processed_optitrack/"
filename.3d   = list.files(path = working_path, pattern = paste('*',".csv",sep=""))
max_sample_no = 2 # maximum amount of samples to have within one bin
bin_size      = 2 # #deg x #deg bins that will have the max amount of samples

# create black data frame
dat_ID = as.data.frame(matrix(nrow = length(filename.3d), ncol = 4))
colnames(dat_ID) = c("species","BirdID","TestID","wing_side")

# Save all the identifying information
for (i in 1:length(filename.3d)){
  dat_ID$species[i]    = paste(tolower(strsplit(filename.3d[i], "_")[[1]][2]),"_",strsplit(filename.3d[i], "_")[[1]][3], sep = "")
  dat_ID$BirdID[i]     = paste(strsplit(filename.3d[i], "_")[[1]][4],"_",str_sub(strsplit(filename.3d[i], "_")[[1]][5],start = 1, end = -2), sep = "")
  dat_ID$TestID[i]     = as.numeric(strsplit(filename.3d[i], "_")[[1]][7])
  dat_ID$wing_side[i]  = str_sub(strsplit(filename.3d[i], "_")[[1]][5],-1)
}

no_species = length(unique(dat_ID$species))

for (i in 1:no_species){
 files_to_read = which(dat_ID$species == unique(dat_ID$species)[i])

 # Read in all optitrack data from the given species
 raw_dat <- read.csv(filename.3d[files_to_read[1]], stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA"))
 for (j in 1:length(files_to_read)){
   raw_dat <- rbind(raw_dat, read.csv(filename.3d[i], stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA")))
 }



}

