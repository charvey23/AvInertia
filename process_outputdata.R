# This script was written to analyse the data that is returned from the birdmoment package

# ------------- Check how the elbow and wrist affects the mass properties of the wing -----------------
clean_dat <- subset(alldat, sweep == 0 & dihedral == 0)[,c(3,4,5,44,45,46)]
full_results <- merge(all_data, clean_dat, all.x = TRUE, by.x = c("species","WingID","TestID","FrameID"), by.y = c("species","WingID","TestID","frameID"))

model_Ixx <- lme4::lmer(value ~ elbow + manus + (1|WingID), data = subset(full_results, component == "wing" & object == "Ixx"))

plot(subset(full_results, component == "wing" & object == "Ixx")$elbow,subset(full_results, component == "wing" & object == "Ixx")$value)
