# Try out the Hurst coefficient and compare to autocorrelation.
# Michael Pyrcz, the University of Texas at Austin, @GeostatsGuy
# Datafile available in the GeoDataSets repository.

# load required libraries
library(pracma)
library(zoo)

# Set the working directory
setwd("C:/Users/pm27995/OneDrive - The University of Texas at Austin/Courses/PGE337_new/R/Bivariate")

# Load the data
mydata = read.csv("corr_formB_realizations.csv") # read csv file
mydata

# Sampling Array in Depth
Depth <- (1:100)
depth_int_df <- data.frame(Depth)

# Extract fines and depth vectors
Fines <- mydata$V1

# Dataset with iterative smoothing and Hurst calculation
par(mfrow=c(4,1)) 

Fines <- (Fines - mean(Fines))/sd(Fines)

plot(Fines,type="l",xlim=c(0,1000),ylim=c(-3.0,3.0))
lines(Fines,xlim=c(0,1000))
hurstexp(Fines,d=50,display = TRUE)

smooth <- rollapply(Fines, width = 8, by = 1, mean)
smooth <- (smooth - mean(smooth))/sd(smooth)
plot(smooth,type="l",xlim=c(0,1000), ylim=c(-3.,3))
lines(smooth,xlim=c(0,1000))
hurstexp(smooth,d=50,display = TRUE)

smooth2 <- rollapply(smooth, width = 8, by = 1, mean)
plot(smooth2,type="l",xlim=c(0,1000),ylim=c(-3,3))
lines(smooth2,xlim=c(0,1000))
hurstexp(smooth2,d=50,display = TRUE)

smooth3 <- rollapply(smooth2, width = 8, by = 1, mean)
plot(smooth3,type="l",xlim=c(0,1000),ylim=c(-3,3))
lines(smooth3,xlim=c(0,1000))
hurstexp(smooth3,d=50,display = TRUE)

# Autocorrelation wit signficance plots.
par(mfrow=c(4,1))

acf(Fines)

acf(smooth)

acf(smooth2)

acf(smooth3)

