# Variogram Analysis in R for Engineers and Geoscientists New to R
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics Undergraduate Class 
# It is assumed that students have no previous R experience.  
# This utilizes the gstat library by Edzer Pedesma, appreciation to Dr. Pedesma for assistance.

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(gstat)                                 # geostatistical methods by Edzer Pebesma
library(sp)                                    # spatial points addition to regular data frames
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 

# Declare functions, this function completes standard normal transform on a dataset

nscore <- function(x) {                        # written by Ashton Shortridge, May/June, 2008
  # Takes a vector of values x and calculates their normal scores. Returns 
  # a list with the scores and an ordered table of original values and
  # scores, which is useful as a back-transform table. See backtr().
  nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/Users/pm27995/OneDrive - The University of Texas at Austin/Courses/PGE337_new/R/variogram")

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("dataset_101b.csv")          # read in comma delimited data file

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console

# The columns are variables with variable names at the top and the rows are samples

# Convert the dataframe to a spatial points dataframe
class(mydata)                                  # confirms that it is a dataframe
coordinates(mydata) = ~X+Y                     # indicate the X, Y spatial coordinates
summary(mydata)                                # confirms that it is now a spatial points dataframe
head(coordinates(mydata))                      # check the first several coordinates

# Normal scores transform of the porosity data
npor.trn = nscore(mydata$porosity)             # normal scores transform
mydata[["NPorosity"]]<-npor.trn$nscore         # append the normal scores transform into the spatial data table
head(mydata)

# Visualize the dataset with a bubble plot
par(mfrow=c(2,2))                              # set up a 2x2 matrix of plots 
# Now let's check the porosity data distribution
hist(mydata$porosity,main="Porosity (%)",xlab="Porosity (%)",nclass = 15) # hist builds a regular frequency histogram
# ecdf makes a cdf object and plot command plots it
plot(ecdf(mydata$porosity),main="Porosity",xlab="Porosity (%",ylab="Cumulative Probability")

# Now let's check the normal score transform of the porosity data distribution
hist(mydata$NPorosity,main="N[Porosity (%)]",xlab="N[Porosity (%)]",nclass = 15) # hist builds a regular frequency histogram
# ecdf makes a cdf object and plot command plots it
plot(ecdf(mydata$NPorosity),main="N[Porosity]",xlab="N[Porosity (%)]",ylab="Cumulative Probability")

# We can plot the porosity data location map
spplot(mydata, "porosity", do.log = TRUE,      # location map of porosity data
       key.space=list(x=0.75,y=0.97,corner=c(0,1)),
       scales=list(draw=T),xlab = "X (m)", ylab = "Y (m)",main ="Porosity (%)")

# We can also plot a simple bubble plot of the porosity data
bubble(mydata, "porosity", fill = FALSE, maxsize = 2, main ="Porosity (%)",identify = FALSE,xlab = "X (m)", ylab = "Y (m)")

# Calculate and fit an isotropic model to the normal score porosity variogram
por.vg.iso = variogram(NPorosity~1,mydata,cutoff = 3000,width =400,alpha = 35.0,tol.hor=45.0) # isotropic N[porosity]
plot(por.vg.iso,scaled = FALSE,var.lines = TRUE,main="Porosity Isotropic Variogram")

por.vm.iso <- vgm(psill = 0.6, "Exp", 800, anis = c(0, 0, 0, 1.0, 1.0),nugget=0.4)
por.vm.iso                                     # check the variogram model parameters

plot(por.vg.iso,por.vm.iso,main="Porosity Isotropic Variogram") # use the built in gstat variogram plot

# Let's calculate variogram maps to identify possible directionality

par(mfrow=c(2,2))                              # set up a 2x32 matrix of plots    
plot(variogram(NPorosity~1,mydata, cutoff=1500, width=400, map=TRUE),main = "Semivariogram Map",max=1.0)
plot(variogram(NPorosity~1,mydata, cutoff=1500, width=400, map=TRUE),main = "Number of Points",np=TRUE)

por.vg.035 = variogram(NPorosity~1,mydata,cutoff = 3000,width =500,alpha = 35.0,tol.hor=22.5)         # Calculate default isotropic variogram of N[Porosity]
por.vg.125 = variogram(NPorosity~1,mydata,cutoff = 3000,width =500,alpha = 125.0,tol.hor=22.5) 

plot(por.vg.035,main="Porosity Anisotropic 035 Variogram")
plot(por.vg.125,main="Porosity Anisotropic 125 Variogram")

# Anisotropy is parameterized as c(azimuth,dip,plunge,hratio,vratio) in #3D and c(azimuth,hratio) in 2D.
por.vm.ani <- vgm(psill = 0.6, "Exp", 800, anis = c(035, 0.5),nugget=0.4)
por.vm.ani                                     # check the variogram model parameters 

plot(por.vg.035,por.vm.ani,main="Porosity Anisotropic 035 Variogram") # use the built in gstat variogram plot
plot(por.vg.125,por.vm.ani,main="Porosity Anisotropic 125 Variogram") # use the built in gstat variogram plot

# Use auto-fitting to try to improve the variogram model, then check model and plot
por.vm.ani.auto <- fit.variogram(por.vg.iso,por.vm.ani,fit.sills = FALSE)
por.vm.ani.auto                                # check the autofit parameters and compare to our fit
plot(por.vg.035,por.vm.ani.auto,main="Porosity Anisotropic 035 Variogram")
plot(por.vg.125,por.vm.ani.auto,main="Porosity Anisotropic 125 Variogram")

# Let's make our own custom variogram plots
# The builtin variogram models do no fit into matrix plots, and are infelxible (e.g. no sill included)

par(mfrow=c(2,2))                              # set up a 2x32 matrix of plots 

name = c("Iso","035","125")                    # make name matrix
color = c("black","blue","red")                # make color matrix

# Isotropic plot

plot(por.vg.iso$dist,por.vg.iso$gamma,main="Porosity Isotropic Variogram",xlab="  Lag Distance (m) ",ylab=" Semivariogram ", col=color[1],ylim=c(0,1.2))
abline(h = 1.0)
lines(por.vm.iso$model)

unit_vector = c(1,0,0)                         # unit vector doesn't matter since isotropic
vm.iso <- variogramLine(por.vm.iso,maxdist=3000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate variogram model
lines(vm.iso$dist,vm.iso$gamma,col=color[1])   # include variogram model

legend(2000,.8,name, cex=0.8, col=color,pch=c(21,21,21),lty=c(1,1,1)) # add legend

# Anisotropic plot

plot(por.vg.035$dist,por.vg.035$gamma,main="Porosity Anisotropic Variogram",xlab="  Lag Distance (m) ",ylab=" Semivariogram ", col=color[2],ylim=c(0,1.2))
points(por.vg.125$dist,por.vg.125$gamma,col=color[3])
abline(h = 1.0)
lines(por.vm.iso$model)

unit_vector = c(sin(35*pi/180),cos(35*pi/180),0) # unit vector for 035 azimuth
vm.ani.035 <- variogramLine(por.vm.ani,maxdist=3000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 035 variogram model
lines(vm.ani.035$dist,vm.ani.035$gamma,col=color[2]) # include variogram model 

unit_vector = c(sin(55*pi/180),-1*cos(35*pi/180),0) # unit vector for 125 azimuth
vm.ani.125 <- variogramLine(por.vm.ani,maxdist=3000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 125 variogram model
lines(vm.ani.125$dist,vm.ani.125$gamma,col=color[3]) # include variogram model

legend(2000,.8,name, cex=0.8, col=color,pch=c(21,21,21),lty=c(1,1,1)) # add legend

