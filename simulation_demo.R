# Nested Variogram Tutorial in R for Engineers and Geoscientists 
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics undergraduate class 
# It is assumed that students have no previous R experience.  
# This utilizes the gstat library by Edzer Pedesma, appreciation to Dr. Pedesma for assistance.

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(gstat)                                 # geostatistical methods by Edzer Pebesma
library(sp)                                    # spatial points addition to regular data frames
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 

# Specify the grid parameters (same as GSLIB / GEO-DAS parameterization in 2D)
nx = 100
ny = 100
xmn = 5.0
ymn = 5.0
xsize = 10.0
ysize = 10.0
xmin = xmn - xsize*0.5
ymin = ymn - ysize*0.5
xmax = xmin + nx * xsize
ymax = ymin + ny * ysize
x<-seq(xmin,xmax,by=xsize)
y<-seq(ymin,ymax,by=ysize)

# Specify the parameters For plotting 
colmap = topo.colors(100)                      # define the color map and descritation
zlim = c(-3,3)                                 # define the property min and max

# Declare functions

# This function completes standard normal transform on a data vector

nscore <- function(x) {                        # written by Ashton Shortridge, May/June, 2008
  # Takes a vector of values x and calculates their normal scores. Returns 
  # a list with the scores and an ordered table of original values and
  # scores, which is useful as a back-transform table. See backtr().
  nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

# This function builds a spatial points dataframe with the locations for estimation / simulation 

addcoord <- function(nx,xmin,xsize,ny,ymin,ysize) { # Michael Pyrcz, March, 2018                      
  # makes a 2D dataframe with coordinates based on GSLIB specification
  coords = matrix(nrow = nx*ny,ncol=2)
  ixy = 1
  for(iy in 1:nx) {
    for(ix in 1:ny) {
      coords[ixy,1] = xmin + (ix-1)*xsize  
      coords[ixy,2] = ymin + (iy-1)*ysize 
      ixy = ixy + 1
    }
  }
  coords.df = data.frame(coords)
  colnames(coords.df) <- c("X","Y")
  coordinates(coords.df) =~X+Y
  return (coords.df)
}  

sim2darray <- function(spdataframe,nx,ny,ireal) { # Michael Pyrcz, March, 2018                      
  # makes a 2D array from realizations spatial point dataframe
  model = matrix(nrow = nx,ncol = ny)
  ixy = 1
  for(iy in 1:ny) {
    for(ix in 1:nx) {
      model[ix,iy] = spdataframe@data[ixy,ireal]  
      ixy = ixy + 1
    }
  }
  return (model)
}  

sim2vector <- function(spdataframe,nx,ny,ireal) { # Michael Pyrcz, March, 2018                      
  # makes a 1D vector from spatial point dataframe
  model = rep(0,nx*ny)
  ixy = 1
  for(iy in 1:ny) {
    for(ix in 1:nx) {
      model[ixy] = spdataframe@data[ixy,ireal]  
      ixy = ixy + 1
    }
  }
  return (model)
} 

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/PGE337")

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("2D_MV_200Wells.csv")        # read in comma delimited data file

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console
# The columns are variables with variable names at the top and the rows are samples

# Convert the dataframe to a spatial points dataframe
class(mydata)                                  # confirms that it is a dataframe
coordinates(mydata) = ~X+Y                     # indicate the X, Y spatial coordinates
summary(mydata)                                # confirms that it is now a spatial points dataframe
head(coordinates(mydata))                      # check the first several coordinates

# Normal scores transform of the porosity data to assist with variogram calculation
npor.trn = nscore(mydata$por)                  # normal scores transform
mydata[["NPorosity"]]<-npor.trn$nscore         # append the normal scores transform into the spatial data table
head(mydata)

# Make the model 2D grid
coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize) # make a dataframe with all the estimation locations
summary(coords)                                # check the coordinates

# Anisotropic variogram
vm.nug1 <- vgm(psill = 0.5, "Sph", 400, anis = c(000, 1.0),nugget=0.5)
vm.nug1
vm.ani1 <- vgm(psill = 1.0, "Exp", 200, anis = c(035, 0.5),nugget=0.0)
vm.ani1
vm.ani2 <- vgm(psill = 1.0, "Sph", 600, anis = c(060, 0.2),nugget=0.0)
vm.ani2

# Gaussian simulation

condsim.nug1 = krige(por~1, mydata, coords, model = vm.nug1, nmax = 500, nsim = 4)
condsim.ani1 = krige(por~1, mydata, coords, model = vm.ani1, namx = 500, nsim = 4)
condsim.ani2 = krige(por~1, mydata, coords, model = vm.ani2, nmax = 500, nsim = 4)
# nmax is the maximum number of data for each kriging solution
# this may take a long time to run, decrease nmax to spped up (e.g. 100)

# Plot all the realizations
par(mfrow=c(2,2))
real1 <- sim2darray(condsim.nug1,nx,ny,1)      # extract a realization and plot
image.plot(real1,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #1", outer=F);box(which="plot")

real2 <- sim2darray(condsim.nug1,nx,ny,2)
image.plot(real2,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #2", outer=F);box(which="plot")

real3 <- sim2darray(condsim.nug1,nx,ny,3)
image.plot(real3,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #3", outer=F);box(which="plot")

real4 <- sim2darray(condsim.nug1,nx,ny,4)
image.plot(real4,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #4", outer=F);box(which="plot")

par(mfrow=c(2,2))
real1 <- sim2darray(condsim.ani1,nx,ny,1)
image.plot(real1,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #1", outer=F);box(which="plot")

real2 <- sim2darray(condsim.ani1,nx,ny,2)
image.plot(real2,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #2", outer=F);box(which="plot")

real3 <- sim2darray(condsim.ani1,nx,ny,3)
image.plot(real3,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #3", outer=F);box(which="plot")

real4 <- sim2darray(condsim.ani1,nx,ny,4)
image.plot(real4,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #4", outer=F);box(which="plot")

par(mfrow=c(2,2))
real1 <- sim2darray(condsim.ani2,nx,ny,1)
image.plot(real1,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #1", outer=F);box(which="plot")

real2 <- sim2darray(condsim.ani2,nx,ny,2)
image.plot(real2,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #2", outer=F);box(which="plot")

real3 <- sim2darray(condsim.ani2,nx,ny,3)
image.plot(real3,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #3", outer=F);box(which="plot")

real4 <- sim2darray(condsim.ani2,nx,ny,4)
image.plot(real4,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #4", outer=F);box(which="plot")

par(mfrow=c(3,2))

# Minimum acceptance checks, check the distributions and variagrams of the realization

# High Nugget Effect
plot(ecdf(condsim.nug1@data[,1]),main="High Nugget",xlab="Gaussian Values",ylab="Cumulative Probability",col="red")
plot(ecdf(condsim.nug1@data[,2]),add=TRUE,col="red")
plot(ecdf(condsim.nug1@data[,3]),add=TRUE,col="red")
plot(ecdf(condsim.nug1@data[,4]),add=TRUE,col="red")

vg.sim1.035 = variogram(sim1~1,condsim.nug1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim1.125 = variogram(sim1~1,condsim.nug1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 
vg.sim2.035 = variogram(sim2~1,condsim.nug1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim2.125 = variogram(sim2~1,condsim.nug1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 
vg.sim3.035 = variogram(sim3~1,condsim.nug1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim3.125 = variogram(sim3~1,condsim.nug1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 
vg.sim4.035 = variogram(sim4~1,condsim.nug1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim4.125 = variogram(sim4~1,condsim.nug1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 

plot(vg.sim1.035$dist,vg.sim1.035$gamma,pch=19,cex=0.1,main="High Nugget",xlab="  Lag Distance (m) ",ylab=" Semivariogram ", col="black",xlim=c(0,1000),ylim=c(0,1.2))
points(vg.sim1.125$dist,vg.sim1.125$gamma,pch=19,cex=0.1,col="red")
points(vg.sim2.035$dist,vg.sim2.035$gamma,pch=19,cex=0.1,col="black")
points(vg.sim2.125$dist,vg.sim2.125$gamma,pch=19,cex=0.1,col="red")
points(vg.sim3.035$dist,vg.sim3.035$gamma,pch=19,cex=0.1,col="black")
points(vg.sim3.125$dist,vg.sim3.125$gamma,pch=19,cex=0.1,col="red")
points(vg.sim4.035$dist,vg.sim4.035$gamma,pch=19,cex=0.1,col="black")
points(vg.sim4.125$dist,vg.sim4.125$gamma,pch=19,cex=0.1,col="red")
abline(h = 1.0)
# Include variogram model
unit_vector = c(0,1,0)                         # unit vector for 000 azimuth
vm.nug1.000 <- variogramLine(vm.nug1,maxdist=1000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 035 variogram model
lines(vm.nug1.000$dist,vm.nug1.000$gamma,pch=19,cex=0.1,col="black") # include variogram model 

unit_vector = c(1,0,0)                         # unit vector for 090 azimuth
vm.ani.090 <- variogramLine(vm.nug1,maxdist=1000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 125 variogram model
lines(vm.ani.090$dist,vm.ani.090$gamma,col="red") 

# Low Anisotropy
plot(ecdf(condsim.ani1@data[,1]),main="Low Aniostropic",xlab="Gaussian Values",ylab="Cumulative Probability",col="red")
plot(ecdf(condsim.ani1@data[,2]),add=TRUE,col="red")
plot(ecdf(condsim.ani1@data[,3]),add=TRUE,col="red")
plot(ecdf(condsim.ani1@data[,4]),add=TRUE,col="red")

vg.sim1.035 = variogram(sim1~1,condsim.ani1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim1.125 = variogram(sim1~1,condsim.ani1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 
vg.sim2.035 = variogram(sim2~1,condsim.ani1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim2.125 = variogram(sim2~1,condsim.ani1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 
vg.sim3.035 = variogram(sim3~1,condsim.ani1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim3.125 = variogram(sim3~1,condsim.ani1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 
vg.sim4.035 = variogram(sim4~1,condsim.ani1,cutoff = 1000,width =20,alpha = 35.0,tol.hor=22.5) 
vg.sim4.125 = variogram(sim4~1,condsim.ani1,cutoff = 1000,width =20,alpha = 125.0,tol.hor=22.5) 

plot(vg.sim1.035$dist,vg.sim1.035$gamma,pch=19,cex=0.1,main="Low Anisotropic",xlab="  Lag Distance (m) ",ylab=" Semivariogram ", col="black",xlim=c(0,1000),ylim=c(0,1.2))
points(vg.sim1.125$dist,vg.sim1.125$gamma,pch=19,cex=0.1,col="red")
points(vg.sim2.035$dist,vg.sim2.035$gamma,pch=19,cex=0.1,col="black")
points(vg.sim2.125$dist,vg.sim2.125$gamma,pch=19,cex=0.1,col="red")
points(vg.sim3.035$dist,vg.sim3.035$gamma,pch=19,cex=0.1,col="black")
points(vg.sim3.125$dist,vg.sim3.125$gamma,pch=19,cex=0.1,col="red")
points(vg.sim4.035$dist,vg.sim4.035$gamma,pch=19,cex=0.1,col="black")
points(vg.sim4.125$dist,vg.sim4.125$gamma,pch=19,cex=0.1,col="red")
abline(h = 1.0)
# Include variogram model
unit_vector = c(sin(35*pi/180),cos(35*pi/180),0) # unit vector for 035 azimuth
vm.nug1.035 <- variogramLine(vm.ani1,maxdist=1000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 035 variogram model
lines(vm.nug1.035$dist,vm.nug1.035$gamma,pch=19,cex=0.1,col="black") # include variogram model 

unit_vector = c(sin(55*pi/180),-1*cos(35*pi/180),0) # unit vector for 125 azimuth
vm.ani.125 <- variogramLine(vm.ani1,maxdist=1000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 125 variogram model
lines(vm.ani.125$dist,vm.ani.125$gamma,col="red") # include variogram model

# High Anisotropy
plot(ecdf(condsim.ani2@data[,1]),main="High Aniostropic",xlab="Gaussian Values",ylab="Cumulative Probability",col="red")
plot(ecdf(condsim.ani2@data[,2]),add=TRUE,col="red")
plot(ecdf(condsim.ani2@data[,3]),add=TRUE,col="red")
plot(ecdf(condsim.ani2@data[,4]),add=TRUE,col="red")

vg.sim1.060 = variogram(sim1~1,condsim.ani2,cutoff = 1000,width =20,alpha = 60.0,tol.hor=22.5) 
vg.sim1.150 = variogram(sim1~1,condsim.ani2,cutoff = 1000,width =20,alpha = 150.0,tol.hor=22.5) 
vg.sim2.060 = variogram(sim2~1,condsim.ani2,cutoff = 1000,width =20,alpha = 60.0,tol.hor=22.5) 
vg.sim2.150 = variogram(sim2~1,condsim.ani2,cutoff = 1000,width =20,alpha = 150.0,tol.hor=22.5) 
vg.sim3.060 = variogram(sim3~1,condsim.ani2,cutoff = 1000,width =20,alpha = 60.0,tol.hor=22.5) 
vg.sim3.150 = variogram(sim3~1,condsim.ani2,cutoff = 1000,width =20,alpha = 150.0,tol.hor=22.5) 
vg.sim4.060 = variogram(sim4~1,condsim.ani2,cutoff = 1000,width =20,alpha = 60.0,tol.hor=22.5) 
vg.sim4.150 = variogram(sim4~1,condsim.ani2,cutoff = 1000,width =20,alpha = 150.0,tol.hor=22.5) 

plot(vg.sim1.060$dist,vg.sim1.060$gamma,pch=19,cex=0.1,main="High Anisotropic",xlab="  Lag Distance (m) ",ylab=" Semivariogram ", col="black",xlim=c(0,1000),ylim=c(0,1.2))
points(vg.sim1.150$dist,vg.sim1.150$gamma,pch=19,cex=0.1,col="red")
points(vg.sim2.060$dist,vg.sim2.060$gamma,pch=19,cex=0.1,col="black")
points(vg.sim2.150$dist,vg.sim2.150$gamma,pch=19,cex=0.1,col="red")
points(vg.sim3.060$dist,vg.sim3.060$gamma,pch=19,cex=0.1,col="black")
points(vg.sim3.150$dist,vg.sim3.150$gamma,pch=19,cex=0.1,col="red")
points(vg.sim4.060$dist,vg.sim4.060$gamma,pch=19,cex=0.1,col="black")
points(vg.sim4.150$dist,vg.sim4.150$gamma,pch=19,cex=0.1,col="red")
abline(h = 1.0)
# Include variogram model
unit_vector = c(sin(60*pi/180),cos(60*pi/180),0) # unit vector for 060 azimuth
vm.nug1.060 <- variogramLine(vm.ani2,maxdist=1000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 035 variogram model
lines(vm.nug1.060$dist,vm.nug1.060$gamma,pch=19,cex=0.1,col="black") # include variogram model 

unit_vector = c(sin(30*pi/180),-1*cos(60*pi/180),0) # unit vector for 125 azimuth
vm.ani.150 <- variogramLine(vm.ani2,maxdist=1000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 125 variogram model
lines(vm.ani.150$dist,vm.ani.150$gamma,col="red") # include variogram model

# On your own try changing the variogram parameters and observe the results.  Also consider kriging with a trend.

# Hope this was helpful,

# Michael

# Appreciation to Edzer Pedesma for the gstat package.
