# Post-processing with Multiple Simulated Realizations Tutorial in R for Engineers and Geoscientists 
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics undergraduate class 
# It is assumed that students have no previous R experience.  
# This utilizes the gstat library by Edzer Pedesma, appreciation to Dr. Pedesma for assistance.

# Uncertainty in the subsurface is represented by multiple models.  We must work jontly with these multiple models to support decision making.  
# This post-processing example demonstrates a simple workflow to summarized over the scenarios and realizations of porosity.
# A summary:
# e-type is the expected value at each location over the models.  We can see where porosity tends to be high or low.
# local standard deviation is a measure of spread / uncertainty at each location in the model.
# local P10 is the P10 over the models at each location. High local P10 indicates surely high as 90% probability of even higher.
# local P90 is the P90 over the models at each location. Low local P90 indicates surely low as 90% probability of even lower.

# Note: We assume each model is equiprobable.  

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(gstat)                                 # geostatistical methods by Edzer Pebesma
library(sp)                                    # spatial points addition to regular data frames
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 
library(fields)                                # required for the image plots

# Specify the grid parameters (same as GSLIB / GEO-DAS parameterization in 2D)
nx = 100;ny = 100; xmn = 5.0; ymn = 5.0; xsize = 40.0; ysize = 40.0
xmin = xmn - xsize*0.5; ymin = ymn - ysize*0.5
xmax = xmin + nx * xsize; ymax = ymin + ny * ysize
x<-seq(xmin,xmax,by=xsize); y<-seq(ymin,ymax,by=ysize) # used for axes on image plots
                   
# Specify the parameters For plotting 
colmap = topo.colors(100)                      # define the color map and descritation

# Declare functions

# This function completes standard normal transform on a data vector
nscore <- function(x) {                        # by Ashton Shortridge, May/June, 2008
  # Takes a vector of values x and calculates their normal scores. Returns 
  # a list with the scores and an ordered table of original values and
  # scores, which is useful as a back-transform table. See backtr().
  nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

# This function builds a spatial points dataframe with the locations for estimation / simulation 
addcoord <- function(nx,xmin,xsize,ny,ymin,ysize) { # by Michael Pyrcz, March, 2018                      
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

sim2darray <- function(spdataframe,nx,ny,ireal) { # by Michael Pyrcz, March, 2018                      
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

sim2vector <- function(spdataframe,nx,ny,ireal) { # by Michael Pyrcz, March, 2018                      
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
mydata = mydata[1:50,]                         # take only the first 50 data (random order so random selection)

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console
# The columns are variables with variable names at the top and the rows are samples

# Convert the dataframe to a spatial points dataframe
class(mydata)                                  # confirms that it is a dataframe
coordinates(mydata) = ~X+Y                     # indicate the X, Y spatial coordinates
summary(mydata)                                # confirms that it is now a spatial points dataframe
head(coordinates(mydata))                      # check the first several coordinates

# Normal scores transform of the porosity data to assist with variogram calculation
npor.trn = nscore(mydata$porosity)             # normal scores transform
mydata[["NPorosity"]]<-npor.trn$nscore         # append the normal scores transform into the spatial data table
head(mydata)

# Make the model 2D grid
coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize)# make a dataframe with all the estimation locations
summary(coords)                                # check the coordinates

# Summary statistics of the property of interest
sill = var(mydata$porosity)
min = min(mydata$porosity)
max = max(mydata$porosity)
zlim = c(0,.30)                                # define the property min and max

# Three variogram scenarios
vm1 <- vgm(psill = 1.0*sill, "Sph", 600, anis = c(000, 1.0),nugget=0.0*sill)
vm1
vm2 <- vgm(psill = 1.0*sill, "Sph", 1200, anis = c(000, 1.0),nugget=0.0*sill)
vm2
vm3 <- vgm(psill = 0.99*sill, "Sph", 300, anis = c(000, 1.0),nugget=0.01*sill)
vm3

# Calculate the Gaussian simulations
nreal = 10; nscen = 3                          # 10 realizations of each 3 realizations
condsim1 = krige(porosity~1, mydata, coords, model = vm1, nmax = 100, nsim = nreal)
condsim2 = krige(porosity~1, mydata, coords, model = vm2, nmax = 100, nsim = nreal)
condsim3 = krige(porosity~1, mydata, coords, model = vm3, nmax = 100, nsim = nreal)
# nmax is the maximum number of data for each kriging solution
# this may take a long time to run, decrease nmax to speed up (e.g. 100)

# Plot four realizations for each scenario
par(mfrow=c(2,2))
real1 <- sim2darray(condsim1,nx,ny,1)      # extract realization #1 to a 2D array and plot
image.plot(real1,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #1", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real2 <- sim2darray(condsim1,nx,ny,2)      # extract realization #2 to a 2D array and plot
image.plot(real2,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #2", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real3 <- sim2darray(condsim1,nx,ny,3)      # extract realization #3 to a 2D array and plot
image.plot(real3,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #3", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real4 <- sim2darray(condsim1,nx,ny,4)      # extract realization #3 to a 2D array and plot
image.plot(real4,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #4", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

par(mfrow=c(2,2))
real1 <- sim2darray(condsim2,nx,ny,1)      # extract realization #1 to a 2D array and plot
image.plot(real1,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #1", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real2 <- sim2darray(condsim2,nx,ny,2)      # extract realization #2 to a 2D array and plot
image.plot(real2,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #2", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real3 <- sim2darray(condsim2,nx,ny,3)      # extract realization #3 to a 2D array and plot
image.plot(real3,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #3", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real4 <- sim2darray(condsim2,nx,ny,4)      # extract realization #3 to a 2D array and plot
image.plot(real4,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #4", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

par(mfrow=c(2,2))
real1 <- sim2darray(condsim3,nx,ny,1)      # extract realization #1 to a 2D array and plot
image.plot(real1,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #1", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real2 <- sim2darray(condsim3,nx,ny,2)      # extract realization #2 to a 2D array and plot
image.plot(real2,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #2", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real3 <- sim2darray(condsim3,nx,ny,3)      # extract realization #3 to a 2D array and plot
image.plot(real3,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #3", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

real4 <- sim2darray(condsim3,nx,ny,4)      # extract realization #3 to a 2D array and plot
image.plot(real4,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "Realization #4", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

#Put all the models in a single 3D array for convenience, indexed by ix, iy, imodel

nmodel = nreal*nscen
models = array(rep(0.0,nx*ny*nmodel),dim=c(nx,ny,nmodel)) # delare the 3D array
for(ireal in 1:nreal) {
  ixy = 1
  for(iy in 1:ny) {
    for(ix in 1:nx) {
      models[ix,iy,ireal] = condsim1@data[ixy,ireal]  
      models[ix,iy,ireal+nreal] = condsim2@data[ixy,ireal]
      models[ix,iy,ireal+2*nreal] = condsim3@data[ixy,ireal]
      ixy = ixy + 1
    }
  }
}

samples = vector('numeric', nmodel)

# Calculate the etype (expected value), and local standard deviation, P10 and P90 at each location over all models

etype = matrix(nrow = nx,ncol = ny)
stdev = matrix(nrow = nx,ncol = ny)
p10 = matrix(nrow = nx,ncol = ny)
p90 = matrix(nrow = nx,ncol = ny)

for(iy in 1:ny) {
  for(ix in 1:nx) {
    aver = 0.0
    count = 0
    for(imodel in 1:nmodel) {
      samples[imodel] = models[ix,iy,imodel]
    }
    etype[ix,iy] = mean(samples)
    stdev[ix,iy] = sd(samples)
    p10[ix,iy] = quantile(samples,0.10)
    p90[ix,iy] = quantile(samples,0.90)
  }
}

# Plot the post-processed summaries 

par(mfrow=c(2,2))

image.plot(etype,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "etype", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

image.plot(stdev,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = c(0,0.05),col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "conditional standard deviation", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

image.plot(p10,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "local P10", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

image.plot(p90,x=x,y=y,xlab="X(m)",ylab="Y(m)",zlim = zlim,col=colmap,legend.shrink = 0.6); mtext(line=1, side=3, "local P00", outer=F);box(which="plot")
points(mydata$X,mydata$Y,pch="+",cex=1.0,col="black") # add well locations

# There is so much we can do to explore our models of uncertainty.  This post-processing approach is very simple, but powerful.

# Hope this was helpful,

# Michael

# Appreciation to Edzer Pedesma for the gstat package.

