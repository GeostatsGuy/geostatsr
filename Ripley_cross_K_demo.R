# Cross and Regular K Functions in R to Access Spatial Point Patterns for Engineers and Geoscientists 

# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics undergraduate class and short courses. It is assumed that students 
# have no previous R experience. Interested in hosting and / or attending a short course contact me at mpyrcz@austin.utexas.edu.  

# This utilizes the spatstat library "Spatial Point Pattern Analysis, Model-Fitting, Simulation, Tests" by many authors.  The docs
# are at https://cran.r-project.org/web/packages/spatstat/spatstat.pdf.  

# For more information check out these references:

# Hajek E.A., Heller P.L., Sheets B.A. 2010. Significance of channel-belt clustering in alluvial basins. Geology, 38, 535-538 that 
# utilized Ripely K to access the stacking patterns of fluvial channels in a Ferris Formantion, South-central Wyoming outcrop.

# Jiang, Z, Shekhar, S., 2017, Spatial Big Data Science, Springer, p 131 a recent book on interesting methods for statistical analysis
# and learning (machine learning) with spatial data.

# To get started, let's load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first

library(spatstat)                              # Spatial statistical analysis by many authors
library(dplyr)                                 # A fast, consistent tool for working with data frame like objects by Hadley Wickam et al.

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/PGE337/Point")                           # choose your local working directory

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("2D_MV_200wells.csv")        # read in comma delimited data file

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console

# Check out the summary statistics for each column
summary(mydata)                                # summary statistics for the multivariate data file

# Now, lets vissualize the histograms of porosity and permeability
par(mfrow=c(1,2))                              # set up a 2x2 matrix of plots 
hist(mydata$porosity,main="Porosity",xlab="Porosity (%)",nclass = 15) # histogram
hist(mydata$permeability,main="Permeability",xlab="Permeability (mD)",nclass = 15) # histogram

# And now, let's see the location maps of porosity and permeability 
por.deciles <- quantile(mydata$porosity, 0:10/10)
cut.por    <- cut(mydata$porosity, por.deciles, include.lowest=TRUE)
plot(mydata$X,mydata$Y, col=grey(10:2/11)[cut.por], pch=20, xlab="X(m)",ylab="Y(m)",main="Porosity")

perm.deciles <- quantile(mydata$permeability, 0:10/10)
cut.perm    <- cut(mydata$permeability, perm.deciles, include.lowest=TRUE)
plot(mydata$X,mydata$Y, col=grey(10:2/11)[cut.perm], pch=20, xlab="X(m)",ylab="Y(m)",main="Permeability")

# Let's convert the wells into an unmarked point set for analysis
pp_unmarked <- ppp(mydata$X,mydata$Y,c(0,4000),c(0,4000), unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units

# Let's see what we have made
pp_unmarked
# A planar point pattern with 200 points in a box including the region X -> [0,4000] and y -> [0,4000] 

# Let's start by checking the regular Ripley's K of all the available sample data
K = Kest(pp_unmarked, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
K
plot(K)
# We use the "isotropic" boundary correction as it is best for rectangular windows (recal our study area is a simple 4000m x 4000m box) 
# We can observe that our well spatial point pattern is identifical to a random (Poisson) process over all isotropic spatial 
# length scales.  This is a reasonable result as the orginal synthetic data set is based on random sampling.

# Let's test out the cross-K function.

# We'll create a new indicator for porosity and permeability 
mydata <- cbind(mydata,mydata[,4])
names(mydata)[7]<-paste("PorInd")
mydata <- cbind(mydata,mydata[,5])
names(mydata)[8]<-paste("PermInd")

# Check to make sure we just made copies of the porosity and permeaiblity columns
head(mydata)

# Use nested conditional statements to make 3 levels of porosity (1 = low, 3 = mid, 2 = high) using the thresholds:

low_threshold = 0.14
high_threshold = 0.18
mydata$PorInd=ifelse(mydata$PorInd < low_threshold,1,ifelse(mydata$PorInd >= high_threshold,2,3))

# Lets look at the vector of the property that we just made
mydata$PorInd
# It should be a bunch of 1's and 2's

# Let's check out the distributions of our classes
par(mfrow=c(1,2))  
hist(mydata$PorInd,main="Porosity Classes",xlab="Porosity Class",nclass = 20, xlim = c(1,3), ylim = c(0,100)) # histogram

# Now let's make a new dataframe with middle values (3) removed, only retain the lows (1) and highs (2)  
filter = filter(mydata, PorInd == "1" | PorInd == "2" )

# We can confirm that we have removed the middle values (3) by looking at the DataFrame
filter

# Let's check out the new filtered distributions of our classes
hist(filter$PorInd,main="Porosity Classes",xlab="Porosity Class",nclass = 20, xlim = c(1,3), ylim = c(0,100)) # histogram

# Now we extract the PorInd column as a 1D vector, we will use this as our 2 factors for our point set
marker <- filter$PorInd

# Now we are ready to make a point set, with locations, bounding box and units
pp <- ppp(filter$X,filter$Y,c(0,4000),c(0,4000), unitname=c("metre","metres"))
# We can also add the marks, we convert our 1 / 2 values as doubles to factors to be recognized as levels 
pp <- pp %mark% factor(marker)

# Let's confirm that we have a planar point pattern with 109 points, multitype with levels 1 and 2
pp
# This can be confirmed in the output below

# We can plot our pointset with the levels indicated
par(mfrow=c(1,1))  
plot(pp, main = "Well Average Porosity, Low (1) and High (2)", xlab = "X(m)", ylab = "Y(m)")

# Now we can calculate the cross K function for the high values (2) relative to the low values (1)
K <- Kcross(pp, "1", "2",correction = "isotropic")

# We can use the built in plot feature to observe the result
plot(K)
# Over all scales we can see a lower density of high values relative to low values for porosity than predicted by
# random.  This suggestions some type of repulsion, or "anti-cross clustering".  This makes sense as the original
# data came from a simulation with spatial correlation (clusters or highs and lows).  The lows and high are repulsing
# eachother by design.  This occurs consistently acorss all scales observed.

# Now we can now reverse and calculate the cross K function for the low values (1) relative to the high values (1)
K <- Kcross(pp, "2", "1",correction = "isotropic")

# Again, we can use the built in plot feature to observe the result
plot(K)

# We get a very similar result.

# On your own try changing the boundary correction method and observe the results.  Also consider changing the low and high 
# porosity thresholds. 

# Hope this was helpful,

# Michael



