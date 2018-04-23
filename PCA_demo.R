# Principle Component Analysis in R for Engineers and Geoscientists New to R
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics Undergraduate Class 
# It is assumed that students have no previous R experience.  We take an multivariate dataset
# and only retain the two variables for a simple demonstration of PCA.

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 
library(ggplot2)                               # for the custom biplot
library(lattice)                               # for the matrix scatter plot
library(corrplot)                              # for the corrplot correlation plot

# Declare functions
# no functions required in this demonstration

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/PGE337/PCA")                         # choose your local working directory

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("unconv_MV.csv")             # read in comma delimited data file

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console

# Check out the summary statistics for each column
summary(mydata)                                # summary statistics for the multivariate data file

# Calculate the correlation matrix 
mydata_noindex <- mydata[,2:length(mydata)]    # remove the first column with the well index
cor_matrix <- round(cor(mydata_noindex),2)     # calculate a mxm matrix with the correlation coeficients
cor_matrix

# Let's use the corrplot package to make a very nice correlation matrix visualization
corrplot(cor_matrix, method = "circle")        # graphical correlation matrix plot

# Now let's view the scatterplot matrices from the lattice Package 
splom(mydata[c(2,3,4,5,6,7,8)],col=rgb(0,0,0,50,maxColorValue=255), pch=19,main = "Unconventional Dataset") 
# This dataset has varriables from 1,000 unconventional wells including well average porosity, log transform 
# of permeability (to linearize the relationships with other variables), accoustic impedance (kg/m2s*106), brittness ratio (%),
# total organic carbon (%), vitrinite reflectance (%), and production (MCFPD)

# Let's start simple with a bivariate (2 variable) problem
mydata_por_perm <- mydata[1:100,2:3]           # new dataframe with only 1st 100 wells of por and logperm
head(mydata_por_perm)                          # check the new dataframe

# Look at a scatter plot of porosity vs. log permeability
plot(LogPerm~Por, mydata_por_perm, main="Log Permeability vs. Porosity", 
     xlab="Porosity (%)", ylab="Log Permeability", col = alpha("black",0.3), pch=19)
# With the log of permeability we have a very nice linear relationship with porosity

# We are ready to perform PCA with porosity and log of permeability
pca <- prcomp(mydata_por_perm, scale=TRUE)     # this does the PCA
# Note, we should scale the data to all have a standard deviation of 1.0.  Otherwise the difference
# between the scale of porosity and permeability would have a significant impact.  We should always
# scale unless the two variables have the same units.

# Let's see what's in the PCA output
names(pca)
# This includes:
# "sdev" is the scaling applied to get all the standard deviations to 1.0 (m)
# "center" is the means that we subtracted from each variable to get them centered (m)
# "rotation" includes the component / factor loadings (mxm)
# "x" includes the principle component scores for all the data x all the principal components (nxm)

# let's look at the average and scaling factor applied to the data before PCA
pca$center                                     # the means substracted from variables
pca$scale                                      # the standarization factor
# These are m sized arrays as they are distinct for each feature

# Now, let's look at the principle component loadings are in the rotation matrix
pca$rotation

# the 'x' matrix provides us with the columns with the kth principal component score vector
dim(pca$x)
# no suprise there is number of wells columns each nvariables long 
# note the porosity and permeabilities are centered first
head(pca$x)

# Let's go head and reverse the PCA with with 2 princple components to check
nComp = 2
Xhat2 = t(t(pca$x[,1:nComp] %*% t(pca$rotation[,1:nComp])) * pca$scale + pca$center)

# Let's plot the original data vs. our back transformation using both PC1 and PC2 
# We can confirm that we got back to the original data; therefore, our method works.
par(mfrow=c(2,1)) 
plot(LogPerm~Por, mydata_por_perm, main="Log Permeability vs. Porosity", 
     xlab="Porosity (%)", ylab="Log Permeability", xlim = c(8,22), ylim = c(0,3),col = alpha("black",0.3), pch=19)
plot(LogPerm~Por, Xhat2, main="LogPerm vs. Porsity (from PC1 and PC2)", 
     xlab="Porosity (%)", ylab="LogPerm", xlim = c(8,22), ylim = c(0,3),col = alpha("black",0.3), pch=19)
# The two cross plots should look exactly the same. If so the method is working.

# Now let's attempt dimensional reduction by only retaining the first principle component
nComp = 1                                  # we only retain PC1
Xhat1 = t(t(pca$x[,1:nComp] %*% t(pca$rotation[,1:nComp])) * pca$scale + pca$center)
# note "%*%" indicate matrix multiplication, otherwise "*" and "+" multiple and add over the columns and
# "t" is the transpose of the matrix, aside matrix math is really easy in R, eh!

# Let's plot original data, principal coordinate scores for 2 and 1 component and backtransform.
par(mfrow=c(2,2)) 
plot(LogPerm~Por, mydata_por_perm, main="Log Permeability vs. Porosity", 
     xlab="Porosity (%)", ylab="Log Permeability", xlim = c(8,22), ylim = c(0,3),col = alpha("black",0.3), pch=19)
plot(LogPerm~Por, Xhat1, main="LogPerm vs. Porsity (from PC1 Only)", 
     xlab="Porosity (%)", ylab="LogPerm", xlim = c(8,22), ylim = c(0,3),col = alpha("black",0.3), pch=19)
plot(PC2~PC1, (pca$x*-1), main="PC2 vs. PC1", 
     xlab="PC1", ylab="PC2", xlim = c(-4,4), ylim = c(-1,1), col = alpha("black",0.3), pch=19)
plot(rep(0,100)~PC1, (pca$x*-1), main="PC1 Only", 
     xlab="PC1", ylab="", xlim = c(-4,4), ylim = c(-1,1), col = alpha("black",0.3), pch=19)
# With one principle component we can describe quite a bit of the variance!

# it is useful to look at the biplot
par(mfrow=c(1,1)) 
biplot(pca,cex = 0.8, scale=0, pch = 21)
# This plot shows the data component scores ploted along with the feature loadings values indicated on right and top axes.

# It may be useful to just look at the factor loadings
theta <- seq(0,2*pi,length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle,aes(x,y)) + geom_path()
loadings <- data.frame(pca$rotation,.names = row.names(pca$rotation))
p + geom_text(data=loadings,mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
# This is a custom ggplot for veiwing the associations of the features relative to each principal
# component. Note the each feature is balanced between the two principal coordinates.

# Calculate and visualize the variance explained
pr.var=pca$sdev^2
pve=pr.var/sum(pr.var)
pve
# new we know we can explain about 90% of the variance with a single principal component

# Let's plot the variance explained by each principle component
plot(pve , xlab="Principal Component ", ylab="Proportion of Variance Explained ", xlim = c(1,5),ylim=c(0,1),pch=19,col="red")
plot(cumsum(pve), xlab="Principal Component ", ylab=" Cumulative Proportion of Variance Explained ",xlim = c(1,5), ylim=c(0,1),pch = 19,col="red")

# The pca object has a built-in variance by principal component plot
plot(pca, type = "l", main = "Variance by Principal Component")

# Also, we can also summarize the PCA model like this 
summary(pca)

# On your own try working with more variables.  How few principal components can we use? How much information 
# do we lose?

# Hope this was helpful,

# Michael
#
# Michael J. Pyrcz, P.Eng., Ph.D.
# Associate Professor, the University of Texas at Austin (@GeostatsGuy)
# mpyrcz@austin.utexas.edu


