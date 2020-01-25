
# Decision Trees for Engineers and Geoscientists New to R
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics Undergraduate Class 
# It is assumed that students have no previous R experience.  We use an 
# multivariate dataset and only retain the three variables for a simple 
# demonstration of decision trees.

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 
library(ggplot2)                               # for the custom biplot
library(lattice)                               # for the matrix scatter plot
library(corrplot)                              # for the corrplot correlation plot
library(tree)                                  # for decision tree

# Declare functions
# no functions required in this demonstration

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("D:/PGE383")                          # choose your local working directory

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("unconv_MV_v2.csv")             # read in comma delimited data file

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
# of permeability (to linearize the relationships with other variables), accoustic impedance (kg/m2s*10^6), brittness ratio (%),
# total organic carbon (%), vitrinite reflectance (%), and production (MCFPD)

# Let's start simple with a trivariate (3 variable) problem
mydata_por <- data.frame(mydata[1:1000,2])      # extract and rename 3 features from the original dataframe
colnames(mydata_por) <- "Por"
mydata_brittle <- data.frame(mydata[1:1000,5])   
colnames(mydata_brittle) <- "Brittle"
mydata_prod <- data.frame(mydata[1:1000,9])
colnames(mydata_prod) <- "Prod" 
mydata_3var <- cbind(mydata_por,mydata_brittle,mydata_prod)
head(mydata_3var)                               # check the new dataframe

par(mfrow=c(1,2))                               # check out the production vs. porosity and brittleness
plot(mydata_3var$Por,mydata_3var$Prod, main="Production vs. Porosity", 
     xlab="Porosity (%)", ylab="Production (MCFPD)", col = alpha("black",0.1), pch=19)
plot(mydata_3var$Brittle,mydata_3var$Prod, main="Production vs. Brittleness", 
     xlab="Brittleness (%)", ylab="Production (MCFPD)", col = alpha("black",0.1), pch=19)

#Let's plot porosity vs. brittleness with production as greyscale
par(mfrow=c(1,2))                               # we save space for the training data set next 
prod.deciles <- quantile(mydata_3var$Prod, 0:10/10)
cut.prod    <- cut(mydata_3var$Prod, prod.deciles, include.lowest=TRUE)
plot(mydata_3var$Por,mydata_3var$Brittle, col=grey(10:2/11)[cut.prod], pch=20, xlab="Porosity (%)",ylab="Brittleness (%)",main="Production (MCFPD)")

# We will use random processes, to ensure repeatability between runs let's set the random seed
set.seed(71071)

# Let's extract a subset of the data to use as training data
train = mydata_3var[sample(nrow(mydata_3var),500,replace = FALSE),]
train_ind <- sample(seq_len(nrow(mydata_3var)), size = 500, replace = FALSE)
train <- mydata_3var[train_ind, ]
test <- mydata_3var[-train_ind, ]

# Let's check the training data to make sure it has the correct format
head(train)                                     # note the indexes are randomized

# Let's plot the training data set and compare to the original data set.
cut.train.prod <- cut(train$Prod, prod.deciles, include.lowest=TRUE)
plot(train$Por,train$Brittle, col=grey(10:2/11)[cut.train.prod], pch=20, xlab="Porosity (%)",ylab="Brittleness (%)",main="Production (MCFPD)")
# With the density and coverage of training data this will not be a difficult prediction exercise
# There is no signiifcant extrapolation

# We can design the controls on the tree growth 
tree.control = tree.control(nobs = 500, mincut = 5, minsize = 10, mindev = 0.01)
# nobs is the number of data in training set, mincut / minsize are minimum node size constraints 
# and mindev is the minimum deviation in a node to allow a split
# These of the defaults in the package. You can change these later and rerun to observe increased or 
# decreased complexity.

# We are ready to run the decision tree on our training data 
tree.prod = tree(Prod~Por+Brittle,train,control = tree.control)

# We can apply the summary command to get some basic information on our tree
summary(tree.prod)                              # note complexity in number of terminal nodes / regions

# The advantage of a decision tree is that it may be view graphically
plot(tree.prod)                                 # plots the decision tree
text(tree.prod,pretty=0)                        # adds the decision rules
# Note: the branch lengths are porportional to the decrease in impurity

# Another good way to visualize our tree is to review the regions and estimates
plot(mydata_3var$Por,mydata_3var$Brittle, col=grey(10:2/11)[cut.prod], pch=20, xlab="Porosity (%)",ylab="Brittleness (%)")
partition.tree(tree.prod, ordvars=c("Por","Brittle"), add=TRUE)

# Let's evaluate if pruning will improve our model
cv.prod = cv.tree(tree.prod,K = 10)             # this runs the k fold cross validation and report RSS
plot(cv.prod$size,cv.prod$dev,type='b')
# K is the number of folds in the cross validation
# This provides the deviance or number of misclassifications as a function of the number of terminal nodes

# We could decide to prune our tree to 6 nodes as we get little improvement with subsequent add in complexity
prune.prod = prune.tree(tree.prod,best = 6)     # we reduce the complexity of the tree to 6 terminal nodes

# Let's look at the new tree side-by-side with the unpruned tree
par(mfrow=c(2,2))                               # compare orginal and pruned tree

plot(tree.prod)
text(tree.prod,pretty=0)

plot(prune.prod)
text(prune.prod,pretty=0)

plot(mydata_3var$Por,mydata_3var$Brittle, col=grey(10:2/11)[cut.prod], pch=20, xlab="Porosity (%)",ylab="Brittleness (%)")
partition.tree(tree.prod, ordvars=c("Por","Brittle"), add=TRUE)

plot(mydata_3var$Por,mydata_3var$Brittle, col=grey(10:2/11)[cut.prod], pch=20, xlab="Porosity (%)",ylab="Brittleness (%)")
partition.tree(prune.prod, ordvars=c("Por","Brittle"), add=TRUE)

# Let's use our pruned tree to make predictions
yhat.prod = predict(prune.prod,newdata = test)  # predict with the tree

# We can plot the predictions vs. the test data set.
par(mfrow=c(1,1))                               # testing data vs. prediction
plot(yhat.prod,test$Prod,xlim = c(0,8000), ylim = c(0,8000))
abline(0,1)
# The binning is due to the fact that we only had specific number of terminal nodes.

# We can calculate the mean square error and the square root of the MSE for our prediction model.
MSE = mean((yhat.prod - test$Prod)^2)
SQRT_MSE = sqrt(MSE)

# We can compare the MSE to the total variance 
var.prod = var(test$Prod)

# On your own try working working with more and less complicated trees.  Also, for improved regression accuracy
# consider bagging, random forest and boosting.

# Hope this was helpful,

# Michael
#
# Michael J. Pyrcz, P.Eng., Ph.D.
# Associate Professor, the University of Texas at Austin (@GeostatsGuy)
# mpyrcz@austin.utexas.edu


