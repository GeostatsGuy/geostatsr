# Linear Regression in R for Engineers and Geoscientists New to R
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# This will be used in my Introduction to Geostatistics Undergraduate Class 
# It is assumed that students have no previous R experience.  

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/PGE337")

# Read the data table from a comma delimited file
mydata = read.csv("Univariate_NonlinearPor_Perm_data.csv") # read csv file

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console
# The columns are variables with variable names at the top and the rows are samples

# Extract the variables of interest from the table into vectors.Well B porosity and permeability. 
# You are welcome to try linear regression with Well A also.
porB <- mydata$PorosityB                       # make a new vector from the PorosityB variable in our data table 
permB <- mydata$PermeabilityB                  # make a new vector from the PemreabilityB variable in our data table 
lnpermB <- log(mydata$PermeabilityB)           # make a new vector from the Natural Log of the PermeabilityB variable in our data table 

# It is a good idea to first explore our variables.
par(mfrow=c(3,2))                              # set up a 2x2 matrix of plots 

# We start with the scatter plot of Perm and Ln(perm) vs. Porosity for Well B. 
plot(permB,porB,main="Well B Permeability vs. Porosity",xlab=" Porosity B (%) ",ylab=" Permeability B ")
plot(lnpermB,porB,main="Well B Ln(permeability) vs. Porosity",xlab=" Porosity B (%) ",ylab=" ln(Permeability B) ")
# From these scatter plots it is clear that we should build a linear model of ln(Perm) vs. Porosity.

# Let's check the variable distributions porosity first and then ln(permeability)
hist(porB,main="Well B Porosity (%)",xlab="Porosity (%)",nclass = 15) # Hist builds a regular frequency histogram
# ecdf makes a cdf object and plot command plots it
plot(ecdf(porB),main="Well B Porosity (%)",xlab="Porosity (%)",ylab="Cumulative Probability")
# Now summary statistics, check the console
summary(porB)

# Now let's check the natural log of permeability distribution
hist(lnpermB,main="Well B Ln(permeability) (ln(mdD))",xlab="ln(permeability) (ln(mD))",nclass = 15) # hist builds a regular frequency histogram
# ecdf makes a cdf object and plot command plots it
plot(ecdf(lnpermB),main="Well B Ln(permeability) (ln(mD))",xlab="ln(permeability B) (ln(mD))",ylab="Cumulative Probability")
# Now summary statistics, check the console
summary(lnpermB)

# The data looks good. Distributions do mot have gaps, spikes nor outliers.

# To simplify our workflow let's make a new data table with just our 2 variables of interest
WellB <- data.frame(cbind(porB,lnpermB))       # command makes a new datafram combining our 2 vectors        
# To check let's look at it
head(WellB)
# and check the summary statistics. 
summary(WellB)

# Now we are ready to calculate a linear regression model to predict ln(permeability) from porosity
lnpermB.lm = lm(lnpermB ~ porB,data=WellB)     # our linear model predicts lmpermB from porB  

# To just see the resulting model use this command
print(lnpermB.lm)

# The linear model object ("lnpermB.lm") is very convenient. It includes all the information we need for hypothesis testing, and 
# confidence and prediction intervals.  In fact the hypothesis tests are completed automatically. 
lm.summary = summary(lnpermB.lm, correlation = TRUE)

# Let's check the statistical significance of the model coefficients.  
print(lm.summary$coefficients)
# Coefficients: You can observed the estimated model parameter, standard error and resulting t-statistic and the two-tailed
# probability of >|t-statistic| or the maximum level that one would reject at.  We reject the null hypothesis that 
# the model coefficients slope (porB) and intercept (intercept) are equal to 0.0.

print(lm.summary$fstatistic)
# F-statistic: is calculate from difference in sum of square error between the model and a constant slope model.
# We clearly reject the reduced, constant slope model in favour of the full model given fcrtical from f(0.95,DF1=1,DF2=103)=3.93

# Let's check the distribution of outliers for bias and outliers
hist(lm.summary$residuals,main="LM Residuals (ln(mD)",xlab="Residuals (ln(mD))",nclass = 15) # Hist builds a regular frequency histogram
summary(lm.summary$residuals)
# The residuals look good, the mean = 0.0 and there are no outliers.

# One could visualize all this output in the console with the summary command
summary(lnpermB.lm, correlation = TRUE) 

# Let's get the confidence intervals in the fitted model parameters
confint(lnpermB.lm,level = 0.95)
# We get the lower and upper 95% confidence intervals for both intercept and slope

# Let's demonstrate prediction intervals
prediction = predict(lnpermB.lm,interval="predict")
head(prediction)
# We now have a data table with the model fit and prediction intervals for all of our porosity data

# Furthermore lm has built in diagnostic plots that are very useful for checking our model
par(mfrow=c(2,2)) 
plot(lnpermB.lm)
# Residual vs Fitted to check linearity and homoscedasticity assumptions.  Check for large residuals / biased expectation.
# Normal Q-Q to check the normality in the residuals assumption.  Points should be on the 45 degree line.
# Scale-Location to check assumption of homoscedasticity.  There should not be a pattern.
# Residual vs. Leverage to identify poorly fit data that have large influence on the model (leverage).


