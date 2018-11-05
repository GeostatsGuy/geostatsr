# 1D example to Demonstrate Overfitting 
# Michael Pyrcz, GeostatsGuy on Twitter, Assoc. Professor at the University of Texas at Austin

# load required libraries
library(car)

# Set the working directory
setwd("C:/PGE337")

# Load the data
mydata = read.csv("Bivariate_model_fit_check_data.csv") # read csv file
mydata

# Sampling Array in Depth
Depth <- (1:100)
depth_int_df <- data.frame(Depth)

# Extract fines and depth vectors
Fines <- mydata$Fines
Fines = (1-Fines)*0.3*100
Depth <- mydata$Depth
mydata$Fines = Fines 

# Matrix for plots
par(mfrow=c(3,2)) 

# Prediction grid
coords <- seq(from=0, to=100, by=0.1)
var(Fines)

# Model 1: Plot for depth and fines 
plot(mydata$Depth,mydata$Fines,xlab=" Depth (m) ",ylab=" Fines (%) ")
title("Model 1: 2nd Order Polynomial")

# Linear Regression Example 
fit <- lm(Fines ~ poly(Depth,2), data=mydata)

# Confidence interval for model

newx <- seq(min(mydata$Depth), max(mydata$Depth), length.out=100)

preds <- predict(fit, newdata = data.frame(Depth=newx),interval = 'confidence')

#plot(Fines ~ Depth, data = mydata, type = 'n')

polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)

lines(newx, preds[ ,3], lty = 'dashed',main="Distribution of Residuals", col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
lines(newx, preds[ ,1], col = 'black')
points(mydata$Depth,mydata$Fines)

#Calculate and plot residuals
resid = residuals(fit)
var(resid)
hist(resid, freq=FALSE, main="Distribution of Residuals",xlim=c(-20,20))


# Model 2: Plot for depth and fines 
plot(mydata$Depth,mydata$Fines,xlab=" Depth (m) ",ylab=" Porosity (%) ")
title("Model 2: 5th Order Polynomial")

# Linear Regression Example 
fit <- lm(Fines ~ poly(Depth,5), data=mydata)

# Confidence interval for model

newx <- seq(min(mydata$Depth), max(mydata$Depth), length.out=100)

preds <- predict(fit, newdata = data.frame(Depth=newx),interval = 'confidence')

#plot(Fines ~ Depth, data = mydata, type = 'n')

polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)

lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
lines(newx, preds[ ,1], col = 'black')
points(mydata$Depth,mydata$Fines)

#Calculate and plot residuals
resid = residuals(fit)
var(resid)
hist(resid, freq=FALSE, main="Distribution of Residuals",xlim=c(-20,20))


# Model 3: Plot for depth and fines 
plot(mydata$Depth,mydata$Fines,xlab=" Depth (m) ",ylab=" Porosity (%) ")
title("Model 3: 8th Order Polynomial")

# Linear Regression Example 
fit <- lm(Fines ~ poly(Depth,8), data=mydata)

# Confidence interval for model

newx <- seq(min(mydata$Depth), max(mydata$Depth), length.out=100)

preds <- predict(fit, newdata = data.frame(Depth=newx),interval = 'confidence')

#plot(Fines ~ Depth, data = mydata, type = 'n')

polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)

lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
lines(newx, preds[ ,1], col = 'black')
points(mydata$Depth,mydata$Fines)

#Calculate and plot residuals
resid = residuals(fit)
var(resid)
hist(resid, freq=FALSE, main="Distribution of Residuals",xlim=c(-20,20))




