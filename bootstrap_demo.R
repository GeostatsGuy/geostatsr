# Bootstrap in R for Engineers and Geoscientists New to R
# This will be used in my Intro to Geostatistics Undergraduate Class 
# It is assumed that Students have no previous R experience.  
# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first
library(boot)                                  # the package for bootstrap by Angelo Canty and Brian Ripley
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/my_directory")

# Let's make up some data. This is a n=5 set of porosity measurements in % 
porosity <- c(7, 9, 10, 11, 15)                # this command concatenates the list of numerical values into a vector
summary(porosity)                              # get basic summary statistics for our dataset 
porosity                                       # visualize our data, if a large dataset we could use head() to preview first m rows

# We could have loaded our data, let's write out our data and read it back in to demonstrate
write.csv(porosity,"Porosity.csv")             # writes out in comma delimited with column titles and 
loaded_porosity = read.csv("Porosity.csv")     # read csv file
loaded_porosity                                # view the loaded array          

# Lets make a PDF and CDF of our dataset
par(mfrow=c(2,2))                              # make an array of plots combine plots 
hist(porosity,main="Porosity (%)")
plot(ecdf(porosity),main="Porosity (%)",xlab="Porosity (%)",ylab="Cumulative Probability")
                                               # ecdf makes a cdf object and plot command plots it

# The bootstrap program requires a function that takes a dataset and array of indexes, and returns a summary statistic

# For simplicity lets use the average.  Here's our function:  
calc_average <- function(d,i=c(1:n)){ 
  d2 <- d[i]
  return(mean(d2))
}
# Let's take it apart and explain what it does

# First Parameter: d is the data vector or array
# Second Parameter: what is i=c(1:n)? i is the index vector expected to be the same size of the data and c(1:n) is a default
i <- c(1:5)
i                                              # so the default is just the indexes in order from 1 to n.  

#What does our function do with the index vector?  It builds a new data vector:
j = c(1,1,1,1,1)
test_array <- porosity[j]
test_array                                     # a vector with 5 values equal to the first element of the porosity vector
# This is resampling from our dataset with replacement, the basic building block of bootstrap

# To check, if we send this list of indexes to our function the average should just be equal to the first element

test_average = calc_average(porosity,j)
test_average

# It works! We have a function that tates our data array, a index vector and then calculates a statistic
# The boot function will take care of the repeated sampling (random index array) and keeps track of statistic realizations

# We are ready to run the bootstrap function (boot) with our dataset and function 
number_boot_samples = 1000                     # this is the number bootstrap samples, should be set to large enough
boot_por_avg <- boot(data=porosity,statistic=calc_average,R=number_boot_samples)

# We can access the bootstrap samples of average porosity directly as the component vector "t" 
boot_por_avg$t

# For visualization we can plot as histogram and CDF
hist(boot_por_avg$t,main="Bootstrap Porosity Average")
plot(ecdf(boot_por_avg$t),main="Bootstrap Porosity Average")

# Could also write out the bootrap realizations to a file, for later use or to load in to another platform 
write.csv(boot_por_avg$t,file="bootstrap_avg_por.csv")

# This was a simple example. There is so much more that can be done.  E.g. bootstrap confidence intervals, data weights, multisample etc.
# https://cran.r-project.org/web/packages/boot/boot.pdf

