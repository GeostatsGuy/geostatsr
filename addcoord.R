# This function builds a spatial points dataframe with the locations for 
# estimation or simulation required by the gstat inverse distance, kriging 
# and simulation methods.  The grid is regular in 2D and parametrized with 
# Geo-DAS mehthod, as used in GSLIB.
# Michael Pyrcz, the University of Texas at Austin, @GeostatsGuy

addcoord <- function(nx,xmin,xsize,ny,ymin,ysize) {                       
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