array2GSLIB <- function(array,nx,ny,file_name,var_name) {       # Michael Pyrcz, April, 2018                      
  # makes a 2D array from realizations spatial point dataframe
  sink(file_name)
  cat(c(var_name,"\n"))
  cat(c("1","\n"))
  cat(c(var_name,"\n"))
  for(ix in 1:nx) {
    for(iy in 1:ny) {  
      cat(c(array[nx - ix + 1,iy],"\n")) 
    }
  }
  sink()
}