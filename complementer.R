complement <- function(mat, loops = FALSE, new_val = mean){
  # Theory from:
  #   http://www.xatlantis.ch/index.php/education/zeus-framework/15-graph-theory
  mat[is.na(mat)] <- 0 # Take care of NAs
  bin <- all(mat%in%c(0, 1)) # Is the matrix binary?
  dir <- any(mat[lower.tri(mat)]!=mat[upper.tri(mat)]) # Is the matrix directed?
  # Determine the new value
  new_val <-   ifelse(bin, 1, # if the network is binary, set to 1
                      # Otherwise, if `new_val` is a function
                      ifelse('function'%in%is(new_val),
                             new_val(mat), # Apply the function,
                             new_val)) # Otherwise it must be a number
  if(as.numeric(new_val)|> is.na())stop('Problem with `new_val`,',
                                        'it does not yield a number!')
  matrix_complement(mat, new_val, bin, dir, loops)
}