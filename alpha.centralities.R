alpha.centrality.sparse <- function(mat, exo = 1, alpha = 1, tol = 1e-07){
  n <- nrow(mat)
  M2 <- Matrix::sparseMatrix(dims = c(n, n), i = 1:n, j = 1:n, x = rep(1, n))
  exo <- cbind(rep(exo, length.out = n))
  M3 <- M2 - alpha * mat
  Matrix::solve(M3, tol = tol, exo)|> as.vector()
}

alpha.centrality.dense <- function(mat, exo = 1, alpha = 1, tol = 1e-07){
  n <- nrow(mat)
  exo <- rep(exo, length.out = n)
  exo <- matrix(exo, ncol = 1)
  id <- matrix(0, nrow = n, ncol = n)
  diag(id) <- 1
  as.vector(solve(id - alpha * mat, tol = tol) %*% exo)
}

alpha.centrality <- function(mat, exo = 1, alpha = 1, tol = 1e-07) {
  # Density
  ({sum(scales::rescale(mat, c(0, 1)))/ # Total normalised value of present ties
      ifelse(all(diag(mat)==0), # If there are no self-loops
                      nrow(mat)*(ncol(mat)-1), # Divide by $n(n-1)$
                      prod(dim(mat))) # Else, divide by $n^2$
    }<=.2)|> # Threshold for sparsity
    # Assign correct function
    ifelse(alpha.centrality.sparse, alpha.centrality.dense) -> f 
  
  # Execute computation
  f(mat, exo = exo, alpha = alpha, tol = tol)
}

alpha.seek <- function(mat, exo = 1, tol = 1e-07,
                       verbose = TRUE, returnNA = FALSE){
  # Determine spectral radius (largest eigenvalue of the matrix)
  spectr_radius <- 1/{eigen(mat)|> unlist()|> abs()|> max()}
  
  # # Find and Remove isolated nodes
  # pos <- {colSums(mat)==0}|{rowSums(mat)==0}
  # if(any(pos))mat <- mat[!pos, !pos]
  
  alpha <- -1 # Initialise alpha centrality
  
  # Search for the alpha parameter
  while(any(alpha<0)){ # As long as there are negative centralities
    # Reduce alpha
    spectr_radius <- spectr_radius-.01
    # Correct for error
    alpha <- tryCatch({
      alpha.centrality(mat, alpha = spectr_radius, exo = exo, tol = tol)
    }, error = function(e){alpha})
    # Stop if alpha was not found
    if(spectr_radius<0)if(returnNA){return(NA)}else{stop('Alpha not found')}
  }
  
  if(verbose)message('Alpha parameter =', spectr_radius, '\n')
  alpha
}
