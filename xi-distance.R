# -----------------------------------------------------------------------------#
# # Run manually                                                               #
# save.exceptFun <- function(file, envir = .GlobalEnv){                        #
#   save(list = ls(envir = envir)[                                             #
#     ls(envir = envir)|> parse(text = _)|> lapply(eval)|>                     #
#       sapply(function(x)!is(x, 'function'))                                  #
#   ], file = file)                                                            #
# }                                                                            #
# if(!require('SPR', quietly = TRUE)){                                         #
#   devtools::install_github("TraME-Project/Shortest-Path-R")                  #
# }                                                                            #
# if(!require('Rcpp', quietly = TRUE))install.packages('Rcpp')                 #
# # Custom C++ code                                                            #
# sourceCpp('./Code/AMtoEL.cpp') # Adjacency matrix to edge list               #
# sourceCpp('Code/TransformPull.cpp') # Transformation for 'pull' ties         #
# sourceCpp('Code/TransformPush.cpp') # Transformation for 'push' ties         #
# sourceCpp('Code/complement.cpp') # Complement of a weighted directed network #
# sourceCpp('Code/table.cpp') # Table function in C++                          #
# # Fast implmentation of alpha centralities                                   #
# source('Code/alpha.centralities.R')                                          #
# # Graph-complementer                                                         #
# source('Code/complementer.R')                                                #
# -----------------------------------------------------------------------------#


# Shannon entropy ####
H <- function(x){
  x <- x[which(x>0)]
  -sum(x*log(x))
}

# Valued distance sequence of a network (VDS) ####
VDS <- function(N, n, tie_type, verbose = TRUE,
                internal = FALSE){
  # Set the right data-transformation
  M <- ifelse(tie_type=='push', TransformPush, TransformPull)
  # Apply the transformation (C + + function takes a vector)
  M <- M(as.vector(N))
  # Get an edge list using a custom C + + function
  EL <- adjMatrixToEdgeList(M)|> as.matrix()
  
  # Compute the Dijkstra distance
  longest_shortest_path <- 1
  if(n>100&verbose)message('Computing distance matrix...\n\t This may take some time')
  distmat <- sapply(1:n, function(i){
    if(verbose&&{n>50&&i%%(n/trunc(10))==1})message(
      paste0('\n\t', i-1, '/', n,' nodes processed')
    )
    d <- SPR::dijkstra(n_nodes = n, source_ind = i, arcs_matrix = EL)
    if(max(d$path_list)>longest_shortest_path)assign(
      'longest_shortest_path', value =  max(d$path_list), pos = -1
    )
    d$min_dist
  })
  
  #---------------------------------------------------------------------------#
  # Better than using a for loop                                              #
  # Unit: milliseconds                                                        #
  # expr      min       lq    mean   median       uq      max neval cld       #
  # sapply 36.91893 41.37102 49.0468 43.48353 49.32844 100.1460   100   a     #
  # for 42.27504 45.76691 50.4920 47.37920 50.91259 119.2321   100   a        #
  #---------------------------------------------------------------------------#
  
  # Unreachable nodes: set distance to sum of all original tie values'
  # absolute values or number of units (whichever is the largest)
  distmat[distmat==Inf] <- max(n, sum(abs(EL[, 3])))
  distmat <- round(distmat, 3)
  
  # Distinguish by destination of use
  if(internal){
    # Compute the probability distribution of the distances
    pdfm <- colMeans(distmat)^{1/3}
    # Compute the normalisation constant
    norm <- sum(pdfm[1:(n-1)]>0, 1)|> c(2)|> max()|> log()
    # Return the VDS as used for both steps of xi
    c(pdfm, max(c(0, H(pdfm)-H(distmat^{1/3})/n))/norm) 
  } else {
    # Find unique values of the distance
    delta <- as.vector(distmat)|> unique()|> sort()
    # Turn into probability distributions
    delta <- as.data.frame(distmat)|> lapply(function(x){
      # as.vector( # Get as a vector
      #   table( # Count the number of occurrences of each value
      #     # Turn into a factor to get all possible values listed
      #     factor(x, levels = delta, ordered = TRUE)
      #   )/n # Normalise to [0, 1]
      # )
      tablecpp(levels = delta, x = x)
    })
    # Compute mean probability distribution
    mu <- Reduce('+', delta) / n
    
    # Compute the Kullback-Leibler divergence
    KL <- sapply(delta, function(x){
      y <- x*log(x/mu)
      # As x approaches 0, the limit of (x)log(x/mu) is 0
      y[x==0] <- 0
      sum(y)
    })
    
    # Compute the information radius
    IRad <- sum(KL)/n
    # Compute VDS
    IRad/longest_shortest_path
  }
}

# xi-distance function ####
xi <-
  function(N, N1, w1 = .45, w2 = .45, w3 = .1,
           tie_type = c('push', 'pull'),
           loops = NULL, verbose = TRUE, 
           alpha_N = NULL, alpha_N1 = NULL,
           exo = 1, tol = 1e-07, return.components = FALSE){
  
    if(return.components){w1 <- w2 <- w3 <- 1/3}
    if(w1+w2+w3!=1)stop('The weights (`w1`, `w2`, and `w3`) must sum to one!')
    # d1 <- d2 <- d3 <- 0
    
    tie_type <- match.arg(tie_type, c('push', 'pull'))
    
    if(nrow(N)!=ncol(N)){
      stop('Network `N` is not one-mode')
    } else if(nrow(N1)!=ncol(N1)){
      stop('Network `N1` is not one-mode')
    }
    
    loops1 <- if(is.null(loops)){
        any(!is.na(diag(N1)))&&any(diag(N1)!=0)
      }
    loops <- if(is.null(loops)){
      any(!is.na(diag(N)))&&any(diag(N)!=0)
    }
    
    # Find number of nodes by network
    n <- nrow(N)
    n1 <- nrow(N1)
    
    # Find maximum number of nodes
    max_nNodes <- max(n, n1)
    
    CmprVDS <- matrix(0, ncol = max_nNodes)
    
    if(w1 + w2>0){
      VDS_N = VDS(N, tie_type = tie_type, n = n, internal = TRUE,
                  verbose = verbose)
      
      CmprVDS[1:(n-1)] = VDS_N[1:(n-1)]
      
      CmprVDS[length(CmprVDS)] <- VDS_N[n]
      
      VDS_N1 = VDS(N1, tie_type = tie_type, n = n1, internal = TRUE,
                   verbose = verbose)
      
      CmprVDS[1:(n1-1)] = CmprVDS[1:(n1-1)] + VDS_N1[1:(n1-1)]
      
      CmprVDS[length(CmprVDS)] <- CmprVDS[length(CmprVDS)] + VDS_N1[n1]
      
      CmprVDS <- CmprVDS/2
      
      d1 <- {(H(CmprVDS)- .5*(H(VDS_N[1:n]) + H(VDS_N1[1:n1])))/log(2)}|>
        c(0)|> max()|> sqrt()
      
      d2 <- abs(sqrt(VDS_N[n + 1])-sqrt(VDS_N1[n1 + 1]))
      
      #------------------------------------------------------------------------#
      # Better than implementing a pipe to a function or to the primitive      #
      # Unit: microseconds                                                     #
      # expr   min     lq      mean median     uq       max neval cld          #
      # - 2.912 4.2635   7.74898  6.686 10.528    25.501   100   a             #
      # Primitive 4.072 5.3650  11.80113  7.658 11.899   156.296   100   a     #
      # fun 4.683 6.6375 262.52363  7.877 13.719 25097.824   100   a           #
      #------------------------------------------------------------------------#
      
    } else {
      d1 <- d2 <- 0
    }
    
    if(w3>0){
      # Check if N is diagonal
      N_IsDiag <- N[-(({0:(n-1)}*n)+{1:n})]
      N_IsDiag <- all(is.na(N_IsDiag))||all(N_IsDiag==0)
      # Check in N1 is diagonal
      N1_IsDiag <- N1[-(({0:(n1-1)}*n1)+{1:n1})]
      N1_IsDiag <- all(is.na(N1_IsDiag))||all(N1_IsDiag==0)
      #----------------------------------------------------------------------#
      # Faster than using an ad-hoc function                                 #
      # Unit: microseconds                                                   #
      # expr    min      lq      mean  median     uq       max neval cld     #
      # fun 24.194 38.9315 611.42112 47.8485 65.158 52877.470   100   a      #
      # code 23.104 37.0435  60.13884 45.9620 64.023   264.938   100   a     #
      #----------------------------------------------------------------------#
      
      
      w3_ToBeHalved <- skip_this_part <- 0
      # If noth both are diagonal
      if(N_IsDiag+N1_IsDiag<=2){
        # Compute alpha centralities
        tryCatch(expr = {
          alpha_N <- if(N_IsDiag){
            rep(0, n)
          } else if(!is.null(alpha_N)){
            alpha.centrality(N, alpha = alpha_N, exo = exo, tol = tol)
          } else {
            tryCatch(expr = {
              alpha.seek(N, verbose = verbose, exo = exo, tol = tol)
            }, error = function(e){
              alpha.centrality(N, alpha = 1/n, exo = exo, tol = tol)
            })
          }
        }, error = function(e){
          skip_this_part <- TRUE
        })
        
        tryCatch(expr = {
          alpha_N1 <- if(N1_IsDiag){
            rep(0, n1)
          } else if(!is.null(alpha_N1)){
            alpha.centrality(N1, alpha = alpha_N1, exo = exo, tol = tol)
          } else {
            tryCatch(expr = {
              alpha.seek(N1, verbose = verbose, exo = exo, tol = tol)
            }, error = function(e){
              alpha.centrality(N1, alpha = 1/n1, exo = exo, tol = tol)
            })
          }
        }, error = function(e){1
          skip_this_part <- TRUE
        })
        # If both alpha centralities were computed
        if(!skip_this_part){
          # Those the alpha centrality of at least one network exist?
          skip_this_part <- all(is.na(alpha_N))||all(is.na(alpha_N1))
        }
        if(!skip_this_part){
          # If one of the networks has no alpha centrality
          # Assign zeroes as if it was diagonal
          if(all(is.na(alpha_N1))){
            alpha_N1 <- rep(0, n1)
          } else if(all(is.na(alpha_N))){
            alpha_N <- rep(0, n)
          } # Then, continue
          
          # Deal with negative centralities
          alpha_N <- abs(alpha_N)
          alpha_N1 <- abs(alpha_N1)
          # Create empty matrices
          alpha_N_mat <- alpha_N1_mat <- matrix(0, ncol = max_nNodes)
          
          # Fill them
          alpha_N_mat[(max_nNodes-n + 1):max_nNodes] <- alpha_N
          alpha_N1_mat[(max_nNodes-n1 + 1):max_nNodes] <- alpha_N1
          
          # Compute the first part of the third component on the distance
          # based on the centrality of nodes computed on present ties
          # (1) Difference between:
          #     - the entropy of the mean-centrality matrix; and
          #     - the mean of the entropies of the centrality vectors.
          # (2) Divided by the natural logarithm of two
          # (3) Square root
          # (4) Divided by two
          d3 <- sqrt({
            H({alpha_N_mat + alpha_N1_mat}/2)-(H(alpha_N) + H(alpha_N1))/2
          }/log(2))/2
          
          #------------------------------------------------------------------------#
          # Better than splitting the operations for readability                   # 
          # Unit: microseconds                                                     #
          # expr    min     lq     mean median     uq     max neval cld            #
          # one  5.924 7.2425 10.16504  7.467 7.7005  83.489   100   a             #
          # more  5.165 6.7685 15.00367  6.977 7.2615 383.371   100   a            #
          #------------------------------------------------------------------------#
        } else {
          # If neither has alpha centrality or at least one was not computed,
          # skip this part and deal with it
          w3_ToBeHalved <- 1
          # Initialise d3
          d3 <- 0
        }
      } else {
        # If both are diagonal
        # (1) Warn the user
        if(verbose)warning(
          paste0('Both networks\' adjacency matrix is diagonal, ',
                 'alpha centrality cannot be computed!\n',
                 '\t Relying only on the complements!')
        )
        # If neither has alpha centrality, skip this part
        w3_ToBeHalved <- w3_ToBeHalved+1
        # Initialise d3
        d3 <- 0
      }
      # Continue with the complements
      skip_this_part <- FALSE
      # Compute complements
      N <- complement(N, loops = loops, new_val = mean)
      N1 <- complement(N1, loops = loops, new_val = mean)
      
      # Complements cannot be diagonal!
      # # Check if the complements are diagonal
      # N_IsDiag <- N[-(({0:(n-1)}*n)+{1:n})]
      # N_IsDiag <- all(is.na(N_IsDiag))||all(N_IsDiag==0)
      # N1_IsDiag <- N1[-(({0:(n1-1)}*n1)+{1:n1})]
      # N1_IsDiag <- all(is.na(N1_IsDiag))||all(N1_IsDiag==0)
      
      # Alpha centralities of the complements
      tryCatch(expr = {
        alpha_N <- tryCatch(expr = {
          alpha.seek(N, verbose = verbose, exo = exo, tol = tol)
        }, error = function(e){
          alpha.centrality(N, alpha = 1/n, exo = exo, tol = tol)
        })
      }, error = function(e){
        skip_this_part <- TRUE
      })
      
      tryCatch(expr = {
        alpha_N1 <- tryCatch(expr = {
          alpha.seek(N1, verbose = verbose, exo = exo, tol = tol)
        }, error = function(e){
          alpha.centrality(N1, alpha = 1/n1, exo = exo, tol = tol)
        })
      }, error = function(e){
        skip_this_part <- TRUE
      })
      
      # If both alpha centralities were computed
      if(!skip_this_part){
        # Those the alpha centrality of at least one network exist?
        skip_this_part <- all(is.na(alpha_N))||all(is.na(alpha_N1))
      }
      if(!skip_this_part){
        # If one the complements has no alpha centrality
        # Assign zeroes as if it was diagonal
        if(all(is.na(alpha_N1))){
          alpha_N1 <- rep(0, n1)
        } else if(all(is.na(alpha_N))){
          alpha_N <- rep(0, n)
        } # Then, continue
        
        # Deal with negative centralities
        alpha_N <- abs(alpha_N)
        alpha_N1 <- abs(alpha_N1)
        
        # Empty matrices
        alpha_N_mat <- alpha_N1_mat <- matrix(0, ncol = max_nNodes)
        
        # Fill them
        alpha_N_mat[(max_nNodes-n + 1):max_nNodes] <- alpha_N
        alpha_N1_mat[(max_nNodes-n1 + 1):max_nNodes] <- alpha_N1
        
        # Compute the first part of the third component on the distance
        # based on the centrality of nodes computed on complementary ties
        # (1) Difference between:
        #     - the entropy of the mean complement-centrality matrix; and
        #     - the mean of the entropies of the complement-centrality vectors.
        # (2) Divided by the natural logarithm of two
        # (3) Square root
        # (4) Divided by two
        # (5) Sum to the first part of the third component
        d3 <- sqrt({
          H((alpha_N_mat + alpha_N1_mat)/2)-(H(alpha_N) + H(alpha_N1))/2
        }/log(2))/max(alpha_N_mat, alpha_N1_mat) + d3
      } else {
        # If neither has alpha centrality, skip this part
        w3_ToBeHalved <- w3_ToBeHalved+1
      }
      
      # If there were issues with the alpha centralities, w3 has to be reduced!
      if(w3_ToBeHalved!=0){
        # If w3_ToBeHalved == 2, all of w3 is redistributed to w1 and w2
        # If w3_ToBeHalved == 1, half of w3 is redistributed to w1 and w2
        div <- ifelse(w3_ToBeHalved==1, 2, 1)
        # Proportion for the redistribution between w1 and w2
        prop <- w1/(w1+w2) # Proportion of w11 on w1+w2=1-w3
        # Increase w_1
        w1 <- w1+prop*(w3/div)
        # Increase w_2
        w2 <- w2+(1-prop)*(w3/div)
        # Reduce w_3
        w3 <- w3 - w3/div
      }
      
    } else {
      d3 <- 0
    }
    
    if(return.components){
      c(d1, d2, d3)
    } else {
      # Weight the components
      w1* d1 + w2*d2 + w3*d3
    }
  }

xi_norm <- function(N, N1, w1 = .45, w2 = .45, w3 = .1,
                    tie_type = c('push', 'pull'),
                    loops = NULL, verbose = TRUE, 
                    alpha_N = NULL, alpha_N1 = NULL,
                    exo = 1, tol = 1e-07, normalise = FALSE){
  
  d <- xi(N, N1, w1, w2, w3, tie_type, loops, verbose, 
          alpha_N, alpha_N1, exo, tol)
  
  max_xi <- max(
    xi(N, complement(N), .45, .45, .1, 'pull', verbose = FALSE),
    xi(N1, complement(N1), .45, .45, .1, 'pull', verbose = FALSE)
  )
  d/max_xi
}
