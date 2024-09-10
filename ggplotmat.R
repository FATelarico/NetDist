# Plot beautified matrix plots
ggplotMat <- function(
    # Arguments ####
    # Data
    M, clu = NULL,
    onemode = !('list'%in%is(clu)&&length(clu)==2),
    loops = NULL,
    # Labels and Axes
    legend.name = '', xlab = '', ylab = '',
    # Values
    print.vals = TRUE, round.vals = 2, mult.val = 1,
    # Appearance
    base.size = 11, cex = 1, cex.val = cex, themify = TRUE,
    col.val = NULL, limits = NULL, show_legend = TRUE
){
  # Check input matrix ####
  if(!'matrix'%in%is(M)){
    tryCatch({
      M <- as.matrix(M)
    }, error = function(e){
      stop('`M` cannot be coerced to a matrix')
    })
    message('Coerced `M` to a matrix')
  }
  # If a clustering was not provided ####
  # Assume singletons
  if(is.null(clu)){
    # If the matrix is square
    clu <- if(nrow(M)==ncol(M)){
      1:nrow(M) # Partition of singletons
    } else { # Otherwise
      # Partition of singletons for rows and cols
      list(1:nrow(M), 1:ncol(M))
    }
  }
  # Count untis by cluster ####
  if(onemode){ # If the network is one mode
    n <- split(clu, clu)|> lengths() # Simply count
    n[1] <- n[1] <- n[1]+.5 # Offset the first
  } else { # Otherwise
    # Count for each dimension
    n <- lapply(clu, function(x)split(x, x)|> lengths())
    # Offset for each dimension
    n <- lapply(n, function(x){
      y <- x
      y[1] <- y[1]+.5
      y
    })
  }
  
  # Store an ordered copy of the matrix
  if(onemode){ # For one mode networks
    mat <- M[clu|> order(), clu|> order()]
  } else {
    # For two-mode
    mat <- M[clu[[1]]|> order(), clu[[2]]|> order()]
  }
  # Create data tibble ####
  # Check column and row-names (cannot be NULL)
  if(is.null(rownames(mat))){
    rownames(mat) <- paste0('net', 1:nrow(mat))
  }
  if(is.null(colnames(mat))){
    colnames(mat) <- paste0('net', 1:ncol(mat))
  }
  # Turn into a tibble (works for one- and two-mode networks)
  M <- tibble::rownames_to_column(as.data.frame(mat),
                                  'from')|>
    tidyr::pivot_longer(-from, names_to = 'to',
                        values_to = 'val')
  
  # Take care of enforcing an order of the units
  M$from <- factor(M$from, ordered = TRUE,
                   levels = rownames(mat)|> rev())
                   # levels = rownames(mat))
  M$to <- factor(M$to, ordered = TRUE,
                 levels = colnames(mat)|> rev())
                 # levels = colnames(mat))
  
  
  
  # Automatic colour-scale selection####
  fill_scale <-
    if(any(mat[!is.na(mat)]<0)&any(mat[!is.na(mat)]>0)){
      # If there are both positive and negative values
      if(is.null(col.val)){
        ggplot2::scale_fill_gradient2(low = 'salmon3', high = 'darkgreen',
                                      mid = 'white', midpoint = 0,
                                      na.value = 'grey77', name = legend.name,
                                      limits = limits)
      } else {
        ggplot2::scale_fill_gradient2(low = col.val[1], high = col.val[2],
                                      mid = 'white', midpoint = 0,
                                      na.value = 'grey77', name = legend.name,
                                      limits = limits)
      }
      
    } else if(is.null(col.val)){
      if(all(mat[!is.na(mat)]<=0)){
        # If all values are negative
        ggplot2::scale_fill_gradient(low = 'darkred', high = 'salmon',
                                     na.value = 'grey77', name = legend.name,
                                     limits = limits)
      } else if(all(mat[!is.na(mat)]>=0)){
        # If all values are positive
        ggplot2::scale_fill_gradient(low = 'lightgreen', high = 'darkgreen',
                                     na.value = 'grey77', name = legend.name,
                                     limits = limits)
      }
    } else {
      ggplot2::scale_fill_gradient(low = col.val[1], high = col.val[2],
                                   na.value = 'grey77', name = legend.name,
                                   limits = limits)
  }
    
  # Draw cluster lines ####
  if(onemode){
    n1 <- sum(n) # Sum of units in the row-mode
                 # (which is the same as the col-mode)
    # Create lines between clusters
    clu_lines <- lapply(seq_along(n), function(i){
      x <- ifelse(i==1, n[i], sum(n[seq(1, i-1)]))
      c(
        # row-clusters
        paste0(
          'ggplot2::geom_segment(ggplot2::aes(y =', x, ', yend = ',
          x, ', x = .5, xend =', n1,'), linewidth = 2/7*', cex,
          ', linetype = 5)'
        )|> parse(text = _)|> eval(),
        # col-clusters
        paste0(
          'ggplot2::geom_segment(ggplot2::aes(x =', x, ', xend = ',
          x, ', y = .5, yend =', n1,'), linewidth = 2/7*', cex,
          ', linetype = 5)'
        )|> parse(text = _)|> eval()
      )
    })|> unlist()|> list()
  } else {
    n1 <- sum(n[[1]]) # Sum of units in the row-mode
    # Row clusters
    clu_lines <- lapply(seq_along(n[[2]]), function(i){
      x <- ifelse(i==1, n[[2]][i], sum(n[[2]][seq(1, i-1)]))
      paste0('ggplot2::geom_segment(ggplot2::aes(y =', x, ', yend = ',
             x, ', x = .5, xend = ', n1, '), linewidth = 2/7*', cex,
             ', linetype = 5)')|> parse(text = _)|> eval()
    })|> unlist()
    # Sum of units in the col-mode
    n2 <- sum(n[[2]])
    # Col clusters
    clu_lines <- c(
      clu_lines,
      lapply(seq_along(n[[1]]), function(i){
        x <- ifelse(i==1, n[[1]][i], sum(n[[1]][seq(1, i-1)]))
        paste0(
          'ggplot2::geom_segment(ggplot2::aes(x =', x, ', xend = ',
          x, ', y = .5, yend =', n2, '), linewidth = 2/7*', cex, ', linetype = 5)'
        )|> parse(text = _)|> eval()
      })|> unlist()
    )
  }
  # Diagonal ####
  # Set values that are zero on the diagonal to NA
  if(onemode&&!isTRUE(loops)){ # For one-mode  networks
    # Locate links on the diagonal
    pos <- which(as.vector(M$to)==as.vector(M$from))
    # If all the diagonal is zero and `loops` is NULL, null it out
    if(is.null(loops)&&{
      length(pos)>0&&all(M$val[pos]==0, na.rm = TRUE)
    }){
      M$val[pos][M$val[pos]==0] <- NA 
    } else if(loops==FALSE){
      M$val[pos][M$val[pos]==0] <- 0
    }
  }
  
  # Build the plot ####
  # Basic elements
  ggplot <- 
    ggplot2::ggplot(M, ggplot2::aes(x = from, y = to,
                                    fill = val))+
    ggplot2::geom_tile(linewidth = 2*cex, show.legend = show_legend)+
    fill_scale+
    ggplot2::xlab(xlab)+ggplot2::ylab(ylab)
  
  ## Thick external borders ####
  ggplot <- if(onemode){
    ggplot+
      ggplot2::geom_segment(
        ggplot2::aes(y = .5, yend = sum(n),
                     x = sum(n), xend = sum(n)),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )+
      ggplot2::geom_segment(
        ggplot2::aes(x = sum(n), xend = .5,
                     y = sum(n), yend = sum(n)),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )+
      ggplot2::geom_segment(
        ggplot2::aes(y = .5, yend = sum(n),
                     x = .5, xend = .5),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )+
      ggplot2::geom_segment(
        ggplot2::aes(x = sum(n), xend = .5,
                     y = .5, yend = .5),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )
  } else {
    ggplot+
      ggplot2::geom_segment(
        ggplot2::aes(y = .5, yend = sum(n[[2]]),
                     x = sum(n[[1]]), xend = sum(n[[1]])),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )+
      ggplot2::geom_segment(
        ggplot2::aes(x = 0.5, xend = sum(n[[1]]),
                     y  = sum(n[[2]]), yend = sum(n[[2]])),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )+
      ggplot2::geom_segment(
        ggplot2::aes(y = .5, yend = sum(n[[2]]),
                     x = .5, xend = .5),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )+
      ggplot2::geom_segment(
        ggplot2::aes(x = 0.5, xend = sum(n[[1]]),
                     y  = .5, yend = .5),
        linewidth = 1/2*cex, linetype = 1,
        lineend = 'round', linejoin = 'bevel'
      )
  }
  
  # Revert axis while keeping the direction of the diagonal
  ggplot <- ggplot+
    ggplot2::scale_x_discrete(limits=rev)+
    ggplot2::scale_y_discrete(limits=rev)
  
  ## Printing values ####
  if(print.vals){
    ggplot <- ggplot+
      ggplot2::geom_text(
        ggplot2::aes(label = round(val, round.vals)*mult.val),
        size = base.size/6*cex.val
      )
  }
  
  ## Theme ####
  ggplot <- if(themify){
    ggplot+
      ggthemes::theme_map(base_family = 'serif',
                          base_size = base.size*cex)
  } else {
    ggplot+
      ggplot2::theme_grey(base_family = 'sans',
                          base_size = base.size*cex)
  }

  ## Add theme and between cluster####
  ggplot+
    ggplot2::theme(
      legend.position = 'right',
      legend.text = ggplot2::element_text(size = base.size/2*cex),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, vjust = 1,
                                          size = 2/5*base.size*cex),
      axis.text.y = ggplot2::element_text(angle = 0, hjust = 1.2, vjust = .5,
                                          size = 2/5*base.size*cex),
      legend.title = ggplot2::element_text(size = 3/8*base.size*cex,
                                           hjust = 0.5, vjust = 0.5),
        legend.location = 'plot', legend.justification = 'center',
    )+
    clu_lines
}
