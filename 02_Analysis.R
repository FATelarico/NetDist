# Setup ====

## (i) Clean environment ####
rm(list = ls()); sink(tempfile()); gc(full = TRUE, verbose = FALSE); sink()

## (ii) Data: Raw partitions ####
res_woutDC <- paste0('./Code/partitions/',
                     'res_HSBM-woutDC-MCMC_projected_partitions.csv')|>
  read.csv()
res_withDC <- paste0('./Code/partitions/',
                     'res_HSBM-withDC-MCMC_projected_partitions.csv')|>
  read.csv()

## (iii) Data: Original networks ####
Ms <- readRDS('./Code/Input/Ms_period.RDS')
Ms_lnBlss <- readRDS('./Code/Input/Ms_lnBlss.RDS')

## (iv) Data: Period labels ####
periods <- rbind(seq(2012, 2023, 4), seq(2012+4-1, 2023, 4))|>
  apply(2, paste, collapse = ' - ')

## (v) Data: Match ISO alpha-3 codes and countries ####
ISOtoCntry <- data.frame(
  iso3c = rownames(Ms$t1),
  cntry = countrycode::countrycode(rownames(Ms$t1),
                                   'iso3c', 'country.name.en',
                                   nomatch = NA, warn = FALSE)
)
# Check for unmatchable ISO alpha-3 codes
stopifnot('There are no unmatchable ISO alpha-3 codes!'=
            any(is.na(ISOtoCntry$cntry)))

## (v) Data: WB indicators ####
if(file.exists('./Analysis/WB.RDS')){
  WB <- readRDS('./Analysis/WB.RDS')
} else {
  WB <- dir('./Analysis/', '_WB.tab', full.names = TRUE)
  WB <- WB|> lapply(read.delim)|> `names<-`(
    WB|> strsplit('//')|> sapply(`[`, 2)|> gsub('_WB.tab', '', x = _)
  )
  
  # Subset to years of interest
  WB <- lapply(WB, function(x){
    colnames(x)[-1:-4] <- paste0('Y', 1960:(1960+ncol(x)-5))
    # offset for metadata columns,
    # target year - first year
    # offset by one to find the right column
    x[, c(2, (4+(2012-1960)+1):(4+(2023-1960)+1))]
  })
  
  # Only relevant countries
  WB <- lapply(WB, function(x){
    x[x$Country.Code%in%ISOtoCntry$iso3c, ]
  })
  
  # Save
  saveRDS(WB, './Analysis/WB.RDS')
}

## (vi) Parameter: Which results to analyse ####
param <- askYesNo('Use results with or without DC?', default = 'TRUE',
                  prompts = c('without', 'with', 'both'))|>  as.character()|>
  switch('TRUE' = 1, 'FALSE' = 2, 'NA' = 3)

## (vii) Function: Custom matrix plots ####
source('./Code/Input/ggplotmat.R')

## (viii) Custom igraph plotting ####

igraph.plotting <- function(G, lyt, angle, vertex.label.cex = 1,
                            vertex.label.color = c(
                              rep('white', 3),
                              rep('black', igraph::gorder(G)-3))
){
  igraph:::plot.igraph(G,
                       mark.border = NA, # No border around the areas
                       # Rescale the edge width based on tie value
                       edge.width = igraph::E(G)$weight|>
                         scales::rescale(c(6, 16)),
                       # Colour edges depending on sign (hue) and value (alpha)
                       edge.color = igraph::E(G)$sign,
                       # Smaller arrows, straight lines
                       edge.arrow.size = .1, edge.curved = FALSE,
                       # layout = igraph::layout_with_sugiyama,
                       # Utilise custom layout
                       layout = lyt, rescale = TRUE,
                       # Vertex size directly from the number of countries
                       vertex.size = igraph::V(G)$size|> as.numeric(),
                       # Colour the vertices
                       vertex.color = igraph::V(G)$community,
                       # Use the frame to mark difference from the background
                       vertex.frame.color = 'white',
                       vertex.frame.width = as.numeric(igraph::V(G)$size)|>
                         max()|> sqrt(),
                       vertex.label.cex = 4*vertex.label.cex,
                       vertex.label.color = vertex.label.color,
                       edge.loop.angle = angle)
}

# 1. Clean-up the partitions ====
if(file.exists('./Code/Output/clu.RData')&&param==3){
  load(file = './Code/Output/clu.RData')
} else {
  if(param!=2){
    # Add a period index
    res_woutDC$t <- strsplit(res_woutDC$node, '_')|> sapply(`[`, 2)
    # Create a by-period matrix for each level
    clu_woutDC <- apply(res_woutDC[, -c(1:2, ncol(res_woutDC))], 2, function(x){
      y <- strsplit(res_withDC$node, '_')|> sapply(`[`, 1)
      split(x, res_woutDC$t)|> list2DF()|> `rownames<-`(y[!duplicated(y)])
    }, simplify = FALSE)# Compute clusters per period within levels
    k_woutDC <- sapply(X = clu_woutDC, FUN = apply,
                       MARGIN = 2, function(x) length(unique(x)))  
    # Identify trivial levels
    pos <- apply(k_woutDC, 2, min)>1
    # Exclude them
    k_woutDC <- k_woutDC[, pos]; clu_woutDC <- clu_woutDC[pos]  
    
    rm(pos)
  }
  if(param!=1){
    # Add a period index
    res_withDC$t <- strsplit(res_withDC$node, '_')|> sapply(`[`, 2)
    # Create a by-period matrix for each level
    clu_withDC <- apply(res_withDC[, -c(1:2, ncol(res_withDC))], 2, function(x){
      y <- strsplit(res_withDC$node, '_')|> sapply(`[`, 1)
      split(x, res_withDC$t)|> list2DF()|> `rownames<-`(y[!duplicated(y)])
    }, simplify = FALSE)
    # Compute clusters per period within levels
    k_withDC <- sapply(X = clu_withDC, FUN = apply,
                       MARGIN = 2, function(x) length(unique(x)))# Identify trivial levels
    pos <- apply(k_withDC, 2, min)>1
    # Exclude them
    k_withDC <- k_withDC[, pos]; clu_withDC <- clu_withDC[pos]  
    
    rm(pos)
  }
  
  # Export the results
  if(param==3){
    save(clu_withDC, clu_woutDC, file = './Code/Output/clu.RData')
    dput(list(withDC = k_withDC, woutDC = k_woutDC),
         file = './Code/Output/k.dput') 
    dput(list(withDC = clu_withDC, woutDC = clu_woutDC),
         file = './Code/Output/clu.dput')
  }
}

# 2. Preparatory analyses ====
# Limit to smaller layers
if(param!=1)clu_withDC <- clu_withDC$level_1
if(param!=2)clu_woutDC <- clu_woutDC$level_1

## 2.1 Reorder clusters by median total outgoing investment ####
if(file.exists('./Code/Output/clu_ord.RData')&&param==3){
  load(file = './Code/Output/clu_ord.RData')
} else {
  if(param!=1){
    ### 2.1.1 Results with DC ####
    ord_withDC <- lapply(seq_along(clu_withDC), function(t){
      x <- Ms[[t]]|> colSums() # Get total outgoing FDI by country in this period 
      # Split the totals by the clusters
      split(x, clu_withDC[[t]])|>
        # Compute the means of the GDPs by cluster, excluding NAs
        sapply(mean, na.rm = TRUE)|>
        # And reorder the clusters accordingly
        order(decreasing = TRUE)
    })
    
    clu_withDC <- sapply(seq_along(clu_withDC), function(t){
      clu <- clu_withDC[[t]]
      clu <- match(clu, sort(unique(clu)))
      # clu <- 
      match(clu, ord_withDC[[t]])
      # split(rownames(clu_withDC), clu)
    })|> `colnames<-`(names(clu_withDC))|>
      `rownames<-`(rownames(Ms$t1))
  }
  
  if(param!=2){
    ### 2.1.1 Results without DC ####
    ord_woutDC <- lapply(seq_along(clu_woutDC), function(t){
      x <- Ms[[t]]|> colSums() # Get total outgoing FDI by country in this period 
      # Split the totals by the clusters
      split(x, clu_woutDC[[t]])|>
        # Compute the means of the GDPs by cluster, excluding NAs
        sapply(mean, na.rm = TRUE)|>
        # And reorder the clusters accordingly
        order(decreasing = TRUE)
    })
    
    clu_woutDC <- sapply(seq_along(clu_woutDC), function(t){
      clu <- clu_woutDC[[t]]
      clu <- match(clu, sort(unique(clu)))
      match(clu, ord_woutDC[[t]])
    })|> `colnames<-`(names(clu_woutDC))|>
      `rownames<-`(rownames(Ms$t1))
  }
  
  ### 2.1.3 Export ####
  if(param==3){
    save(clu_woutDC, clu_withDC, file = './Code/Output/clu_ord.RData')
    dput(list(withDC = clu_withDC, woutDC = clu_woutDC),
         file = './Code/Output/clu_ord.dput')
    dput(list(withDC = ord_withDC, woutDC = ord_woutDC),
         file = './Code/Output/ord.dput')
  }
}

# 3. Country by cluster ====
if(file.exists('./Code/Output/cntr_all.RDS')&&param==3){
  cntr_all <- readRDS('./Code/Output/cntr_all.RDS')
} else {
  ## 3.1 Split ISO codes according to the clusters ####
  cntr_all <- list(
    # Using ISO Alpha-3 codes
    iso3c = lapply(switch(param, list(woutDC = clu_woutDC),
                          list(withDC = clu_withDC),
                          list(withDC = clu_withDC, woutDC = clu_woutDC)),
                   function(clu){
                     apply(clu, 2, function(Z){
                       split(rownames(clu), Z)
                     })|> `names<-`(paste0('t', 1:ncol(clu)))
                   })
  )
  
  ## 3.2 Split country names according to the clusters ####
  cntr_all$cntry <- lapply(cntr_all$iso3c, function(clu){
    lapply(clu, function(t)lapply(t, function(x){
      ISOtoCntry$cntry[match(x, ISOtoCntry$iso3c)]
    }))
  })
  
  ## 3.3 Export the results ####
  if(param==3)saveRDS(cntr_all, './Code/Output/cntr_all.RDS')
  sapply(seq_along(cntr_all$cntry), function(mthd){
    sapply(seq_along(cntr_all$cntry[[mthd]]), function(t){
      # File name
      md <- paste0('./Analysis/clu_', names(cntr_all$cntry)[mthd],
                   '_t', t, '.md')
      cat('', file = md, append = FALSE) # equivalent to `> file.md` in bash
      sapply(seq_along(cntr_all$cntry[[mthd]][[t]]), function(k){
        cat('Cluster ', k, '\n\n- ', sep = '', file = md, append = TRUE)
        cat(cntr_all$cntry[[mthd]][[t]][[k]], sep = '\n- ',
            file = md, append = TRUE)
        rep('\n', ifelse(k<=length(cntr_all$cntry[[mthd]][[t]]), 2, 1))|>
          cat(sep = '', file = md, append = TRUE)
      })
      TRUE
    })
  })#|> `colnames<-`(names((cntr_all$cntry)))
}
# 4. Compare clusterings using ARI ====
if(file.exists('Analysis/ARI.RDS')&&param==3){
  ARIs <- readRDS(file = 'Analysis/ARI.RDS')
} else {
  ARIs <- list()
  if(param==3){
    ARIs$byDC <- sapply(1:ncol(clu_withDC),
                        function(t)aricode::ARI(clu_withDC[, t],
                                                clu_woutDC[, t]))
  }
  ARIs$byT <- list()
  if(param!=1)ARIs$byT$withDC <- clu_withDC
  if(param!=2)ARIs$byT$woutDC <- clu_woutDC
    
  ARIs$byT <- lapply(ARIs$byT, function(x){
    sapply(1:ncol(x), function(i)sapply(1:ncol(x), function(j){
      aricode::ARI(x[, i], x[, j])
    }))|> `diag<-`(NA)
  })
  
  if(param==3)saveRDS(ARIs, file = 'Analysis/ARI.RDS')
}
# 5. Describe clusters using WB data ====

## 5.1  Clean WB data frames ####
WB <- if(file.exists('./Analysis/WB_clean.RDS')){
  readRDS('./Analysis/WB_clean.RDS')
} else {
  # Limit WB data to countries in the network
  WB <- lapply(WB, function(x)x[x$Country.Code%in%rownames(Ms$t1), ])
  # Which countries in the network have no data?
  pos <- lapply(WB, function(x)which(rownames(Ms$t1)%in%x$Country.Code))
  # Reorder WB data to match the order in the network
  WB <- lapply(WB, function(x)x[match(rownames(Ms$t1), x$Country.Code), ])
  # Export
  saveRDS(WB, './Analysis/WB_clean.RDS')
}

## 5.2 Compute summary statistics by cluster and period ####
if(file.exists('./Analysis/WB_tbls.RDS')){
  WB_tbls <- readRDS('./Analysis/WB_tbls.RDS')
} else {
  # Functions to be computed
  f <- list(Mean = function(x)mean(x, na.rm = TRUE),
            Median = function(x)median(x, na.rm = TRUE),
            SD = function(x)sd(x, na.rm = TRUE),
            IQR = function(x)IQR(x, na.rm = TRUE),
            min = function(x)min(x, na.rm = TRUE),
            max = function(x)max(x, na.rm = TRUE))
  
  # Results as summary tables
  WB_tbls <- lapply(
    switch(params, list(woutDC = clu_woutDC), list(withDC = clu_withDC),
           list(woutDC = clu_woutDC, withDC = clu_withDC)),
    function(x){
      lapply(WB, function(x)lapply(1:ncol(x), function(t){
        # Starting column (four columns for each previous period) + (one column)
        i <- (t-1)*4+1
        # End column: From (i) to (Three further columns) and (shift by one metadata)
        i <- i:(i+3) + 1
        # Compute country-means in this period
        dat <- x[, i]|> rowMeans()
        # Split means by cluster
        dat <- split(dat, x[, t])
        # Compute the functions by cluster
        sapply(f, function(foo)sapply(dat, function(y)foo(y)))
      })|> `names<-`(paste0('t', 1:ncol(x))))
    }
  )
  
  # Export
  saveRDS(WB_tbls, './Analysis/WB_tbls.RDS')
}

# 6. Create image matrices ====
#| Note: [M]_{i,j} Shows the YoY change in i's gross-debt position vis-a-vis j

if(file.exists('./Analysis/IMs.RData')){
  load(file = './Analysis/IMs.RData')
} else {
  ## 6.1 With DC ####
  if(params!=1){
    IM_withDC <- lapply(seq_along(Ms), function(t){
      lapply(list(raw = Ms, logBlss = Ms_lnBlss),
             function(M){
               lapply(list(mean = function(x)mean(x, na.rm = TRUE),
                           median = function(x)median(x, na.rm = TRUE),
                           sd =  function(x)sd(x, na.rm = TRUE)),
                      function(f){
                        blockmodeling::funByBlocks(M[[t]], clu_withDC[, t], fun = f)
                      })
             })
    })|> `names<-`(names(Ms))
  }
  
  ## 6.2 Without DC ####
  if(params!=2){
    IM_woutDC <- lapply(seq_along(Ms), function(t){
      lapply(list(raw = Ms, logBlss = Ms_lnBlss),
             function(M){
               lapply(list(mean = function(x)mean(x, na.rm = TRUE),
                           median = function(x)median(x, na.rm = TRUE),
                           sd =  function(x)sd(x, na.rm = TRUE)),
                      function(f){
                        blockmodeling::funByBlocks(M = M[[t]],
                                                   clu = clu_woutDC[, t],
                                                   FUN = f)
                      })
             })
    })|> `names<-`(names(Ms))
  }
  
  
  ## 6.3 Export ####
  if(params==3){
    save(IM_withDC, IM_woutDC, file = './Analysis/IMs.RData')
    dput(IM_withDC, file = './Analysis/IM_withDC.dput')
    dput(IM_woutDC, file = './Analysis/IM_woutDC.dput')
  }
}

if(param==3)param <- askYesNo('Analyse results with or without DC?',
                              default = 'TRUE',
                              prompts = c('without', 'with', 'cancel'))
x <- switch(param, IM_woutDC, IM_withDC, stop('Invalid `param`'))

if(FALSE)sapply(seq_along(x), function(t){
  if(t!=1)cat('\n\n')
  cat('Period', t, '\n\n')
  cat('Mean (log-BI)\n')
  knitr::kable(x[[t]]$logBlss$mean)|> cat(sep = '\n')
  cat('\n\nMedian (log-BI)\n')
  knitr::kable(x[[t]]$logBlss$median)|> cat(sep = '\n')
  cat('\n\nMean (USD YoY difference)\n')
  knitr::kable(x[[t]]$raw$mean)|> cat(sep = '\n')
  cat('\n\nSD (USD YoY difference)\n')
  knitr::kable(x[[t]]$raw$sd)|> cat(sep = '\n')
  cat('\n')
  if(askYesNo('Type `yes` to continue with the next period'))flush.console()
})

# 7. Sankey plots ====

clu <- switch(param, clu_woutDC, clu_withDC, stop('Invalid `param`'))

plt_snk <- lapply(2:ncol(clu), function(t){
  df <- data.frame(t1 = clu[, t-1], t2 = clu[, t])|> `rownames<-`(rownames(Ms$t1))
  
  df <- ggsankey::make_long(df, t1, t2)
  df$x <- periods[df$x|> gsub('t', '', x = _)|> as.integer()]
  
  plt <- ggplot2::ggplot(df,
                         ggplot2::aes(x = x, next_x = next_x, node = node,
                                      next_node = next_node, fill = factor(node),
                                      label = node)) +
    ggsankey::geom_sankey(flow.alpha = .7, node.color = "gray30") +
    ggsankey::geom_sankey_text(size = 8, color = "white") +
    # ggplot2::scale_fill_viridis_d(drop = FALSE, option = 'H') +
    ggplot2::scale_fill_manual(values = viridis::turbo(15)[c(2:4, 7:9, 12:15)])+
    ggsankey::theme_sankey(
      base_size = 32
    )+
    ggplot2::coord_flip()+
    ggplot2::xlab('')+
    ggplot2::theme(legend.position = 'none',
                   plot.title = ggplot2::element_text(hjust = .5),
                   axis.line.x = ggplot2::element_blank(), 
                   axis.line.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 20, face = 'bold'),
                   plot.caption = ggplot2::element_text(size = 16, 
                                                        colour = '#434343',
                                                        margin = ggplot2::margin(
                                                          0, 2, 0, 2
                                                        )),
                   plot.margin = grid::unit(c(0, 10, 20, 10), 'pt')
    )
  
  paste0('./Figures/Sankey_t', t-1, 't', t,'.pdf')|>
    pdf(width = 7*2, height = 7*2*3/4)
  print(plt)
  dev.off()
  plt
})|> `names<-`(paste0('t', (1:ncol(clu))[-ncol(clu)],
                      '_t', (1:ncol(clu))[-1]))

# 8. Mesoscopic graphs ====

plt_meso <- lapply(1:ncol(clu), function(t){
  ## 8.1 Create mesoscopic networks ####
  # t <- 1
  M <- blockmodeling::funByBlocks(Ms_lnBlss[[t]]-attr(Ms_lnBlss[[t]], 'offset'),
                                  clu[, t], fun = 'mean')
  G <- igraph::graph_from_adjacency_matrix(M, diag = TRUE,
                                           mode = 'directed', weighted = TRUE)
  
  
  ## 8.2 Aesthetics for the plot ####
  
  ### 8.2.1 Colour the edges based on the sign ####
  igraph::E(G)$sign <- viridis::turbo(n = 10)[c(2, 9)][ # Select two colours
    as.integer(sign(igraph::E(G)$weight)<0)+1 # Detect negative and positive ties
  ]|> gsub('FF$', '', x = _) # Remove the fourth bit (alpha)
  
  ### 8.2.2 Set edges' transparency to reflect edge value ####
  alpha_vals <- igraph::E(G)$weight # Copy the edge weights
  alpha_range <- range(abs(alpha_vals)) # Determine a common range
  if(.25>alpha_range[1]){
    alpha_range[2] <- alpha_range[2]+(.25-alpha_range[1])
    alpha_range[1] <- .25
  }
  # Separately rescale negative and positive values to the common range
  alpha_vals[alpha_vals<0] <-
    abs(alpha_vals[alpha_vals<0])|> scales::rescale(alpha_range)
  alpha_vals[alpha_vals>0] <-
    abs(alpha_vals[alpha_vals>0])|> scales::rescale(alpha_range)
  
  # Turn into alpha-channel values
  igraph::E(G)$sign <- 
    # Set alpha in the range $[55, 100]\subseteq\mathbb{N}$ on $[0, 255]$
    scales::rescale(alpha_vals, c(55, 100))|> round(0)|>
    as.hexmode()|> toupper()|> # Turn into hex and capitalise for consistency
    paste0(igraph::E(G)$sign, alpha = _) # Append to the edges' colour
  
  ### 8.2.3 Make all weights positive for the width of the edges ####
  igraph::E(G)$weight <- abs(igraph::E(G)$weight)
  
  ### 8.2.4 Set vertex size to the number of countries in the cluster ####
  igraph::V(G)$size <- table(clu[, t])|> unname()|>
    as.numeric()|> scales::rescale(to = c(15, 50))
  
  ### 8.2.5 Set colors as in the sankey ####
  igraph::V(G)$community <- viridis::turbo(15)[c(2:4, 7:9, 12:15)][
    as.integer(igraph::V(G)$name)
  ]
  # # -=-=-=-=-=-=-=-=-=
  # igraph::tkplot(G)
  # lyt2 <- igraph::tk_coords(3)
  # dput(lyt2, file = './Analysis/lyt_t3.dput')
  # # -=-=-=-=-=-=-=-=-=
  
  ## 8.3 Load coordinates ####
  lyt <- paste0('./Analysis/lyt_t', t, '.dput')|> dget(file = _)
  
  ## 8.4 Set other parameters ####
  angle <- switch(t, .5, .5, .1, stop())
  
  ## 8.5 Plot by period ####
  png(width = 768*1.5, height = 1024*2, bg = "transparent", type = "cairo-png", 
      filename = paste0('./Figures/orig_t', t, '.png'))
  igraph.plotting(G, lyt, angle, vertex.label.cex = 1.5)
  dev.off()
  
  list(G = G, lyt = lyt, edge.loop.angle = angle, M = M)
})

# ## 8.6 Plot all periods ####
# png(width = 768*1.5*3+124, height = 1024,
#     bg = "transparent", type = "cairo-png", 
#     filename = './Figures/meso_all.png')
# layout(matrix(nrow = 1, 1:3))
# for(t in plt_meso){print(
#   igraph.plotting(t$G, t$lyt, t$edge.loop.angle,
#                   vertex.label.cex = 1.25)
# )}
# dev.off()


## 8.7 Export data structures ####
saveRDS(plt_meso, file = paste0('./Analysis/plt_meso',
                                switch(param, 'woutDC', 'withDC'),
                                '.RDS'))

# 9. Matrix plots (aligned with the mesoscopic-graphs) ====

## 9.1 Load data ####
plt_meso <- readRDS(paste0('./Analysis/plt_meso',
                           switch(param, 'woutDC', 'withDC'), '.RDS'))
# Extract mesoscopic matrices
M_meso <- lapply(plt_meso, `[[`, 'M')|>
  `names<-`(paste0('t', seq_along(plt_meso)))

## 9.2 Plot matrices ####
plt_mesoM <- lapply(seq_along(M_meso), function(t){
  x <- ggplotMat(M = M_meso[[t]]-attr(Ms_lnBlss[[t]], 'offset'),
            onemode = TRUE, loops = TRUE,
            print.vals = TRUE, xlab = '', ylab = '', round.vals = 2,
            themify = TRUE, show_legend = FALSE, cex.val = 6)+
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 32),
      axis.text.y = ggplot2::element_text(size = 32)
    )
  pdf(width = 7*2, height = 7*1.9,
      file = paste0('./Figures/mesoM_t', t, '.pdf'))
  print(x)
  dev.off()
  x
})


saveRDS(plt_mesoM, file = paste0('./Analysis/plt_mesoM',
                                 switch(param, 'woutDC', 'withDC'),
                                 '.RDS'))
# 10. Export everything ====
save(list = ls(), file = paste0('./Analysis/Analysis',
                                switch(param, 'woutDC', 'withDC'),
                                '.RData'))
