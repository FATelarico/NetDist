# Setup ====
## R packages ####
lapply(base::setdiff(c('countrycode', # Utilities
                       # Data manipulation
                       'data.table', 'dplyr', 'purrr', 'tidyr',
                       # Data visualisation
                       'ggplot2', 'ggthemes',
                       # Blockmodeling
                       'igraph'),
                     installed.packages()[, 1]), function(pkg){
                       install.packages(pkg)
                     })


## Functions ####
ggplot2_theme <- ggthemes::theme_tufte()+
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1),
    axis.line = ggplot2::element_line(linewidth = .3),
    panel.background = ggplot2::element_rect(fill = '#efefefef',
                                             colour = '#ffffffff'),
    strip.text = ggplot2::element_text(family = 'sans', face = 'bold')
  )

ggplot2_palette <- c(lightgrey = tryCatch(ggsci::pal_lancet()(9)[8],
                                          error = function(e)'#ADB6B6FF'),
                     cyan = tryCatch(gssci::pal_lancet()(9)[4],
                                     error = function(e)'#0099B4FF'))
## Data ####
# IMF DIP Data
#| Format: Observation per row
#| Dimensions:
  #| (a) Metadata:
  #|    - Country,
#|    - Counterpart Country,
#|    - Unit of Measure 
#| (b) Country:
    #|    - All countries, no aggregates
    #| (c) Indicator:
        #|    - Inward Direct investment       {Why? the narrative is about external
        #|                                      poles competing for African embedment}
        #|    - Liabilities (gross)            {Why? Preserves direction}
        #|    - Debt instruments               {Why? Data limitation for All Entities}
        #|    - All entities
        #| (d) Counterpart Country:
            #|    - All countries, no aggregates
            #| (e) Period:
                #|    - From 01/01/2011
                #|    - To 12/31/2023
                
                if(file.exists('./Code/Input/dat_raw.RDS')){
                  dat_raw <- readRDS('./Code/Input/dat_raw.RDS')
                } else {
                  dat_raw <- read.csv('./Code/Input/dataset.csv')
                  saveRDS(dat_raw, './Code/Input/dat_raw.RDS')
                  gc()
                }
                
                ## Limit to countries (excl. aggregates) ####
                if(file.exists('./Code/Input/dat_clean.RData')){
                  load('./Code/Input/dat_clean.RData')
                } else {
                  dat_clean <- split(dat_raw, dat_raw$TIME_PERIOD)|> lapply(function(x){
                    # cntry <- base::intersect(x$COUNTRY.ID, x$COUNTERPART_COUNTRY.ID) # OLD
                    cntry <- base::intersect(x$COUNTRY, x$COUNTERPART_COUNTRY)
                    pos <- countrycode::countrycode(origin = 'country.name.en',
                                                    destination = 'iso3c',
                                                    # origin = 'iso3c', # OLD
                                                    # destination = 'country.name.en', # OLD
                                                    sourcevar = cntry, nomatch = NA, warn = 0)
                    # Known countries and territories with unrecognised names
                    pos[grep('Hong Kong', cntry)] <- 'HKG'
                    pos[grep('Kosovo', cntry)] <- 'KHX'
                    pos[grep('Montserrat', cntry)] <- 'MST'
                    pos[grep('Anguilla', cntry)] <- 'AGL'
                    pos[grep('Netherlands Antilles', cntry)] <- 'ANT'
                    pos[grep('Sint Maarten', cntry)] <- 'STM'
                    pos[grep('Macao', cntry)] <- 'MCA'
                    pos[grep('Aruba', cntry)] <- 'ABW'
                    pos[grep('Curaçao', cntry)] <- 'CUR'
                    
                    cntry <- cntry[!is.na(pos)]
                    # cntry <- cntry[!is.na(cntry)] # OLD
                    pos <- pos[!is.na(pos)]
                    y <- x[(x$COUNTRY%in%cntry)&(x$COUNTERPART_COUNTRY%in%cntry), ]
                    # x[(x$COUNTRY.ID%in%cntry)&(x$COUNTERPART_COUNTRY.ID%in%cntry), ] # OLD
                    y$COUNTRY.ID <- pos[match(y$COUNTRY, cntry)]
                    y$COUNTERPART_COUNTRY.ID <- pos[match(y$COUNTERPART_COUNTRY, cntry)]
                    y
                  })
                  
                  ## Ensure each period only includes countries that are active in all years ####
                  # Countries by year
                  cntry <- lapply(dat_clean, function(x){
                    lapply(c('COUNTRY.ID', 'COUNTERPART_COUNTRY.ID'), function(i)x[[i]])|>
                      unlist()|> unique()
                  })
                  
                  
                  # OLD
                  ## Limit to sovereign states (excl. territories) ####
                  # svrgnty <- unlist(cntry)|> unique()|>
                  #   data.frame(what = _, is = c(
                  #     F, T, F, T, T, T, T, F, T, T, T, T, F, F, T, T, T, T, T, F, T, T, T, T, T,
                  #     T, T, T, T, T, T, T, T, T, T, T, T, F, T, T, F, F, T, F, T, T, T, T, T, T,
                  #     T, T, F, T, T, T, T, T, T, T, T, F, T, F, T, T, T, T, T, T, T, T, T, T, F,
                  #     T, T, T, T, T, T, T, F, T, T, T, T, T, F, T, T, T, F, T, T, T, F, T, F, T,
                  #     F, T, T, T, T, T, T, T, T, T, T, T, T, T, F, T, T, T, F, F, T, T, T, F, F,
                  #     T, T, T, T, T, F, T, T, T, T, T, T, T, T, T, F, T, F, T, T, T, T, F, T, T,
                  #     F, T, T, F, T, T, F, T, T, F, T, F, T, F, T, T, F, T, T, T, F, T, T, T, T,
                  #     T, T, F, T, F, T, T, T, T, T, F, F, T, T, T, T, T, F, T, T, T, T, T, T, T,
                  #     T, T, T, T, T, T, T, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T,
                  #     F, T, T, T, T, F, T, F, T, T, T, T, T, T
                  #   ))
                  # 
                  # dat_clean <- lapply(dat_clean, function(x){
                  #   pos <- svrgnty$what[svrgnty$is]
                  #   x[(x$COUNTRY.ID%in%pos)&(x$COUNTERPART_COUNTRY.ID%in%pos), ]
                  # })
                  # 
                  # cntry <- lapply(dat_clean, function(x){
                  #   lapply(c('COUNTRY.ID', 'COUNTERPART_COUNTRY.ID'), function(i)x[[i]])|>
                  #     unlist()|> unique()
                  # })
                  
                  # if all years have the same countries
                  if(lapply(cntry[-1], `%in%`, cntry[[1]])|> sapply(all)|> all()){
                    cntry <- unlist(cntry)|> unique()|> sort()
                  } else {
                    stop('Not all years have the same countries')
                  }
                  
                  save(dat_clean, cntry,
                       # svrgnty, # OLD
                       file = './Code/Input/dat_clean.RData')
                  gc()
                }
                
                # 1. Construct matrices ====
                if(file.exists('./Code/Input/Ms_lnBlss.RDS')){
                  Ms_lnBlss <- readRDS('./Code/Input/Ms_lnBlss.RDS')
                } else if(file.exists('./Code/Input/Ms_lnBlss.txt')){
                  Ms_lnBlss <- dget('./Code/Input/Ms_lnBlss.txt')
                } else {
                  ## 1.1 By year ####
                  Ms_year <- lapply(seq_along(dat_clean), function(t){
                    x <- dat_clean[[t]] # Data for this year
                    # Empty matrix
                    M <- matrix(NaN, length(cntry), length(cntry), dimnames = list(cntry, cntry))
                    i <- match(x$COUNTRY.ID, cntry) # Mark the rows
                    j <- match(x$COUNTERPART_COUNTRY.ID, cntry) # Mark the columns
                    M[cbind(i, j)] <- x$OBS_VALUE # Fill them
                    diag(M) <- 0 # No data for loops
                    
                    # Save information on the structural zeroes
                    attr(M, 'true0') <- which(as.vector(M)==0)
                    
                    # Mark the unreported values
                    attr(M, 'unreported') <- which(is.nan(as.vector(M)))
                    M[is.nan(M)] <- 0 # Overwrite them with zeroes
                    
                    # Mark the values reported missing
                    attr(M, 'missing') <- which(is.na(as.vector(M)))
                    M[is.na(M)] <- 0 # Overwrite them with zeroes
                    
                    M # Return the matrix with attributes
                  })|> `names<-`(names(dat_clean)); gc()
                  
                  # Turn stocks into flows
                  Ms_year <- lapply(seq_along(Ms_year)[-1], function(t){
                    out <- Ms_year[[t]]-Ms_year[[t-1]] # Compute differences
                    # Consider the attributes
                    for(ttrbt in c('true0', 'missing', 'unreported')){
                      # If same status in both, carry on the status
                      attr(out, ttrbt) <- base::intersect(attr(Ms_year[[t]], ttrbt),
                                                          attr(Ms_year[[t-1]], ttrbt))
                      # If only in $t-1$, mark it as 'prev'
                      attr(out, paste0(ttrbt, '_prev')) <-
                        base::setdiff(attr(Ms_year[[t-1]], ttrbt),
                                      attr(Ms_year[[t]], ttrbt))
                      # If only in $t$, mark it as 'pres'
                      attr(out, paste0(ttrbt, '_pres')) <-
                        base::setdiff(attr(Ms_year[[t]], ttrbt),
                                      attr(Ms_year[[t-1]], ttrbt))
                      
                    }
                    out
                  })|> `names<-`(names(Ms_year)[-1])
                  
                  saveRDS(Ms_year, './Code/Input/Ms_year.RDS')
                  
                  ## 1.2 By period ####
                  Ms_period <- lapply(seq(1, length(Ms_year), 4), function(t){
                    # Sum flows over the years in this period
                    M <- paste0('Ms_year[[', t:(t+3), ']]', collapse = ' + ')|>
                      parse(text = _)|> eval()
                    M <- M/4 # Take the average
                    
                    # Specially marked values
                    attrs <- for(ttrbt in c('true0', 'missing', 'unreported')){
                      x <- rep(FALSE, prod(dim(M))) # Vector of values
                      x <- sapply(0:3, function(tt){ # For each year in this period
                        x[attr(Ms_year[[t+tt]], ttrbt)] <- TRUE # Which units have `ttrbt`?
                        x
                      })
                      attr(M, ttrbt) <- which(rowSums(x)==4) # Units that always have `ttrbt`
                      x <- x[rowSums(x)<4, ]
                      for(i in 1:ncol(x)){ # Units that have `ttrbt`, but not always
                        attr(M, paste0(ttrbt, '_year', i)) <- which(x[, i])
                      }
                      # Remove first-differencing attributes
                      attr(M, paste0(ttrbt, '_prev')) <- NULL 
                      attr(M, paste0(ttrbt, '_pres')) <- NULL
                    }
                    M # Return matrix with attributess
                  })|> `names<-`(paste0('t', 1:3))
                  
                  saveRDS(Ms_period, './Code/Input/Ms_period.RDS')
                  
                  ## 1.3 Balassa normalisation ####
                  
                  Ms_Blss <- lapply(Ms_period, function(M){
                    r_tot <- rowSums(M)
                    c_tot <- colSums(M)
                    tot <- sum(M)
                    x <- sapply(seq_along(c_tot), function(j)sapply(seq_along(r_tot), function(i){
                      (M[i, j]/r_tot[i])*(tot/c_tot[j])
                    }))
                    x[is.nan(x)] <- 0 # 0/0 is NaN
                    stopifnot(all(!is.na(x)))
                    attributes(x) <- attributes(M)
                    attributes(x)$is_neg <- which(as.vector(M)<0)
                    x
                  })
                  
                  ## 1.4 Logged Balassa ####
                  Ms_lnBlss <- lapply(Ms_Blss, function(M){
                    pos <- which(M==0)
                    positivise <- abs(min(M))+1e-6
                    x <- M
                    x[-pos] <- x[-pos]+positivise
                    x <- log(x)
                    x[pos] <- 0 # Correct for ln(0)=-Inf
                    # Correct for ln(1)=0 and ln(0<x<1)<0
                    offset <- abs(trunc(min(x[-pos])))+1
                    x[-pos] <- x[-pos] + offset
                    x[pos] <- 0
                    attributes(x) <- attributes(M)
                    attributes(x)$offset <- offset
                    attributes(x)$neg_corr <- positivise
                    attributes(x)$is_neg_Blss <- which(as.vector(M)<0)
                    x
                  })
                  attr(Ms_lnBlss, 'metadata') <- c(list(
                    true0 = 'Stocks were reported as zero in years $t$ and $t-1$ year, for each year of the period.',
                    unreported = 'Stocks were unreported in years $t$ and $t-1$, for each year of the period.',
                    missing = 'Stocks were reported as `NA` the years $t$ and $t-1$, for each year of the period.',
                    true0_yearX = 'Flows computed on the reported FDI stocks were zero in the $X^\\mbox{th}$ year of this period.',
                    unreported_yearX = 'The FDI stocks for year $t$ or $t-1$ (where $t$ is the $X^\\mbox{th}$ year of this period) were unreported.',
                    missing_yearX = 'The FDI stocks for year $t$ or $t-1$ (where $t$ is the $X^\\mbox{th}$ year of this period) were reported as `NA`.',
                    is_neg = 'The average flows computed on the reported FDI stocks were negative.',
                    offset = 'All non-zero cells were offset by this quantity because by construction $\\ln(\\mbox{Balassa}(m_{i,j})):=0$ if $\\mbox{Balassa}(m_{i,j})=0$, but mathematically also $\\mbox{Balassa}(m_{i,j})=1\\Rightarrow \\ln(\\mbox{Balassa}(m_{i,j}))=0$.',
                    neg_corr = 'Correction applied to ensure all non-zero cells are positive',
                    is_neg_Blss = 'Cells with a negative Balassa index.'
                  ))
                  saveRDS(Ms_lnBlss, './Code/Input/Ms_lnBlss.RDS')
                  dput(Ms_lnBlss, './Code/Input/Ms_lnBlss.txt')
                }
                
                
                # 2. Descriptive statistics ====
                dir.create('./Code/Output', FALSE)
                if(file.exists('./Code/Output/desc_stats.dput')){
                  desc_stats <- dget('./Code/Output/desc_stats.dput')
                } else {
                  desc_stats <- list(
                    df = sapply(c(Ms_period, Ms_lnBlss), function(M){
                      pos <- base::union(which(as.vector(M)!=0), attributes(M)$true0)
                      c(Units = nrow(M),
                        Ties = sum(M!=0),
                        Density = round(sum(M!=0)/nrow(M)/(ncol(M)-1), 2),
                        Median = median(M[pos])|> round(2),
                        Mean = mean(M[pos])|> round(2),
                        SD = sd(M[pos])|> round(2),
                        Max = max(M[pos])|> round(2)) -> M
                      M <- as.character(M)|> `names<-`(names(M))
                    })|> as.data.frame()
                  )
                  
                  desc_stats$kbl <- desc_stats$df|> knitr::kable(
                    col.names = paste(seq(2012, 2023, 4),
                                      seq(2012+3, 2023, 4), sep = '-')|> rep(2),
                    format = 'latex',
                    align = 'c', caption = paste(
                      'Descriptive statistics for the global FDI network',
                      'across all thre periods, reported for both the nominal',
                      'bilateral flow matrix and its log-Balassa',
                      'transformation. The table shows the number of units,',
                      'observed ties, density, and the distribution of',
                      'non-structural zeroes (unreported flows) in each',
                      'period.'
                    ), label = 'DescrStats', digits = 0
                  )|> kableExtra::kable_styling(
                    latex_options = c('HOLD_position', 'scale_down', 'striped'),
                    row_label_position = 'r',
                    # bootstrap_options = c('striped', 'condensed')
                  )|> kableExtra::add_footnote(
                    c('Raw data is in million USD',
                      paste('All indicators computed considering non-zero',
                            'values and structural zeroes')),
                    notation = 'symbol'
                  )|> kableExtra::add_header_above(
                    c(' ' = 1, 'Raw data*' = 3,  'ln-Balassa' = 3),
                    bold = TRUE, italic = TRUE
                  )
                  
                  dput(desc_stats, './Code/Output/desc_stats.dput')
                  
                  zeroes_stats <- list(
                    df = sapply(Ms_year, function(M){
                      x <- c(
                        length(M|> attr('true0')),
                        length(M|> attr('unreported')),
                        length(M|> attr('missing'))
                      )/(
                        sum(M==0)-nrow(M)
                      )
                      x <- round(x*100, 2)
                      x[1] <- x[1]+(100-sum(x))
                      x
                    })|> as.data.frame()|>
                      `rownames<-`(c('true0', 'unreported', 'missing'))
                  )
                  
                  zeroes_stats$kbl <- zeroes_stats$df|> knitr::kable(
                    format = 'latex',
                    align = 'c', caption = paste(
                      'Detailed breakdown of the percentage of zeroes in the',
                      'data matrices by \'type\' of observation acros',
                      'all years.'
                    ), label = 'ZeroesStats', digits = 0
                  )|> kableExtra::kable_styling(
                    latex_options = c('HOLD_position', 'scale_down', 'striped'),
                    row_label_position = 'r',
                    # bootstrap_options = c('striped', 'condensed')
                  )|> kableExtra::add_footnote(
                    c('Columns may not sum to 100 due to rounding'),
                    notation = 'symbol'
                  )|> kableExtra::add_header_above(
                    c(1, rep(4, 3))|> `names<-`(c(' ',
                                                  paste('Period', 1:3))),
                    bold = TRUE, italic = TRUE
                  )
                  
                  dput(zeroes_stats, './Code/Output/zeroes_stats.dput')
                }
                
                if(file.exists('./Code/Output/box_plt.dput')){
                  box_plt <- dget('./Code/Output/box_plt.dput')
                } else {
                  box_plt <- lapply(1:4, function(i){
                    x <- list(Ms_period, Ms_period, Ms_Blss, Ms_lnBlss)[[i]]
                    lapply(seq_along(x), function(t){
                      M <- x[[t]]
                      pos <- base::union(which(as.vector(M)!=0), attributes(M)$true0)
                      tab <- data.frame(
                        x = as.vector(M[pos])*switch(i, 1e6, 1e6, 1, 1),
                        t = t, i = c('Raw data (USD)', 'logged',
                                     'Balassa', 'logged-Balassa')[i]
                      )
                      if(i==2){
                        tab <- tab[tab$x>0, ]
                        tab$x <- log(tab$x)
                      }
                      tab
                    })|> data.table::rbindlist()
                  })|> data.table::rbindlist()|> as.data.frame()
                  
                  box_plt$t <- factor(x = box_plt$t,
                                      levels = 1:3,
                                      labels = paste(seq(2012, 2023, 4), seq(2012+3, 2023, 4),
                                                     sep = '-'),
                                      ordered = TRUE)
                  box_plt$i <- gsub('logged', 'Logged', box_plt$i,
                                    ignore.case = FALSE)
                  box_plt$i <- factor(box_plt$i, levels = unique(box_plt$i),
                                      ordered = TRUE)
                  
                  
                  
                  dput(box_plt, './Code/Output/box_plt.dput')
                }
                dir.create('./Figures', FALSE)
                box_plt2 <- box_plt[box_plt$i!='logged', ]
                pdf('./Figures/01_BoxPlots.pdf')
                ggplot2::ggplot(box_plt2, ggplot2::aes(x = x, y = 1))+
                  ggplot2::geom_boxplot(
                    fill = '#ffffffff', color = ggsci::pal_uchicago()(9)[2],
                    linewidth = .25, size = .5*1.25
                  )+
                  ggplot2::geom_violin(alpha = .25,
                                       fill = ggsci::pal_uchicago()(9)[3],
                                       color = ggsci::pal_uchicago()(9)[2],
                                       linewidth = .25*1.25)+
                  ggplot2::facet_grid(t ~ i, scales = 'free_x')+
                  ggplot2::xlab('')+
                  ggplot2::scale_y_continuous('', breaks = NULL)+
                  ggthemes::theme_tufte(base_size = 11*1.25)+
                  ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90,
                                                        hjust = 1, vjust = .5)
                  )
                dev.off()
                
                # 3. Prepare data for `graph-tool` ====
                if(!exists('Ms_lnBlss')){
                  Ms_lnBlss <- if(file.exists('./Code/Input/Ms_lnBlss.RDS')){
                    readRDS('./Code/Input/Ms_lnBlss.RDS')
                  } else if(file.exists('./Code/Input/Ms_lnBlss.txt')){
                    dget('./Code/Input/Ms_lnBlss.txt')
                  }
                  message('Data loaded')
                }
                
                ## 3.1 Prepare edge lists ####
                ELs <- lapply(Ms_lnBlss, function(M)data.frame(
                  i = rep(rownames(M), times = ncol(M)),
                  j = rep(colnames(M),  each = nrow(M)),
                  weight = as.vector(M)
                ))
                ELs <- lapply(seq_along(ELs), function(t){
                  x <- ELs[[t]]
                  x$i <- paste0(x$i, '_', t)
                  x$j <- paste0(x$j, '_', t)
                  x
                })|> `names<-`(names(ELs))
                ELs <- c(ELs,
                         lapply(seq_along(ELs)[-1], function(t){
                           data.frame(
                             i = paste0(colnames(Ms_lnBlss[[t-1]]), '_', t-1),
                             j = paste0(colnames(Ms_lnBlss[[t]]), '_', t),
                             weight = 1
                           )
                         })|> `names<-`(paste0('D', seq_along(ELs)[-1])))
                
                
                EL <- data.table::rbindlist(ELs)
                
                ## 3.2 Turn into an igraph object####
                G <- igraph::graph_from_data_frame(EL, directed = TRUE)
                
                # All weights are non-negative
                stopifnot(all(igraph::E(G)$weight>=0))
                
                # Add the period as a vertex attribute
                igraph::V(G)$t <- as.integer(substr(igraph::V(G)$name, 5, 6))
                
                ## 3.3 Export to file ####
                igraph::write_graph(G, file = './Code/Input/lnBlss.gml', format = 'graphml')