# More updated than `BlockmodelingITN` (2024.12.24)
.create_subfigures <- function(files, width_prop, height,
                               captions, labels, hspace){
  sapply(seq_along(files), function(i)paste0(c(
    paste0('\\begin{subfigure}[t]{', width_prop[i], '\\textwidth}'),
    paste0('\t\\includegraphics[height=', height[i], 'in]{', files[i], '}'),
    paste0('\t\\caption{', captions[i], '}\\label{', labels[[i]], '}'),
    '\\end{subfigure}%',
    ifelse(i!=length(files), ifelse(hspace, '~', '\\thinspace'), '%')
  ), collapse = '\n'))|> paste(collapse = '\n')
}

include_subfigures <- function(files, valign = c(top = 't', bottom = 'b',
                                          special_page = 'p'),
                               float = rep('h', length(files)),
                               float_override = FALSE, hspace = FALSE,
                               width_prop = c(.5, .5), height = c(1.2, 1.2),
                               captions = c('Sub-fig 1', 'Sub-fig 2'),
                               caption = 'Figure with two sub-figures',
                               labels = c('fig:subfig1', 'fig:subfig2'),
                               label = 'fig:withsubfig', layout = NULL){
  align <- paste0(float, valign, ifelse(float_override, '!', ''),
                  collapse = '')|> gsub(pattern = ' ', replacement = '')
  if(is.null(layout)){
    layout <- matrix(seq_along(files), nrow = 1)
  }
  
  if(nrow(layout)==1){
    subfigures <- .create_subfigures(files, width_prop, height,
                                     captions, labels, hspace)
  } else {
    subfigures <- sapply(seq_len(nrow(layout)), function(i){ # For each row
      pos <- layout[i, ][!is.na(layout[i, ])]
      sapply(pos, function(j){
        .create_subfigures(files[j], width_prop[j], height[j],
                           captions[j], labels[j], hspace)
      })|> paste(collapse = '\n')
    })|> paste(collapse = '\n\\vspace{0.5cm}\n')
  }
  paste(c(paste0('\\begin{figure}', '[', align, ']'),
          '\\centering',
          subfigures,
          paste0('\\caption{', caption, '}'),
          paste0('\\label{', label, '}'),
          '\\end{figure}'), collapse = '\n')
}

# More updated than the one for `BlockmodelingITN`
ColourBox <- function(
    colour = color, color = 'Yellow', # Color of the background for box and title
    margin_lx = '5pt', # Left indent
    margin_rx = '5pt', # Right indent
    margin_tb = '4pt', # Top and bottom indent
    title = 'A box', # Title of the box
    label = 'box:a_box', # Label of the box
    text = 'This is a box', # Text inside the box
    full_page = FALSE, # If TRUE, the box will be the full page width
    size = 'normalsize' # Font size
){
  if(size!=''&&size!='normalsize'){
    text <- paste0('\\', size, '{', text, '}')
  }
  paste0('\\', ifelse(full_page, 'colorboxtextFullpage', 'colorboxtext'),
         '[', colour, ']{', margin_lx, '}{', margin_rx, '}{',
        margin_tb, '}{', title, '}{', label, '}{')|>
    paste(sep = '\n', paste0('\t', text), '}')|> cat()
}

# Absent in `BlockmodelingITN`
set_notation <- function(set_letter, elements_letter = tolower(set_letter),
                         one_subscript = TRUE, fun = 'mathcal',
                         delim = c('\\{', '\\}'), wrap = TRUE,
                          i = if(one_subscript){'i'}else{c('i', 'j')},
                          card = if(one_subscript){'n'}else{c('n', 'm')},
                          order_matters = ifelse(one_subscript, FALSE, TRUE)){
  stopifnot(length(i)>=ifelse(one_subscript, 1, 2))
  if(one_subscript)order_matters <- FALSE
  subscripts <- if(one_subscript){
    list(1:2)
  } else {
    if(order_matters){
      list(paste(c(1, 1), 1:2, sep = ',\\ '),
           paste(c(2, 2), 1:2, sep = ',\\ '))
    } else {
      list(paste0(c(1, 1), 1:2))
    }
  }
  
  paste0('\\', fun,'{', set_letter, '} = \\left', delim[1],
         rep(elements_letter, 2)|>
           paste0('_{', subscripts[[1]],'}')|> paste(collapse = ',\\ '),
         if(order_matters){
           paste0('\\ldots,\\ ', 
                  rep(elements_letter, 2)|>
                    paste0('_{', subscripts[[2]], '}')|> paste(collapse = ',\\ '))
         },
         '\\ldots,\\ ', paste0(elements_letter, '_{', paste(i, collapse = ',\\ '), '}'),
         if(order_matters){paste0(
           '\\ldots,\\ ',
           paste0(elements_letter, '_{', paste(rev(i), collapse = ',\\ '), '}')
         )},
         '\\ldots,\\ ', paste0(elements_letter, '_{', paste(card, collapse = ',\\ '), '}'),
         '\\right', delim[2]) -> x
  if(wrap){
    paste0('$', x, '$')
  } else {
    x
  }
}

matrix_latex <- function(mat, type = 'array', title = NULL, wrap = NULL,
                         label = NULL,
                         align = ifelse(type=='array', 'auto', '')){
  align <- if(type == 'array'){
    if(align == 'auto'){
      align <- rep('c', ncol(mat))|> paste0(collapse = '')
    }
    paste0('{', align, '}')
  } else {
    ''
  }
  title <- if(is.null(title)){
    c('', '')
  } else {
    c(paste0('\\begin{array}{c}\n\\textbf{', title,'}\\\\\n'), '\n\\end{array}')
  }
  wrap <- if(is.null(wrap)){
    c('', '')
  } else {
    if(is.null(label))label <- ''
    c(paste0('\\begin{', wrap, '}',
             ifelse(label=='', '', paste0('\\label{', label, '}')), '\n'),
      paste0('\n\\end{', wrap, '}'))
  }
  apply(mat, 1, function(x){
    paste0('x[', seq_along(x), ']')|> sapply(function(y)parse(text = y)|> eval())|>
      paste('\t', x = _)|> paste(collapse = ' & ')
  })|> paste(collapse = ' \\\\\n')|> 
    paste0(wrap[1], title[1], 
           '\\begin{', type, '}', align, '\n', x= _, '\n\\end{', type, '}',
           title[2], wrap[2])
}
