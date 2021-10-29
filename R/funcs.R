#' @export
scrptflt <- function(fls, category){
  
  box::use(
    dplyr[...]
  )
  
  out <- fls %>%
    filter(category == !!category) %>%
    mutate(
      Index = as.character(1:nrow(.)),
      Link = paste0('[link](https://github.com/FloridaSEACAR/SEACAR/blob/master/', File, ')')
      ) %>% 
    select(Index, File, Link)
  
  return(out)
  
}

#' @export
scrptinp <- function(fl){
  
  box::use(
    dplyr[...],
    here[...]
  )

  tmp <- readLines(here(fl))
  tmp <- grep('read\\.csv|fread|st\\_read', tmp, value = T)
  tmp <- tmp[!grepl('\\#', tmp)] # remove commented lines
  
  out <- tmp %>% 
    gsub('^.*\"(.*)\".*$', '\\1', .) %>% 
    basename %>% 
    # sapply(tmp, function(x) x %>% basename %>% gsub('")$', '', .)) %>% 
    paste(collapse = ', ')

  if(nchar(out) == 0)
    out <- 'none'
  
  return(out)
  
}

#' @export
scrptout <- function(fl){
  
  box::use(
    dplyr[...],
    here[...]
  )

  tmp <- readLines(here(fl))
  tmp <- grep('^pdf|write\\.csv|^sink\\(', tmp, value = T)
  tmp <- tmp[!grepl('\\#', tmp)] # remove commented lines

  out <- sapply(tmp, function(x) x %>% strsplit('"') %>% .[[1]] %>% .[grepl('\\.pdf$|\\.csv$|\\.txt$', .)]) %>% 
    unlist
  
  if(length(out) == 0){
    return('none')
  }

  # hyperlink if files in output
  outfls <- list.files(here('docs/output'))
  fnd <- out %in% outfls
  out[fnd] <- paste0('[', out[fnd], '](https://FloridaSEACAR.github.io/SEACAR/output/', out[fnd], '){target="_blank"}')
  
  # concatenate
  out <- paste(out, collapse = ', ')
  
  return(out)
  
}

#' @export
scrptlib <- function(fl){
  
  box::use(
    dplyr[...],
    here[...]
  )
  
  tmp <- readLines(here(fl))
  tmp <- grep('^library\\(', tmp, value = T)
 
  out <- sapply(tmp, function(x) x %>% gsub('library\\(|\\)', '', .)) %>% 
    unlist %>% 
    unique %>% 
    sort %>% 
    paste(collapse = ', ')
  
  if(nchar(out) == 0)
    out <- 'none'
  
  return(out)
  
}