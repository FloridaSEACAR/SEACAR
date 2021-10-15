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
  tmp <- grep('read\\.csv', tmp, value = T)
  tmp <- tmp[!grepl('\\#', tmp)] # remove commented lines
  
  out <- sapply(tmp, function(x) x %>% basename %>% gsub('")$', '', .)) %>% 
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
    unlist %>% 
    paste(collapse = ', ')
  
  if(nchar(out) == 0)
    out <- 'none'
  
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