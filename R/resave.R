#' Updates a .Rdata file
#'
#' Updates an .Rdata file by adding the variables specified in ellipsis to the rdata file.
#' No longer exported.
#' @param ... Objects in work space.  Will overwrite objects that already exist
#'  in file with the same name.
#' @param list Option to pass objects in work space as a character vector.
#' @param file An (potentially already existing) rdata file to save objects
#'  specified by ellipsis into.
#' @keywords save, resave, rdata
resave <- function(..., list = character(), file) {
  list <- union(list, as.character(substitute((...)))[-1L])
  # loads saved objects into new environment
  e <- new.env()
  if(file.exists(file)){
    load(file, e)
  }
  # copies objects from global environment
  # assumes they all exist
  for(n in list) assign(n, get(n, .GlobalEnv), e)
  save(list = ls(e), envir = e, file = file)
}
