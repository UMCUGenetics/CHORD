#' Imitates the bash mkdir -p
#'
#' @param path Path
#' @param recursive Analogous to -p option
#' @param ... 
#'
#' @return NULL
#' @export
mkdir <- function(path, recursive=T, ...){
  if(!dir.exists(path)){
    dir.create(path, recursive=recursive, ...)
  }
}