#' Updates the simplex used in a local search of EPIC algorithm
#'
#' @param iter - list containing all the information about the current EPIC algorithm iteration
#' @param x - simplex vertice in the decision space
#' @param y - scalarized response in the objective space corresponding to the \code{x}
#' @return an updated list with the information about the current EPIC algorithm iteration

update_simplex <- function(iter,x,y){
  iter$simplex.x[iter$simplex.w,] <- x
  iter$simplex.y[iter$simplex.w] <- y
  iter$simplex.status <- "done"
  iter$simplex.centroid <- NULL
  iter$simplex.w <- NULL
  iter$simplex.s <- NULL
  iter$simplex.b <- NULL
  iter$simplex.ref.y <- NULL
  iter$simplex.ref.x <- NULL
  iter$simplex.exp.x <- NULL
  iter$simplex.exp.y <- NULL
  iter$simplex.contr.x <- NULL
  iter$simplex.contr.y <- NULL
  return(iter)
}
