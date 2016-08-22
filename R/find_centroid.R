#' Finds the centroid of the simplex
#'
#' @param iter - variable containing all the information related with the current EPIC iteration
#' @return iter - updated information of the the current EPIC iteration

find_centroid <- function(iter){

  # Find the highest (worst) y value
  tmp <- sort(iter$simplex.y,index.return = TRUE,decreasing = TRUE)
  iter$simplex.w <- tmp$ix[1] # index of worst value
  iter$simplex.s <- tmp$ix[2] # index of second worst value
  iter$simplex.b <- tmp$ix[nrow(iter$simplex.x)] #index of best value

  # calculate the centroid
  # sum all vertices except the worst one and divide by number of included vertices
  p <- colSums (iter$simplex.x, na.rm = FALSE, dims = 1)
  p <- p - iter$simplex.x[iter$simplex.w,]
  p <- p /(nrow(iter$simplex.x)-1)
  iter$simplex.centroid <- p
  return(iter)
}
