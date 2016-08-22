#' Checks if a vector V is feasible in the decision space constrained by function con
#'
#' @param V - a decision vector
#' @param con - constraints function ( g <= 0)
#' @return True, if vector V is feasible; FALSE, otherwise

is_feasible_con <- function(con,V){
  if (all(con(V)>=0)){
    k <- TRUE
  }else{
    k<- FALSE
  }
  return(k)
}
