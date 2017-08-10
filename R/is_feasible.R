#' Checks if a vector V is feasible in the decision space bounded by box constraints
#'
#' @param V - a decision vector
#' @param L - lower bound, a vector
#' @param U - upper bound, a vector
#' @return True, if vector V is feasible, FALSE, otherwise
#' @examples
#' is_feasible(c(3,5), c(2,1), c(5,4))
#' is_feasible(c(3,5), c(2,1), c(5,6))
#' @export

is_feasible <- function(V,L,U){
  if (all(V<=U) & all(V>=L)){
    k <- TRUE
  }else{
    k<- FALSE
  }
return(k)
}
