#' This function performs simplex shrink operation. It replaces all but the best vertice by newly calculated
#  x_j = x_l + \delta(x_j-x_l) and f_j = f(x+j), for j=0,.,n , with j!=l.
#' @param iter - a list conatinig information related with the current algorithm iteration
#' @param delta - a constant used to shrink the simplex edges
#' @param con - constraints, an analytical function cheap to evaluate, in a form g_i(x)>=0, if available; otherwise it is equal to FALSE
#' @param L,U - row vectors of lower and upper bounds of the design space
#' @return iter -  an updated list
shrink_simplex <- function(iter,delta,con,L,U){
  large <- 999999999
  l <- iter$simplex.b
  xl <- iter$simplex.x[l,]
  i <- 0
  var <- ncol(iter$simplex.x)
  new.x <- matrix(NA,ncol = var)
  # the best vertice is left
  iter$simplex.x[1,] <- xl
  iter$simplex.y[1] <- iter$simplex.y[l]
  # indicator for assigning the follwing vertice if infeasible
  k <- 1
  for (j in 1:nrow(iter$simplex.x)){
    if (j!=l){
      #calculate new vertice
      new <- xl + delta*(iter$simplex.x[j,]-xl)
      # check if the obtained vertice is feasible
      if (is.function(con)){
        if (is_feasible(new,L,U) & is_feasible_con(con,new)){
          i <- i+1
          if(i==1){
            new.x[i,] <- new
          }else{
            new.x <- rbind(new.x,new)
          }
        }else{
          k <- k+1
          iter$simplex.x[k,] <- new
          iter$simplex.y[k] <- large
        }
      }else{
        if (is_feasible(new,L,U)){
          i <- i+1
          if(i==1){
            new.x[i,] <- new
          }else{
            new.x <- rbind(new.x,new)
          }
        }else{
          k <- k+1
          iter$simplex.x[k,] <- new
          iter$simplex.y[k] <- large
        }
      }
    }
  }
  if(i>0){
    iter$simplex.x[(k+1):(var+1),] <- new.x
    iter$x <- as.matrix(new.x)
  }else{
    stop('Shrinked simplex is infeasible') # for testing purpose
  }
  iter$ys <- NULL
  iter$simplex.status <- "done"
  iter$simplex.centroid <- NULL
  iter$simplex.w <- 0
  iter$simplex.s <- 0
  iter$simplex.b <- 0
  return(iter)
}
