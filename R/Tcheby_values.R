#' Converts a multiobjective problem into a scalarized single objective using weighted augmented Chebyshev approach
#' @param y - objective values of a multiobjective problem in a matrix form
#' @param gwv - a weigthing vector
#' @param absmin - minimal values of objectives
#' @param absmax - maximal values of objectives
#' @return k - objective value of a scalarized problem

Tcheby_values <- function(y,gwv,absmin,absmax){
  if (length(absmin) < 1){
    absmin <- apply(y, 2, min)
  }
  if (length(absmax) < 1){
    absmax <- apply(y, 2, max)
  }
  # check the dimensions of y and transpose if needed
  if(dim(y)[2] != length(absmin)){
    y <- t(y)
  }
    # normalize objective values to [0 1]
    yy <- sweep(y,2,absmin,"-")
    yy <- sweep(yy,2,absmax - absmin, "/")

    # multiply by weigths
    w.y <-sweep(yy,2,gwv,"*")

    # calculate the first term, i.e. find a maximum weighted value
    term1 <- apply(w.y,1,max)

    # calculate the second term of augmented weighted Tchebycheff equation
    rho <- 0.05 #default value
    term2 <- rowSums(w.y)*rho

    k <- term1+term2
    return(k)
  }
