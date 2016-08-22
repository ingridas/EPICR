#' Generates normalized weight vector of k dimensions and s divisions, so that each weight vector is near
#' the previous one. It uses the method of generating reflected k-ary Gray codes.
#'
#' @param y - objective values corresponding to \code{x}
#' @param x - decision values
#' @return yp - Pareto front in the objective space, xp - Pareto set in the decision space,
#' ip - logical indicator, TRUE if solution is nondominated, FALSE if solution is dominated
pareto_filter <- function(y,x) {
  d <- ncol(y)
  n <- nrow(y)
  is.optimal <- rep(TRUE, n)
  for(i in 1:(n-1)) {
    for (j in i:n) {
      if (i != j && (is.optimal[i] || is.optimal[j])) {
        yi <- y[i,]
        yj <- y[j,]
        if (all(yi <= yj) && any(yi < yj)) { ## i dominates j
          is.optimal[j] <- FALSE
        } else if (all(yj <= yi) && any(yj < yi)) { ## j dominates i
          is.optimal[i] <- FALSE
        }
      }
    }
  }
  yp <- y[is.optimal,,drop = FALSE]
  xp <- x[is.optimal,,drop = FALSE]
  ip <- is.optimal
  results <- list("P" = ip,"Yp" = yp,"Xp" = xp)
  return(results)
}
