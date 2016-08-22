#' Generates normalized weight vector of k dimensions and s divisions, so that each weight vector is near
#' the previous one. It uses the method of generating reflected k-ary Gray codes.
#'
#' @param k - dimension
#' @param s - divisions
#' @return The normalized weight vectors
#' @examples
#' generate_weights(10, 2)
#' generate_weights(5, 3)

generate_weights <- function(s,k){

  d <- k-1
  m <- s
  i <- 0:(s-1)
  wv <- matrix(0,s^k,k)
  wv[i+1,d+1] <- i # integer variable

  while(m < s^k){
    wv <- reverse(k,m,s,wv)
    d <- d-1
    m <- s*m
    i <- 0:(m-1)
    wv[i+1,d+1] <- floor(i/(m/s)) # integer variable
  }
  dwv <- wv/(s-1)
  sum2 <- rowSums(wv)
  ind <- which(sum2 == s-1)
  normwv <- dwv[ind,]

return(normwv)
}

## this function is used to ....
reverse <- function(k,n,s,wv){
  for (i in 0:(n-1)){
    for (h in 0:(s-2)){
      for (j in (0:(k-1))){
        if (h %% 2 == 0){
          wv[h*n+n+i+1,j+1] <- wv[n-i,j+1]
        } else {
          wv[h*n+n+i+1,j+1] <- wv[i+1,j+1]
        }
      }
    }
  }
  return(wv)
}
