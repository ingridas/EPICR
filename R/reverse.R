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
