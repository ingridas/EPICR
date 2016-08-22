#' Set of 9 DTLZ test functions for the case D -> D-9
#'
#' @param x is a vector of length n such that 0 < x_i < 1, i=1,...,n
#' @return a vector of length M = n - 9

dtlz1 <- function(x) {
  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- numeric(M)

  g <- 100 * (k + sum((x[M:n] - 0.5) ^ 2 - cos(20 * pi * (x[M:n] - 0.5))))

  for (id in 1:M){
    y[id] <- 1/2 * (1 + g)
  }
  for (id in 1:(M-1)){
    for (idd in 1:(M-id)){
      y[id] <- y[id] * x[idd]
    }
    if (id > 1){
      aux <- M - id + 1
      y[id] <- y[id] * (1 - x[aux])
    }
  }
  y[M] <- 1/2 * (1 - x[1]) * (1 + g)
  return(y)
}

dtlz2 <- function(x) {
  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- numeric(M)

  g <- sum((x[M:n] - 0.5)^2)

  for (id in 1:M){
    y[id] <- (1 + g)
  }
  for (id in 1:(M-1)){
    for (idd in 1:(M-id)){
      y[id] <- y[id]  *  cos(x[idd] * pi / 2)
    }
    if (id > 1){
      aux <- M - id + 1
      y[id] <- y[id] * sin(x[aux] * pi / 2)
    }
  }
  y[M] <- (1 + g) * sin(x[1] * pi / 2)
  return(y)
}

dtlz3 <- function(x){
  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- numeric(M)

  g <- 100 * (k + sum((x[M:n] - 0.5) ^ 2 - cos(20 * pi * (x[M:n] - 0.5))))

  for (id in 1:M){
    y[id] <- (1 + g)
  }
  for (id in 1:(M-1)){
    for (idd in 1:(M-id)){
      y[id] <- y[id]  *  cos(x[idd] * pi / 2)
    }
    if (id > 1){
      aux <- M - id + 1
      y[id] <- y[id] * sin(x[aux] * pi / 2)
    }
  }
  y[M] <- (1 + g) * sin(x[1] * pi / 2)
  return(y)
}

dtlz4 <- function(x){
  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- numeric(M)
  x <- x^100

  g <- sum((x[M:n] - 0.5)^2)

  for (id in 1:M){
    y[id] <- (1 + g)
  }
  for (id in 1:(M-1)){
    for (idd in 1:(M-id)){
      y[id] <- y[id]  *  cos(x[idd] * pi / 2)
    }
    if (id > 1){
      aux <- M - id + 1
      y[id] <- y[id] * sin(x[aux] * pi / 2)
    }
  }
  y[M] <- (1 + g) * sin(x[1] * pi / 2)
  return(y)
}

dtlz5 <- function(x) {

  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- numeric(M)
  theta <-  x[1:(M-1)]

  g <- sum((x[M:n] - 0.5)^2)
  t <- pi / (4 * (1 + g))
  theta[1] <- x[1] * pi / 2

  theta[2:(M-1)] <- t * (1 + 2 * g * x[2:(M-1)])

  for(id in 2:(M-1)){
    theta[id] <- t * (1 + 2 * g * x[id])
  }

  for (id in 1:M){
    y[id] <- (1 + g)
  }
  for (id in 1:(M-1)){
    for (idd in 1:(M-id)){
      y[id] <- y[id]  *  cos(theta[idd] * pi / 2)
    }
    if (id > 1){
      aux <- M - id + 1
      y[id] <- y[id] * sin(theta[aux] * pi / 2)
    }
  }
  y[M] <- (1 + g) * sin(x[1] * pi / 2)
  return(y)
}

dtlz6 <- function(x) {

  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- numeric(M)
  theta <-  x[1:(M-1)]

  g <- sum(x[M:n]^0.1)
  t <- pi / (4 * (1 + g))
  theta[1] <- x[1] * pi / 2

  theta[2:(M-1)] <- t * (1 + 2 * g * x[2:(M-1)])

  for(id in 2:(M-1)){
    theta[id] <- t * (1 + 2 * g * x[id])
  }

  for (id in 1:M){
    y[id] <- (1 + g)
  }
  for (id in 1:(M-1)){
    for (idd in 1:(M-id)){
      y[id] <- y[id]  *  cos(theta[idd] * pi / 2)
    }
    if (id > 1){
      aux <- M - id + 1
      y[id] <- y[id] * sin(theta[aux] * pi / 2)
    }
  }
  y[M] <- (1 + g) * sin(x[1] * pi / 2)
  return(y)
}

dtlz7 <- function(x){

  k <- 10 # recommended value
  n <- length(x)
  M <- n - k + 1
  y <- x[1:M]
  h <- 0

  g <- 1 + 9 * sum(x[M:n])/k

  for (id in 1:(M-1)){
    h <- h + y[id] / (1 + g) * (1 + sin(3 * pi  * y[id]))
  }
  y[M] <- (1 + g) * (M - h)
  return(y)
}

# this problem is constrained
dtlz8 <- function(x){

  n <- length(x)
  M <- n / 10  # recommended value
  y <- numeric(M)

  for (j in 1:M){
    id1 <- max(1,floor((j - 1) * (n / M)))
    id2 <- floor(j * n / M)
    y[j] <- 1 / (n / M) * sum(x[id1:id2])
  }
  return(y)
}

dtlz8_con <- function(x){

  n <- length(x)
  M <- floor(n / 10)  # recommended value
  y <- numeric(M)
  g <- numeric(M)

  for (j in 1:M){
    id1 <- max(1,floor((j - 1) * (n / M)))
    id2 <- floor(j * n / M)
    y[j] <- 1 / (n / M) * sum(x[id1:id2])
    }

  for (j in 1:(M-1)){
    g[j] <- y[M] + 4 * y[j] - 1
  }
  g[M] <- 2 * y[M] - 1

  min_val <- -999999
  for (i in 1:(M-1)){
    for (j in 1:(M-1)){
      if ((i != j) & ((y[i] + y[j]) < min_val)){
        min_val <- y[i] + y[j]
      }
    }
  }
  g[M] <- 2 * y[M] - 1 + min_val
  return(g)
}

# this problem is constrained
dtlz9 <- function(x){

  n <- length(x)
  M <- n / 10  # recommended value
  y <- numeric(M)

  for (j in 1:M){
    id1 <- max(1,floor((j - 1) * (n / M)))
    id2 <- floor(j * n / M)
    y[j] <- 1 / (n / M) * sum(x[id1:id2]^0.1)
  }
  return(y)
}

dtlz9_con <- function(x){

  n <- length(x)
  M <- floor(n / 10)  # recommended value
  y <- numeric(M)
  g <- numeric(M-1)

  for (j in 1:M){
    id1 <- max(1,floor((j - 1) * (n / M)))
    id2 <- floor(j * n / M)
    y[j] <- 1 / (n / M) * sum(x[id1:id2]^0.1)
  }

  for (j in 1:(M-1)){
    g[j] <- y[M]^2 + y[j]^2 - 1
  }
  return(g)
}
