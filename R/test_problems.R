#' Set of multiobjective benchamrk test including ZDT test functions, OKA2, Viennet,9 of DTLZ problems and TF4 (a constrained problem)
#' @param x - is a decision vector of lenght \code{n}
#' @return a vector of function values
#' @examples
#' oka2(c(1,2,3))
#' zdt1(c(2,2,2))

zdt1 <- function (x) {
  f <- numeric(2)
  n <- length(x)
  f[1] <- x[1]
  g <- 1 + 9 * mean(x[2:n])
  f[2] <- g * (1 - sqrt(f[1] / g))
  return(f)
}

zdt2 <- function (x) {
  f <- numeric(2)
  n <- length(x)
  f[1] <- x[1]
  g <- 1 + 9 * mean(x[2:n])
  f[2] <- g * (1 - (f[1] / g)^2)
  return(f)
}

zdt3 <- function (x) {
  f <- numeric(2)
  n <- length(x)
  f[1] <- x[1]
  g <- 1 + 9 * mean(x[2:n])
  f[2] <- g * (1 - sqrt(f[1]/g) - f[1]/g * sin(10 * pi * f[1]))
  return(f)
}

zdt4 <- function (x) {
  f <- numeric(2)
  n <- length(x)
  f[1] <- x[1]
  g <- 1 + 10*(n-1) + sum(x[2:n]^2 - 10 * cos(4 * pi * x[2:n]))
  f[2] <- g * (1 - sqrt(f[1] / g))
  return(f)
}

zdt6 <- function (x) {
  f <- numeric(2)
  n <- length(x)
  f[1] <- 1 - exp(-4 * x[1]) * (sin (6 * pi * x[1]))^6
  g <- 1 + 9 * mean(x[2:n])^0.25
  f[2] <- g * (1 - (f[1]/g)^2)
  return(f)
}

oka2 <- function (x) {
  f <- numeric(2)
  f[1] <- x[1]
  f[2] <- 1 - (1./(4.*pi^2))*(x[1] + pi)^2 + abs(x[2] - 5*cos(x[1]))^(1./3) + abs(x[3] - 5.*sin(x[1]))^(1./3)
  return(f)
}

viennet <- function (x) {
  f <- numeric(3)
  f[1] <- 0.5 * (x[1]^2 + x[2]^2) + sin(x[1]^2 + x[2]^2)
  f[2] <- 0.125 * (3*x[1] - 2*x[2] + 4)^2 + (1.0/27.0) * (x[1] - x[2] + 1)^2 + 15
  f[3] <- 1.0/(x[1]^2 + x[2]^2 + 1) - 1.1*exp(-(x[1]^2 + x[2]^2))
  return(f)
}

kursawe <- function (x) {
  n <- length(x)
  f <- numeric(2)
  a <- 0.8
  b <- 3
  f[1] <- sum(-10*exp(-0.2*sqrt(x[1:n-1]^2 + x[2:n]^2)))
  f[2] <- sum(abs(x[1:n])^a + 5*sin(x[1:n]^b))
  return(f)
}

# this function is constrained
TF4 <- function (x) {
  y <- numeric(2)
  y[1] <- x[1]^2 -x[2]
  y[2] <- -0.5*x[1]-x[2]-1
  return(y)
}

TF4_con <- function (x) {
  g <- numeric(3)
  g[1] <- 6.5 - x[1]/6 - x[2]
  g[2] <- 7.5 - 0.5*x[1] - x[2]
  g[3] <- 20 - 5*x[1] - x[2]
  return(g)
}

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
    id1 <- floor((j - 1) * (n / M)) + 1
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
    id1 <- floor((j - 1) * (n / M)) + 1
    id2 <- floor(j * n / M)
    y[j] <- 1 / (n / M) * sum(x[id1:id2])
  }

  for (j in 1:(M-1)){
    g[j] <- y[M] + 4 * y[j] - 1
  }
  g[M] <- 2 * y[M] - 1

  min_val <- 999999
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
    id1 <- floor((j - 1) * (n / M)) + 1
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
    id1 <- floor((j - 1) * (n / M)) + 1
    id2 <- floor(j * n / M)
    y[j] <- 1 / (n / M) * sum(x[id1:id2]^0.1)
  }

  for (j in 1:(M-1)){
    g[j] <- y[M]^2 + y[j]^2 - 1
  }
  return(g)
}
