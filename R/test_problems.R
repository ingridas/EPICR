#' Set of multiobjective benchamrk test including ZDT test functions, OKA2, Viennet and TF4 (a constrained problem)
#' @param x - is a decision vector of lenght \code{n}
#' @return function output as a vector of lenght 2 or 3
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
