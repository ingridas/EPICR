#' @importFrom e1071 svm
NULL

#' Performs cross validation and returns the best parameters values of SVM
#'
#' @param x - classifications data
#' @param y - decision values
#' @param v - number of folds
#' @param w - weight
#' @return The best values of gamma and C

best_parameters_cross <- function (x, y, v, w){

# This function returns the optimal C, gamma values found using the cross validation.
# ###################################################################
# cross validation scale 1
# This is the big scale (rough)
# ###################################################################
step_size <- 1
log2c_list <- seq(-5,15,step_size)
log2g_list <- seq(-15,3,step_size)

m <- length(log2c_list)
n <- length(log2g_list)
C <- matrix(rep(log2c_list,each = n),nrow = n)
gam <- matrix(rep(log2g_list,m),nrow = n)

# grid search, and cross-validation
cv_acc <- 0
bestcv <- -10
for (i in 1:length(C)){
mod <- e1071::svm(x, y, gamma = 2^gam[i], cost = 2^C[i], class.weights = c(dominated = 1,nondominated = w), cross = v, probability = TRUE)
cv_acc <- mod$tot.accuracy
if (cv_acc > bestcv){
  bestcv <- cv_acc
  bestLog2c <- C[i]
  bestLog2g <- gam[i]
  }
}
# ###################################################################
# cross validation scale 2
# This is the medium scale
# ###################################################################
prev_step_size <- step_size
step_size <- prev_step_size/2
log2c_list <- seq(bestLog2c - prev_step_size,bestLog2c + prev_step_size,step_size)
log2g_list <- seq(bestLog2g - prev_step_size,bestLog2g + prev_step_size,step_size)

m <- length(log2c_list)
n <- length(log2g_list)
C <- matrix(rep(log2c_list,each = n),nrow = n)
gam <- matrix(rep(log2g_list,m),nrow = n)


cv_acc <- 0
for (i in 1:length(C)){
  mod <- e1071::svm(x, y, gamma = 2^gam[i], cost = 2^C[i], class.weights = c(dominated = 1,nondominated = w), cross = v, probability = TRUE)
  cv_acc = mod$tot.accuracy
  if (cv_acc > bestcv){
    bestcv <- cv_acc
    bestLog2c <- C[i]
    bestLog2g <- gam[i]
  }
}
# ###################################################################
# cross validation scale 3
# This is the small scale
# ###################################################################

prev_step_size <- step_size
step_size <- prev_step_size/2
log2c_list <- seq(bestLog2c - prev_step_size,bestLog2c + prev_step_size,step_size)
log2g_list <- seq(bestLog2g - prev_step_size,bestLog2g + prev_step_size,step_size)

m <- length(log2c_list)
n <- length(log2g_list)
C <- matrix(rep(log2c_list,each = n),nrow = n)
gam <- matrix(rep(log2g_list,m),nrow = n)


cv_acc <- 0
for (i in 1:length(C)){
  mod <- e1071::svm(x, y, gamma = 2^gam[i], cost = 2^C[i], class.weights = c(dominated = 1,nondominated = w), cross = v, probability = TRUE)
  cv_acc <- mod$tot.accuracy
  if (cv_acc > bestcv){
    bestcv <- cv_acc
    bestLog2c <- C[i]
    bestLog2g <- gam[i]
  }
}

c.best <- 2^bestLog2c
g.best <- 2^bestLog2g
results <- list("c.best" = c.best,"g.best" = g.best)
return(results)
}
