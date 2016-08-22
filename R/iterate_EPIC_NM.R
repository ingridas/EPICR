#' Performs one iteration of EPIC algorithm
#'
#' @param wv - a weighting vector
#' @param i - iteration number
#' @param con - constraints, an analytical function cheap to evaluate, in a form g_i(x)>=0, if available; otherwise it is equal to FALSE
#' @param L, U - row vectors of lower and upper bounds of the design space
#' @param p.pareto - indicator used to assign points to one of the set: dominated or nondominated
#' @param p.next - probability used to select a next decision vector
#' @param strategy - strategy used to select the nex vector to evaluate
#' @param iter - a list conatinig information related with the current algorithm iteration
#' @param absmin - the estimate of minimum values of objective functions (if not provided we will use the min value of Y)
#' @param absmax - the estimate of maximum values of objective functions (if not provided we will use the max value of Y)
#' @param local.search - indicates what EPIC version will be run. TRUE - enhanced EPIC version,FALSE -  the original one. DEFAULT = TRUE.
#' @return iter -  an updated list

iterate_EPIC_NM <- function(wv,i,con,L,U,p.next,p.pareto,strategy,iter,absmin,absmax,local.search){

  ns <- nrow(wv)
  nobj <- ncol(iter$Ye)
  # --------------- check the optimization status ----------------------
  if (iter$status == "simplex"){
    # check how many evaluations have been done with a local search
    if (nrow(iter$Xe) - iter$marker >= 5){
      # switch to EPIC
      print(paste('Running EPIC at function evaluation ', (nrow(iter$Xe)+1)))
      iter$marker <- nrow(iter$Xe) #change the marker
      iter$status <- "EPIC"
    }else{
      # continue with simplex
      iter$curent.P <- NULL
      # call function performing the following simplex operation
      iter$ys <- Tcheby_values(iter$y,iter$w,absmin,absmax) # scalarize MOP
      iter <- transform_simplex(iter,L,U,con)
    }
  }

  if (iter$status != "simplex"){
    # find the current Pareto optimal set
    pareto <- pareto_filter(iter$Ye,iter$Xe)
    iter$current.P <- pareto$P
    current.front <- pareto$Yp
    # if EPIC version is selected, no need to run a local search
    if (local.search & (nrow(iter$Xe) - iter$marker >= 5)){
      # find the Pareto front 5 iterations ago
      k <- nrow(iter$Ye)
      pareto <- pareto_filter(iter$Ye[1:(k-5),],iter$Xe[1:(k-5),])
      previous.front <- pareto$Yp
      impr <- check_improvement(previous.front,current.front)
    }else{
      impr <- TRUE # to continue with EPIC
    }
    if (impr){
      # continue with EPIC search
      # Find labels to the initial set: nondominated or dominated
      iter$classE <- matrix("dominated",nrow(iter$Xe),1)
      iter$classE[iter$current.P] <- "nondominated";
      ## ------------ selecting number of folds -------------------
      # first 10 iterations the number of folds equals to the number of vectors to apply leave-one-out crossvalidation
      if (nobj <= 3){
        if (nrow(iter$classE) <= 10){
          v <- nrow(iter$classE)
        } else {
          v <- 10 #v-fold cross validation is applied
        }
      } else {
        v <- 10
      }
      ## ------------------- Train an SVM classifier with the svmtrain function -------------------
      k1 <- length(which(iter$classE=="nondominated"))
      k2 <- length(which(iter$classE=="dominated"))
      w <- k2/k1; # weigth is proportional to the ratio of examples in both classes

      # finding the best C and gamma parameters values for a new classifier and train the SVM with the best values
      classE <- factor(iter$classE)
      param <- best_parameters_cross(iter$Xe,classE,v,w)
      c <-param$c.best
      g <-param$g.best
      model <- svm(x=iter$Xe,y=classE, gamma = g, cost = c, class.weights = c(dominated = 1,nondominated = w), cross = v, probability = TRUE)
      # print(model$tot.accuracy)

      # ------------------ Training set - evaluated vectors ----------------------
      pred <- predict(model,iter$Xe, decision.values = TRUE, probability = TRUE)
      prob <- attr(pred, "probabilities")  # probabilities
      prob.non <- subset(prob,select = nondominated)
      nondom <- which(prob.non >= 0.5)

      # the following is not used anywere, just to check the model prediction
      classP <- matrix("dominated",nrow(prob),1)
      classP[nondom] <- "nondominated"
      p <- table(pred,classE)
      accuracy_class1 = (p[3] / (p[2]+p[3]))*100 # nondominated class
      accuracy_class2 = (p[1] / (p[1]+p[4]))*100 # dominated class
      accuracy_overall = (p[1]+p[3])*100/sum(p)

      ## ------------------ Make prediction on unevaluated vectors ----------------
      pred <- predict(model,iter$Xu, decision.values = TRUE, probability = TRUE)
      prob <- attr(pred, "probabilities")  # probabilities
      prob.non <- subset(prob,select = nondominated)
      nondom <- which(prob.non >= 0.5)
      classP <- matrix("dominated",nrow(prob),1)
      classP[nondom] <- "nondominated"

      ## --------------------- Call selection_strategy function -------------------------------
      # to suggest the next decision vector
      m <- 1 # currently only 1, To Do m > 1!!!!
      I <- selection_strategy(strategy,iter$Xu,iter$Xe,prob.non,p.next,i,m)
      iter$x <- iter$Xu[I,]
      iter$Xu <- iter$Xu[-I,]
      }else{
        # --------------------------------------------------------
        # we start a local search using Neldear-Mead simplex approach
        print(paste('Starting a local search at function evaluation ', (nrow(iter$Xe)+1)))
        # select weigthing vector for scalarizing objective function
        iter$marker <- nrow(iter$Xe)
        iter$kk <- iter$kk %% ns + 1
        iter$w <- wv[iter$kk,] # or we can use wv[kk,]/k, where k = max(Yp)-min(Yp) and we won't need normalization
        # find corresponding scalarized response
        Ys <- Tcheby_values(iter$Ye,iter$w,absmin,absmax)
        # check if there is enough evaluated vectors to compose initial simplex
        if (nrow(iter$Xe) < ncol(iter$Xe)+1){
          stop('There is not enough evaluated vectors to compose a simplex!')
        }
        # call init simplex selection function
        temp <- create_init_simplex(iter$Xe,iter$Ye,Ys,iter$current.P)
        iter$simplex.x <- temp[1:(nvar+1),1:nvar]
        iter$simplex.y <- temp[1:(nvar+1),nvar+1] # scalarized value
        iter$status <- "simplex"
        iter$simplex.status <- "init"
        iter <- transform_simplex(iter,L,U,con)
      }
    }
  return(iter)
}
