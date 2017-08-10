#' @importFrom stats predict
#' @importFrom utils read.table write.table
#' @importFrom lhs augmentLHS
NULL

#' This is a main function which iterates an EPIC algorithm function, saves data
#'
#' @param maxEval - maximum number of evaluations
#' @param nobj - number of objectives
#' @param nvar - number of decision variables
#' @param L,U - row vectors of lower and upper bounds of the design space
#' @param func - TRUE, if there is availble function, FALSE if we have only a table of evaluated data (experiments)
#' @param local.search - indicates what EPIC version will be run. TRUE - enhanced EPIC version,FALSE -  the original one. DEFAULT = TRUE.
#' @param problem.name - the string of problem corresponding to its function, NULL - if problem function is unavailable
#' @param p.next -  the propbability value used to select a decision vector, defalut value is 0.6
#' @param p.pareto -  indicator used to assign points to one of the set: dominated or nondominated, default value is 0.7
#' @param strategy - strategy id used to select nex vector to evaluate, default value is 1
#' @param repN - a size of decision space representation sample
#' @param absmin - the estimate of minimum values of objective function (if not provided, calculated minimum value is used)
#' @param absmax - the estimate of maximum values of objective functions (if not provided we will use the max value of Y)
#' @param con - constraints, an analytical function cheap to evaluate, in a form g_i(x)>=0, if available; otherwise it is equal to FALSE
#' @param stopping - an indicator used to stop an algorithm at every iteration; TRUE - if algorithm is stopping, otherwise FALSE
#' @return Pareto optimal solutions of evaluated vectors as well as all evaluated solutions.

run_epic_nm <- function(maxEval,nobj,nvar,L,U,func,problem.name = NULL,local.search = TRUE,p.next = 0.6,p.pareto = 0.7,strategy = 1,repN = 500,absmin = NULL,absmax = NULL,con = FALSE,stopping = FALSE){

  if (missing(maxEval) | maxEval<=0)
    stop("Need to specify number of evaluation as a positive integer")

  if (missing(nobj) | nobj<=0)
    stop("Need to specify number of objective as a positive integer")

  if (missing(nvar) | nvar<=0)
    stop("Need to specify number of evaluation as a positive integer")

  if (missing(L) | length(L)<nvar)
    stop("Need to specify lower bound as a vector with length equal to a number of variables")

  if (missing(U) | length(U)<nvar)
    stop("Need to specify upper bound as a vector with length equal to a number of variables")

  if (missing(func) | (!isTRUE(func) & !identical(func, FALSE)))
    stop("Need to specify func value as TRUE or FALSE")

  if (isTRUE(func)){
    if (is.null(problem.name)){
      stop("Need to specify problem function name")
    }else{
      fobj <- get(problem.name) # covert problem to a function writen as string, e.g. "dtz1"
    }
  }

  check.min <- FALSE
  if (is.null(absmin) | length(absmin) < nobj){
    print("You didn't specify correct estimated minimum objective values; minimum values will be used")
    check.min <- TRUE
    }

  check.max <- FALSE
  if (is.null(absmax) | length(absmax) < nobj){
    print("You didn't specify correct estimated maximum objective values; maximum values will be used")
    check.max <- TRUE
    }

  if (is.character(con)){
    con <- get(con)  # if con exist, convert to function
  }else{
    con <- FALSE #constraints do not exist
    print("You didn't spcified constraints")
  }


  np <- 1              # number of points for simultaneous evaluations

  # output files storing the results
  temp.decision.file <- "decision.txt"
  temp.response.file <- "response.txt"
  init.sample.file <- "initsample.txt"

  if (stopping){
    # file for writing and reading intermediate results, when stopping = TRUE
    file.int <- "Intermediate_results.txt"
    if (file.exists(file.int)){
      # read intermediate results from text file
      load(file.int)
      i <- temp[[1]]   # first element in the list
      wv <- temp[[2]]  # second element in the list
      # read evaluated response from the file
      temp.resp <- as.matrix(utils::read.table(temp.response.file, sep = "\t"))
      # and remove that temporary file
      file.remove(temp.response.file)
      colnames(temp.resp) <- NULL
      if (ncol(temp.resp) == nobj){
        iter$y <- temp.resp
        temp.resp <- NULL
        # the response is obtained, so we can update our list
        iter$Xe <- rbind(iter$Xe,iter$x)
        iter$Ye <- rbind(iter$Ye,iter$y)
        do.init <- FALSE
      }else{
        stop("Format of response is incorrect: number of columns doesn't match number of objectives")
        }
      }else{
        # no intermedate reults exist; initialization is required
        do.init <- TRUE
      }
  }else{
    do.init <- TRUE
    }

  ## ---------------------INITIALIZATION ----------------------------------
  if (do.init){
    # -------------- Reading an initial sample from a file ------------------------
    initSample <- as.matrix(utils::read.table(init.sample.file, sep = "\t"))
    colnames(initSample) <- NULL
    initEval <- nrow(initSample)
    # separate X and Y
    Xe <- initSample[,1:nvar]
    Ye <- initSample[,(nvar+1):(nvar+nobj)]

    # -------------- Generate design space representation -------------
    # normalize to interval [0 1]
    initDesign <- apply(Xe, MARGIN = 1, FUN = function(X) ((X-L) / (U-L)))
    initDesign <- t(initDesign)
    # use an adaptive sampling to add new points
    #x <- matrix(NA,nrow = repN, ncol = nvar)
    XX <- lhs::augmentLHS(initDesign,repN)
    # denormalize
    XX <- apply(XX, MARGIN = 1, FUN = function(X) (X * (U-L)+L))
    XX <- t(XX)
    if(is.function(con)){
      # we have to evaluate constraints as well
      Y_con <- apply(XX, MARGIN = 1, FUN = con)
      Y_con <- t(Y_con)
      # check if they are feasible
      if(!all(Y_con>=0)){
        # find infeasible decisions and remove unfeasible ones
        k <- which(apply(Y_con,1,min)<0)
        XX <- XX[-k,]
      }
    }
    # remove already evaluated points,i.e. rows from Xe matching rows in XX, unevaluated set
    k<-prodlim::row.match(data.frame(Xe),data.frame(XX),nomatch = 0)
    # remove evaluated points from the set of unevaluated ones (design space representation)
    Xu <- XX[-k,] # unevaluated designs in the decision space

    # generate the weighting vector
    # the formula is: number of weight vectors = (N+k-1)!/N!(k-1)! where k is number of objectives
    if (nobj == 2){
      N <- 10
      }else if (nobj == 3) {
        N <- 5
        }else if (nobj == 4){
          N <- 4
          } else {
            N <- 3
            }
    wv <- generate_weights(N,nobj)

    #----------Initialization of a list storing all the information require to perform one iteration -------
    # iter is a list storing info related with a current iteration
    iter <- vector("list",0)      # list initialization
    iter$status <- 'init'         # initial iteration of EPIC-NM algorithm
    iter$simplex.status <- NA     # stores the performed simplex operation
    iter$x <- NA                  # stores a decicision to be evaluated
    iter$y <- NA                  # stores a corresponding response
    iter$ys <- NA                 # stores a scalarized objective value
    iter$x.pareto <- NA           # decision predicted to be nondominated at the current iteration
    iter$classE <- NA             # class labels of evaluated vector given for SVM
    iter$marker <- initEval       # marker is used to indicate iteration ID before running EPIC method, it is updated every time we run simplex method
    iter$kk <- 2                  # used to select different weigthing vector
    iter$Ye <- Ye                 # initially evaluated responses
    iter$Xe <- Xe
    iter$Xu <- Xu                 #  unevaluated decision vectors

    i <- 0         # current number of iterations
    }
  T <- maxEval   # amount of iterations to be run
  ## ----------------------- MAIN CYCLE ------------------------------------------
  while (nrow(iter$Xe) < T){  # repeat until max evaluation is reached
    i <- i+1
    print(paste('Iteration : ',i))
    # check the min and max values of responses if they aren't provided
    if(isTRUE(check.min)){
      absmin <- apply(iter$Ye,2,min) - 1
    }
    if(isTRUE(check.max)){
      absmax <- apply(iter$Ye,2,max) + 1
    }
    # performs one iteration of the EPIC-NM
    iter <- iterate_EPIC_NM(wv,i,con,L,U,p.next,p.pareto,strategy,iter,absmin,absmax,local.search)
    # if function is available, evaluate its resoponse, otherwise ask the user to evaluate selected vector
    if (func){
      if (is.array(iter$x)){
        temp.resp <- t(apply(iter$x, MARGIN = 1, FUN = fobj))
        }else{
          temp.resp <- fobj(iter$x)
          }
      }else if(stopping){
        # ask user to evaluate next decision vector
        if (is.array(iter$x)){
          utils::write.table(iter$x,temp.decision.file,sep="\t", col.names = F, row.names = F)
        }else{
          utils::write.table(t(as.matrix(iter$x)),temp.decision.file,sep="\t", col.names = F, row.names = F)
        }
        print(paste('Please evaluate the decision vector(s) written to the file: ',temp.decision.file))
        # save intermediate results and halt the algorithm
        temp <- list(i,wv)
        save(iter,temp,file = file.int)
        print(paste('Write evaluated response(s) to file:',temp.response.file,'and restar the algorithm'))
        return()
        }else{
          print('Please evaluate the following decision vector(s): ')
          print(iter$x)
          # ask to input the response and
          # wait until the response will be typed
          temp <- NULL
          while(is.null(temp)){
            temp <- readline(prompt="Enter the response(s) separated by SPACE: ")
            temp.resp <- as.numeric(unlist(strsplit(temp," ")))
            temp.resp <- matrix(temp.resp,nrow = nobj, byrow = TRUE)
            # if only 1 vector was evaluated, it gives us a matrix of 1 column
            # thus we have to transpose the response to the matric of 1 row
            if(ncol(temp.resp)==1){
              temp.resp <- t(temp.resp)
            }
            Sys.sleep(1)
            }
          }
    iter$y <- temp.resp
    temp.resp <- NULL
    # the response is obtained, so we can update our list
    iter$Xe <- rbind(iter$Xe,iter$x)
    iter$Ye <- rbind(iter$Ye,iter$y)
    }
  # inform user that the max number of function evaluations was reached
  print("Max number of function evaluations was reached!")
  # find Pareto optimal set among evaluated vectors
  pareto <- pareto_filter(iter$Ye,iter$Xe)
  P <- pareto$P
  final.Yep <- iter$Ye[P,]
  final.Xep <- iter$Xe[P,]
  final <- as.matrix(cbind(final.Xep,final.Yep))
  final_all <- as.matrix(cbind(iter$Xe,iter$Ye))
  # and print the results to the file
  utils::write.table(final,file="Results_nondominated_solutions",sep="\t", col.names = F, row.names = F)
  utils::write.table(final_all,file="Results_all_evaluated_solutions",sep="\t", col.names = F, row.names = F)
}
