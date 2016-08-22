#' Generates a simplex from the evaluated decicion vectors and coresponding scalarized objective value.
#'
#' @param Xe - an k x N matrix of evaluated vectors in the decison space, where N is a number of decision variables
#' @param Ye - an k x M matrix of evaluated vectors in the objective space, where M is a number of objective functions
#' @param Ys - a vector of lenght k containing scalarized objective values
#' @param P - a list returned by function pareto_filter
#' @return initial simplex used in a local search; it's a matrix (N+1) x (N+1), where rows correspond to vertices, first N columns are decision variables, last column - scalarized function value

create_init_simplex <- function(Xe,Ye,Ys,P) {
  nvar <- ncol(Xe)
  # if number of evaluated vectors matches the number of simplex vertices, use them all to form a simplex
  if(nrow(Xe) == nvar+1){
    start.x <- Xe
    start.y <- Ys
  }else{
    nond1 <- sum(P)
    Yp <- Ye[P,]
    Xp <- Xe[P,]
    Ysp <- Ys[P]
    Yd <- Ye[!P,]
    Xd <- Xe[!P,]
    Ysd <- Ys[!P]
    # sort the values from smallest (best for this scalarization) to largest
    if (length(Ysp)>1){
      tmp <- sort(Ysp,index.return=TRUE)
      Ysp <- tmp$x
      id <- tmp$ix
      Xp <- Xp[id,]
      Yp <- Yp[id,]
      }else{
        Yp <- t(as.matrix(Yp))
        Xp <- t(as.matrix(Xp))
        }
    # --------------- select an initial simplex ---------------------
    # compose initial simplex from nondominated vectors
    if (nrow(Xp) >= nvar+1) {
      start.x <- Xp[1:(nvar+1),]
      start.y <- Ysp[1:(nvar+1)]
      } else {
        # if number of nondominated solution is smaller than number of simplex vertices,
        # find 2nd nondominated set
        pareto <- pareto_filter(Yd,Xd)
        P2 <- pareto$P
        nond2 <- sum(P2)
        Yp2 <- Yd[P2,]
        Xp2 <- Xd[P2,]
        Ysp2 <- Ysd[P2]
        # sort scalarized problem output if its length > 1
        if (length(Ysp2)>1){
          tmp <- sort(Ysp2,index.return=TRUE)
          Ysp2 <- tmp$x
          id <- tmp$ix
          Xp2 <- Xp2[id,]
          Yp2 <- Yp2[id,]
          }else{
            # convert vector to matrix and transpose
            Yp2 <- t(as.matrix(Yp2))
            Xp2 <- t(as.matrix(Xp2))
            }
        if (nond1+nond2 >= nvar+1){
          start.x <- rbind(Xp[1:nond1,], Xp2[1:(nvar+1-nond1),])
          start.y <- c(Ysp[1:nond1], Ysp2[1:(nvar+1-nond1)])
          } else {
            # of nondominated is smaller than # of simplex vertices; find 3rd nondominated set
            Yd2 <- Yd[!P2,]
            Xd2 <- Xd[!P2,]
            Ysd2 <- Ysd[!P2]
            pareto <- pareto_filter(Yd2,Xd2)
            P3 <- pareto$P
            nond3 <- sum(P3)
            Yp3 <- Yd2[P3,]
            Xp3 <- Xd2[P3,]
            Ysp3 <- Ysd2[P3]
            if (length(Ysp3)>1){
              tmp <- sort(Ysp3,index.return=TRUE)
              Ysp3 <- tmp$x
              id <- tmp$ix
              Xp3 <- Xp3[id,]
              Yp3 <- Yp3[id,]
              }else{
                Yp3 <- t(as.matrix(Yp3))
                Xp3 <- t(as.matrix(Xp3))
                }
            if (nond+nond2+nond3 >= nvar+1){
              start.x <- rbind(Xps[1:nond1,], Xps2[1:nond2,], Xps3[1:(nvar+1-nond1-nond2),])
              start.y <- c(Ysp[1:nond1], Ysp2[1:nond2], Ysp3[1:(nvar+1-nond1-nond2)])
              }else{
                dom <- nvar+1-nond1-nond2-nond3
                Yd3 <- Yd2[-P3,]
                Xd3 <- Xd2[-P3,]
                Ysd3 <- Ysd2[-P3]
                tmp <- sort(Ysd3,index.return=TRUE)
                Ysd3 <- tmp$x
                id <- tmp$ix
                Xd3 <- Xd3[id,]
                Yd3 <- Yd3[id,]
                start.x <- rbind(Xp[1:nond1,], Xp2[1:nond2,], Xp3[1:nond3,], Xd3[1:dom,])
                start.y <- c(Ysp[1:nond1], Ysp2[1:nond2], Ysp3[1:nond3], Ysd3[1:dom])
              }
          }
      }
    }
  # combining to a matrix
  return(cbind(start.x,as.matrix(start.y)))
}
