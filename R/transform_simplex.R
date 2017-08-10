#' Updates the simplex with obtained new values and proceed with the following simplex operation untill
#' the new decision vector is generated and requires a function evaluation
#'
#' @param iter - a list conatining all the information related with the current algorithm iteration
#' @param L - row vectors of lower bounds of the design space
#' @param U - row vectors of upper bounds of the design space
#' @param con - constraints, an analytical function cheap to evaluate, in a form g_i(x)>=0, if available; otherwise it is equal to FALSE
#' @return an updated list with the information about the current EPIC algorithm iteration

transform_simplex <- function(iter,L,U,con){
  alpha <- 1
  beta <- 0.5
  gamma <- 2
  delta <- 0.5
  large <- 999999999
  # check what operation was performed
  while(1){
  if ((iter$simplex.status == "init") | (iter$simplex.status == "done")){
    if((iter$simplex.status == "done")&(nrow(iter$simplex.x) > length(iter$simplex.y))){
      # add new calculated verteces values after simplex shrinkage
      iter$simplex.y <- c(iter$simplex.y,iter$ys)
    }
    # calculate centroid
    iter <- find_centroid(iter)
    #---------------- calculate reflection through the centroid -------------------
    ref.p <- iter$simplex.centroid +(iter$simplex.centroid - iter$simplex.x[iter$simplex.w,])*alpha
    # check if other constraints exist
    if (is.function(con)){
      # check all constraints
      if (is_feasible(ref.p,L,U) & is_feasible_con(con,ref.p)){
        iter$simplex.ref.x <- ref.p
        iter$simplex.status <- "reflection"
        iter$x <- ref.p
        iter$ys <- NULL # we need to evaluate it
        iter$simplex.ref.y <- NULL # we need to evaluate it
        return(iter)
      }else{
        # variable is infeasible, assign a large value
        iter$simplex.ref.x <- ref.p
        iter$simplex.status <- "reflection"
        iter$x <- ref.p
        iter$ys <- large # we assigned a large number
        iter$simplex.ref.y <- large # we assigned a large number
      }
    }else{
      # check only box constraints
      if (is_feasible(ref.p,L,U)){
        iter$simplex.ref.x <- ref.p
        iter$simplex.status <- "reflection"
        iter$x <- ref.p
        iter$ys <- NULL # we need to evaluate it
        iter$simplex.ref.y <- NULL # we need to evaluate it
        return(iter)
      }else{
        # variable is infeasible, assign a large value
        iter$simplex.ref.x <- ref.p
        iter$simplex.status <- "reflection"
        iter$x <- ref.p
        iter$ys <- large # we assigned a large number
        iter$simplex.ref.y <- large # we assigned a large number
      }
    }
    # ------------ end of centroid calculation ------------------------------
    }else if(iter$simplex.status == "reflection"){
      #check if reflection point is evaluated
      if(is.null(iter$simplex.ref.y)){
        iter$simplex.ref.y <- iter$ys
      }
      # compare the reflection point value with the current simplex vertices
      # If fb<=fr<fs, accept reflection point
      if(all(iter$simplex.y[iter$simplex.b] <= iter$simplex.ref.y) & all(iter$simplex.ref.y < iter$simplex.y[iter$simplex.s])){
        # update the current simplex by replacing the worst vertice by a new one
        iter <- update_simplex(iter,iter$x,iter$ys)
      }else if(all(iter$simplex.y[iter$simplex.b] > iter$simplex.ref.y)){
        #If fr<fb, compute the expansion point
        exp.p <- iter$simplex.centroid +(iter$x - iter$simplex.centroid)*gamma
        iter$simplex.status <- "expansion"
        # check if other constraints exist
        if (is.function(con)){
          # check all constraints
          if (is_feasible(exp.p,L,U) & is_feasible_con(con,exp.p)){
            iter$simplex.exp.x <- exp.p
            iter$x <- exp.p
            iter$ys <- NULL # we need to evaluate it
            iter$simplex.exp.y <- NULL # we need to evaluate it
            return(iter)
          }else{
            iter$simplex.exp.x <- exp.p
            iter$x <- exp.p
            iter$ys <- large # we need to evaluate it
            iter$simplex.exp.y <- large # we assigned a large number
          }
        }else{
          if (is_feasible(exp.p,L,U)){
            iter$simplex.exp.x <- exp.p
            iter$x <- exp.p
            iter$ys <- NULL # we need to evaluate it
            iter$simplex.exp.y <- NULL # we need to evaluate it
            return(iter)
          }else{
            iter$simplex.exp.x <- exp.p
            iter$x <- exp.p
            iter$ys <- large
            iter$simplex.exp.y <- large  # we assigned a large number
          }
        }
      }else if(all(iter$simplex.y[iter$simplex.s] <= iter$simplex.ref.y)){
        # If fr>= fs, calculate an contraction point by using the better of the two points xw and xr
        if(all(iter$simplex.y[iter$simplex.w] > iter$simplex.ref.y)){
          # If fs<=fr<fw, compute an outside contraction point
          contr.p <- iter$simplex.centroid +(iter$simplex.ref.x - iter$simplex.centroid)*beta
          iter$simplex.status <- "outside contraction"
        }else{
          #If fr>=fw, compute an inside contraction point
          contr.p <- iter$simplex.centroid +(iter$simplex.x[iter$simplex.w,] - iter$simplex.centroid)*beta
          iter$simplex.status <- "inside contraction"
        }
        if (is.function(con)){
          # check all constraints
          if (is_feasible(contr.p,L,U) & is_feasible_con(con,contr.p)){
            iter$simplex.contr.x <- contr.p
            iter$x <- contr.p
            iter$ys <- NULL # we need to evaluate it
            iter$simplex.contr.y <- NULL # we need to evaluate it
            return(iter)
          }else{
            iter$simplex.contr.x <- contr.p
            iter$x <- contr.p
            iter$ys <- large
            iter$simplex.contr.y <- large # we assigned a large number
            }
        }else{
          if (is_feasible(contr.p,L,U)){
            iter$simplex.contr.x <- contr.p
            iter$x <- contr.p
            iter$ys <- NULL # we need to evaluate it
            iter$simplex.contr.y <- NULL # we need to evaluate it
            return(iter)
          }else{
            iter$simplex.contr.x <- contr.p
            iter$x <- contr.p
            iter$ys <- large
            iter$simplex.contr.y <- large # we assigned a large number
            }
          }
      } # ------------ end of reflection -------------------------------------------
    }else if(iter$simplex.status == "expansion"){
      # check if an expansion point is evaluated
      if(is.null(iter$simplex.exp.y)){
        iter$simplex.exp.y <- iter$ys # assign a scalarized objective value
      }
      # Compare expantion point value with the reflection
      # If fe<fr, accept xe and terminate the iteration.
      if(iter$simplex.exp.y < iter$simplex.ref.y){
        iter <- update_simplex(iter,iter$simplex.exp.x,iter$simplex.exp.y)
      }else{
        # Otherwise (if fe>=fr), accept xr and terminate the iteration
        iter <- update_simplex(iter,iter$simplex.ref.x,iter$simplex.ref.y)
      } # ------------ end of expansion --------------------------
    }else if(iter$simplex.status == "outside contraction"){
      # check if the value to contraction point is evaluated
      if (is.null(iter$simplex.contr.y)){
        iter$simplex.contr.y <- iter$ys
      }
      #If f_cont<=fr, accept x_contr
      if(iter$simplex.contr.y <= iter$simplex.ref.y){
        iter <- update_simplex(iter,iter$simplex.contr.x,iter$simplex.contr.y)
      }else{
        # perform a shrink transformation
        iter <- shrink_simplex(iter,delta,con,L,U)
        return(iter)
      }
      # ------------ end of outside contraction -------------------------
    }else if(iter$simplex.status == "inside contraction"){
      # check if the value to contraction point is evaluated
      if (is.null(iter$simplex.contr.y)){
        iter$simplex.contr.y <- iter$ys
      }
     # If fcontr<fw
      if(all(iter$simplex.contr.y < iter$simplex.y[iter$simplex.w])){
        iter <- update_simplex(iter,iter$simplex.contr.x,iter$simplex.contr.y)
      }else{
        # perform shrink transformation
        iter <- shrink_simplex(iter,delta,con,L,U)
        return(iter)
      }
     # ------------ end of inside contraction -------------------------
    }
  }
}
