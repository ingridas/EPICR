#' @importFrom emoa dominated_hypervolume
NULL

#' This function checks if there is an improvement based on HV metric;
# compares the current Pareto front with the one obtained before 5 evaluations
#'
#' @param previous.front - nondominated solutions in the objective space of the solution set obtained 5 evaluations ago
#' @param current.front - nondominated solutions in the objective space of the current solution set
#' @return TRUE or FALSE

check_improvement <- function (previous.front,current.front){

  # treschold value for HV improvement
  HV.improvement <- 1
  previous.front <- t(previous.front)
  current.front <- t(current.front)
  # compare hypervolume of previous Pareto set and the current one
  allpoints <- cbind(previous.front,current.front)
  # take a reference point as a maximum in each coordinate
  ref <- apply(allpoints, 1, max) + 1
  previous.HV <- emoa::dominated_hypervolume(previous.front, ref)
  current.HV <- emoa::dominated_hypervolume(current.front,ref)
  # calculate the change in HV metric
  HV.change <- 100*(current.HV - previous.HV)/current.HV
  # if HV metric change is lower than threshold value, the improvement is FALSE
  if (HV.change < HV.improvement){
    improvement <- FALSE
    }else{
      improvement <- TRUE
      }
  return(improvement)
  }
