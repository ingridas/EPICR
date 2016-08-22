#' Selects a decision vector from a set representing the decision space
#' @param  s - strategy ID, from 1 to 4
#' @param  Xu - unevaluated decision vectors
#' @param  Xe - evaluated decision vectors
#' @param prob.non - probabilities of unevaluated vectors to be nondominated
#' @param p.next - probability of nondominance the candidate vector should have to be selected to evaluate next
#' @param iter - iteration number
#' @param m - number of vectors to be selected selected for evaluation
#' @return I - index of the selectec decision vector from Xu
##
selection_strategy <- function (s,Xu,Xe,prob.non,p.next,iter,m){
  switch(s,
         # Strategy 1: select vectors whith probability closest to p.next
         {
           tmp <- sort(abs(prob.non - p.next),index.return=TRUE); # sorting the propbabilities according to closest to selected p.next
           idx <- tmp$ix
           if (length(idx)> m){
             I <- idx[1:m]
           } else{
             I <- idx
           }
         },
         # Strategy 2: Having evaluated vectors, we have to select one vector from unevaluated which is the farthest apart from the rest
         # and closest to the probability p.next
         # First check the range of probabilities to be nondominated
         {
           I = which(prob.non >= p.next)
           if (length(I) < 1){
             tmp <- sort(abs(prob.non - p.next),index.return=TRUE); # sorting the propbabilities according to closest to selected p.next
             idx <- tmp$ix
             mm <- ceiling(length(I)/10)
             I <- idx[1:mm]
           } else {
             d <- matrix(, nrow = length(I), ncol = nrow(Xe))
             for (j in 1:length(I)){
               for (i in 1:nrow(Xe)){
                 # calculate euclidean distance
                 d[j,i] <- sqrt(sum((Xe[i,] - Xu[I[j],])^2))
                # or use the following
                # d[j,i] <- dist(rbind(Xe[i,],Xu[I[j],]))
               }
             }
             D <- rowMeans(d)
             k <- which.max(D)
             I <- I[k]
           }

         },
         # Strategy 3:  as Strategy 2 except that every 5th iteration we select a random vector
         {
           if (iter %% 5 == 0){
             I <- sample(x = 1:nrow(Xu), size = 1,replace = FALSE)
           }
           else {
             # ----- Strategy 2 ----------------------
             I <- which(prob.non >= p.next)
             if (length(I) < 1){
               # sorting the propbabilities according to closest to selected p.next
               tmp <- sort(p.next - prob.non,index.return=TRUE)
               idx <- tmp$ix
               mm <- ceiling(length(I)/10)
               I <- idx[1:mm]
             } else {
               d <- matrix(, nrow = length(I), ncol = nrow(Xe))
               for (j in 1:length(I)){
                 for (i in 1:nrow(Xe)){
                   d[j,i] <- sqrt(sum((Xe[i,] - Xu[I[j],])^2))
                   # or use the following
                   # d[j,i] <- dist(rbind(Xe[i,],Xu[I[j],]))
                 }
               }
               D <- rowMeans(d)
               k <- which.max(D)
               I <- I[k]
             }

             # ---------------------------------------
           }
         },
         {
           # Strategy 4: At the very first iteration we select any random vector which has probability of nondominance higher than pnext.
           # Later, we select randomly a vector with a probability from the interval [pnext ptr]. Probability threshold decreases every
           # 10th iteration by 0.1. If there are no such vectors, we select one with closest probability to pnext.
            diff <- 1 - p.next
            step <- diff/10
            p_tr <- seq(1,1-10*step,-step)
            i <- trunc(iter/5)
            if (i < 11 ) {
              II <- which (p.next <= prob.non & prob.non <= p_tr[i+1])
              if (length(II)>0){
              I <- II[sample(x = 1:length(I), size = 1,replace = FALSE)]
              } else {
                idx <- which.min(abs(p.next-prob.non))
                I <- idx
              }
            } else {
              idx <- which.min(abs(p.next-prob.non))
              I <- idx
            }
         },
{
  print('There is no such strategy!')
}
  )
  return(I)
}
