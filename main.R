# p2_mle_ratio function returns the MLE of pi2 based on (x1,n1,x2,n2,theta) when pi1 / pi2 = theta

# n1 = sample size from the first population -- it should be a positive integer
# n2 = sample size from the second population -- it should be a positive integer

# x1 = number of successes in the sample drawn from the first population -- it should be a nonnegative integer less than or equal to n1
# x2 = number of successes in the sample drawn from the second population -- it should be a nonnegative integer less than or equal to n2


# Index populations 1 & 2 such that pi1 / pi2 = theta >= 1
# MLE by formula in closed form case - failing which by method of bisection
# convergence in typical case is set below

con_acc <- 10^(-10)	# in terms of accuracy of MLE p

p2_mle_ratio <- function (x1,n1,x2,n2,theta) {
  
  dlogL <- function(p){
    sum <- 0
    if (x1 + x2> 0) { sum <- sum + (x1+x2)/p}
    if (n1 - x1 > 0)  {  sum <- sum - (n1-x1)* theta/(1-p * theta) }
    if (n2 - x2 > 0) { sum <- sum -  (n2-x2)/ (1-p) }
    return(sum)
  }   # end dlogL function
  
  L <- 0
  U <- 1 / theta   
  
  if( theta == 1 ) { mle <-  (x1+x2)/(n1+n2) } 
  
  else if (x1-n1 == 0)  { 		#Case 2
    mle <- min(U, (n1+ x2)/(n1+n2) ) }
  
  else if ( (x1 - n1 < 0) & (x2 - n2 ==  0) ) { mle <- (x1+ n2)/ (theta * (n1+n2)  ) } #Case 3
  
  
  else if ( (x1 == 0) & (x2 == 0) ) { mle <- 0 }  	# Case 4
  
  
  else {   # all other  cases where nonlinear equation is to be solved
    # begin else loop 1
    
    pt <- matrix(seq(L, U, by=0.01 *U),ncol=1)
    val <- lapply(pt, dlogL)
    val_neg <- val[val <0]
    val_pos <- val[val > 0]
    pt_pos <- pt[val >0 ]
    pt_neg <- pt[val <0 ]
    
    min_pos <- pt_pos[which.min(val_pos)]
    max_neg <- pt_neg[which.max(val_neg)]
    val_min_pos <- val_pos[which.min(val_pos)]
    val_max_neg <- val_neg[which.max(val_neg)]
    
    
    U <- min_pos
    L <- max_neg
    
    
    while (abs(U-L) > con_acc )
    {	#begin while loop 1
      est <- (L+U)/2
      v <- dlogL(est)
      if (v > 0) { U <- est}
      if (v <= 0) { L <- est}
    }	#end while loop 1
    
    est <- (L+U)/2
    mle <- est 
    # break 
    
    #}  	# end else loop 2
    
    
  }		# end else loop 1 # all other  cases where nonlinear equation is to be solved
  
  
  p2_mle <- mle
  p1_mle <- mle * theta
  
  out <- list(p1=p1_mle, p2=p2_mle)
  
  return( out )
  
}		# end p2_mle_ratio



#############################

#########################
# This function calculates the estimate of standard error of p1 - p2 where p1, p2 are restricted MLE of pi1 and pi2 under the binding constraint pi1 / pi2 = theta

SE_theta <- function (x1,n1,x2,n2,delta) {
  
  p2 <- p2_mle_ratio(x1,n1,x2,n2,delta)$p2
  p1 <-  p2_mle_ratio(x1,n1,x2,n2,delta)$p1
  
  SE <- ( p1 * (1-p1)/n1 + p2 *(1-p2)/n2  )^0.5
  
  return(SE) 
}
#######


