# MLE-2-sample-proportion-ratio-constraint
Find MLE of 2 sample population proportions when their ratio is known
# p2_mle_ratio function returns the MLE of pi2 based on (x1,n1,x2,n2,theta) when pi1 / pi2 = theta

# n1 = sample size from the first population -- it should be a positive integer
# n2 = sample size from the second population -- it should be a positive integer

# x1 = number of successes in the sample drawn from the first population -- it should be a nonnegative integer less than or equal to n1
# x2 = number of successes in the sample drawn from the second population -- it should be a nonnegative integer less than or equal to n2


# Index populations 1 & 2 such that pi1 / pi2 = theta >= 1

#Examples
p2_mle_ratio(45, 125, 35, 75, 1)
p2_mle_ratio(40, 40, 5, 10, 1.2)
p2_mle_ratio(40, 40, 2, 20, 1.2)
p2_mle_ratio(20, 40, 20, 20, 1.5)
p2_mle_ratio(0, 40, 0, 20, 1.1)
p2_mle_ratio(30, 40, 10, 20, 1.1)
p2_mle_ratio(30, 40, 10, 20, 2)

# To find the standard error of p1 - p2 when p1, p2 are MLE of Pi1 and Pi2 under the constraint that Pi1 / Pi2 = theta
# use SE_theta function
# Example
SE_theta(30, 40, 10, 20, 5)
