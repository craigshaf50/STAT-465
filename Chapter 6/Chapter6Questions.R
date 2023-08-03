library(Bolstad)

#6.1
theta<-c(0.0,0.2,0.4,0.6,0.8,1.0)
theta.prior<-c(1/6,1/6,1/6,1/6,1/6,1/6)
#function differs from book
results<-binodp(3,8,uniform=F,pi=theta,pi.prior=theta.prior,ret=T)
results

#a. 
# Conditional distribution of x given pi and  n:
#   
#   0      1      2      3      4      5      6      7      8
# 0   1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
# 0.2 0.1678 0.3355 0.2936 0.1468 0.0459 0.0092 0.0011 0.0001 0.0000
# 0.4 0.0168 0.0896 0.2090 0.2787 0.2322 0.1239 0.0413 0.0079 0.0007
# 0.6 0.0007 0.0079 0.0413 0.1239 0.2322 0.2787 0.2090 0.0896 0.0168
# 0.8 0.0000 0.0001 0.0011 0.0092 0.0459 0.1468 0.2936 0.3355 0.1678
# 1   0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 1.0000

#b.
# Column 4 contains the likelihoods

#c.
# Joint distribution:
#   
#   0      1      2      3      4      5      6      7      8
# [1,] 0.1667 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
# [2,] 0.0280 0.0559 0.0489 0.0245 0.0076 0.0015 0.0002 0.0000 0.0000
# [3,] 0.0028 0.0149 0.0348 0.0464 0.0387 0.0206 0.0069 0.0013 0.0001
# [4,] 0.0001 0.0013 0.0069 0.0206 0.0387 0.0464 0.0348 0.0149 0.0028
# [5,] 0.0000 0.0000 0.0002 0.0015 0.0076 0.0245 0.0489 0.0559 0.0280
# [6,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1667
# found by multiplying the conditional probability Y|pi by the prior prob of pi

#d.
# Marginal distribution of x:
#   
#           0      1      2      3      4      5      6      7      8
# [1,] 0.1975 0.0722 0.0908 0.0931 0.0927 0.0931 0.0908 0.0722 0.1975
# found by summing the probabilities down the columns

#e. 
# posterior probabilities are found by multiplying posterior and likelihood, 
# summing the column, and dividing the product of posterior X likelihood by the sum 
# of the column. This gives the posterior