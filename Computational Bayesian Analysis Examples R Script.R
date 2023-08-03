#Computing from Exact Posterior: Beta(6, 16)
a=6
b=16
m = a/(a+b)
v = a*b/((a+b+1)*(a+b)^2)
med = qbeta(0.50, a, b)
l = qbeta(0.025, a, b)
u = qbeta(0.975, a, b)

#Computing from a Sample taken from the Exact Posterior: Beta(6, 16)
n=100
n=500
n=1000
n=10000
n=100000
sample = rbeta(n, a, b)
hist(sample, breaks=25)           
m = mean(sample)
med = median(sample)
v = var(sample)
summary(sample)
CI = quantile(sample, c(0.025, 0.975))


#Inverse Probability Sampling
#Using the exponential distribution ~ exp(lambda=2)
#[Transforming from unif(0,1) to exp(2)]
lambda=2
u = runif(1000, 0, 1)
inv = -log(1-u)/lambda  #inv of CDF of exp

hist(inv, breaks=25, prob=TRUE)
lines(density(inv), col="blue")
curve(dexp(x, rate = lambda), from=0, to=5, add=TRUE, col="red")


#Acceptance-Rejection Sampling
#target posterior density: f~beta(2,2)
a=2
b=2
#candidate/proposal density: g(x)=1
#M=0.25
#M*g(x) = Mg = 0.25
Mg=0.25

accept = numeric(0)
n=10000
i=1
while(i<=n){
  x = runif(1, 0, 1) #g(x) = candidate
  u = runif(1, 0, 1)
  f = x*(1-x) #f is proportional to x(1-x)   
  if(u<f/Mg){
    accept[i] = x
    i=i+1
  }
}
hist(accept, breaks=25, freq=FALSE)
curve(dbeta(x, a, b), from=0, to=1, col="red", add=TRUE)


#Acceptance-Rejection Sampling
#To see acceptance rate:
a=2
b=2
Mg = 0.25

n=100000
accepted = numeric(0)
x.ALL = runif(n, 0, 1) #all candidates from g
j=1
for(i in 1:n){
  u = runif(1, 0, 1)
  f = x.ALL[i]*(1-x.ALL[i]) #f is proportional to x(1-x)  #f = target 
  if(u<f/Mg){
    accepted[j] = x.ALL[i]
    j=j+1
  }
}
acceptance.rate = length(accepted)/n

hist(accept, breaks=25, freq=FALSE)
curve(dbeta(x, a, b), from=0, to=1, col="red", add=TRUE)



#Sampling Importance Resampling
#target posterior density: f~beta(2,8)
a=2
b=8
#importance density: h ~ normal(mean=0.2, var approx 0.015)
n1=1000000
n2=100000
theta = rnorm(n1, 0.2, 0.2) #from importance density # mean & sd

p = theta*(1-theta)^7 #target f proportional to this
h = exp(-((theta-0.2)^2)/(2*0.2)) #importance density
w = p/h #importance weights
w[which(theta<=0 | theta>=1 )] = 0
#
s=sum(w)
norm_w = w/s #normalizing importance weights
#
re_theta = sample(theta, n2, prob=norm_w, replace=TRUE) #resampling
hist(re_theta, breaks=25, freq=FALSE)
curve(dbeta(x, a, b), from=0, to=1, col="red", add=TRUE)


#Importance Sampling to estimate P(Z>5)
#Z~N(0,1)
1-pnorm(5, 0,1) #Exact value

x = rexp(n1, 1)
y = x+5 #importance density = shifted exponential

p = dnorm(y, 0, 1)
h = exp(5-y)
prob.est = sum(p/h)/n1
#No re-sampling here


######################################################
######################################################

#MCMC Metroplis Hastings 
#Candidate density = RW (random walk) 
#Shape of q(theta, theta'):
#q<-exp(-(1/2)*(theta-theta.prime)^2) 

#Posterior is proportional to mixture of: 
#N(0,1), N(3, 0.5^2), N(-3, 0.5^2)
#unscaled posterior g(theta | y):
#g <- 0.7*exp((-1/2)*theta^2)+
#     0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta-3)^2)+
#     0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta+3)^2)

#Initial value of theta = 2

#Graph of initial densities:
theta.prime.all<-seq(-6,6, by = 0.001)

q<-exp(-(1/2)*(2-theta.prime.all)^2) #2 plugged in for theta.prime
plot(theta.prime.all,q, type="l", col="red", lty=2) #shape of candidate density 

g <- 0.7*exp((-1/2)*theta.prime.all^2)+
     0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime.all-3)^2)+
     0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime.all+3)^2)
lines(theta.all,g) #adding shape of posterior mixture density to plot 


#Starting MCMC Algorithm!
#For n=1:5
theta<-numeric(0)

#n=1 ~ Finding theta.1
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q
theta.prime

#Plot (like in book)
plot(theta.prime.all,q, type="l", col="red", lty=2) #Initial q
points(2,0, pch=18, col="red") #Initial theta = 2
lines(theta.prime.all,g) #g
points(theta.prime,0, pch=16, col="orange") #Drawn theta.prime

#Computing acceptance probability for theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*2^2)+      #Initial theta = 2
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(2-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(2+3)^2)

acceptance.prob = min(1, g.theta.prime/g.theta)
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta[1]<-theta.prime #theta.prime becomes theta_1
} else{
  theta[1] <- 2 #No move ~ theta.1 = initial theta (theta.0) = 2
}
theta

#Compare new and previous q's:
plot(theta.prime.all,q, type="l", col="red", lty=2) #shape of candidate density when n=0
q<-exp(-(1/2)*(theta[1]-theta.prime.all)^2) #new q using theta.1
lines(theta.prime.all,q, type="l", col="blue", lty=3) #shape of new q
lines(theta.all,g) #adding shape of posterior mixture density to plot 


#n=2 ~ Finding theta.2
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q (using new q based off of theta.1)
theta.prime

#Plot (like in book)
plot(theta.prime.all,q, type="l", col="blue", lty=3) #New q using theta.1
points(theta[1],0, pch=18, col="blue") #theta.1
lines(theta.prime.all,g) #g
points(theta.prime,0, pch=16, col="orange") #Drawn theta.prime

#Computing acceptance probability for new theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta[1]^2)+      #theta.1, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[1]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[1]+3)^2)

acceptance.prob = min(1, g.theta.prime/g.theta)
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta[2]<-theta.prime #theta.prime becomes theta.2
} else{
  theta[2] <- theta[1] #No move 
}
theta

#Compare new and previous q's:
plot(theta.prime.all,q, type="l", col="blue", lty=3) #shape of candidate density when n=1
q<-exp(-(1/2)*(theta[2]-theta.prime.all)^2) #new q using theta.2
lines(theta.prime.all,q, type="l", col="salmon", lty=4) #shape of new q
lines(theta.all,g) #adding shape of posterior mixture density to plot 


#n=3 ~ Finding theta.3
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q (using new q based off of theta.2)
theta.prime

#Plot (like in book)
plot(theta.prime.all,q, type="l", col="salmon", lty=4) #New q using theta.2
points(theta[2],0, pch=18, col="salmon") #theta.2
lines(theta.prime.all,g) #g
points(theta.prime,0, pch=16, col="orange") #Drawn theta.prime

#Computing acceptance probability for new theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta[2]^2)+      #theta.2, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[2]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[2]+3)^2)

acceptance.prob = min(1, g.theta.prime/g.theta)
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta[3]<-theta.prime #theta.prime becomes theta.3
} else{
  theta[3] <- theta[2] #No move 
}
theta

#Compare new and previous q's:
plot(theta.prime.all,q, type="l", col="salmon", lty=4) #shape of candidate density when n=2
q<-exp(-(1/2)*(theta[3]-theta.prime.all)^2) #new q using theta.3
lines(theta.prime.all,q, type="l", col="green", lty=5) #shape of new q
lines(theta.all,g) #adding shape of posterior mixture density to plot 


#n=4 ~ Finding theta.4
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q (using new q based off of theta.3)
theta.prime

#Plot (like in book)
plot(theta.prime.all,q, type="l", col="green", lty=5) #New q using theta.3
points(theta[3],0, pch=18, col="green") #theta.3
lines(theta.prime.all,g) #g
points(theta.prime,0, pch=16, col="orange") #Drawn theta.prime

#Computing acceptance probability for new theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta[3]^2)+      #theta.3, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[3]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[3]+3)^2)

acceptance.prob = min(1, g.theta.prime/g.theta)
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta[4]<-theta.prime #theta.prime becomes theta.4
} else{
  theta[4] <- theta[3] #No move 
}
theta

#Compare new and previous q's:
plot(theta.prime.all,q, type="l", col="green", lty=5) #shape of candidate density when n=3
q<-exp(-(1/2)*(theta[4]-theta.prime.all)^2) #new q using theta.4
lines(theta.prime.all,q, type="l", col="purple", lty=6) #shape of new q
lines(theta.all,g) #adding shape of posterior mixture density to plot 


#n=5 ~ Finding theta.5
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q (using new q based off of theta.3)
theta.prime

#Plot (like in book)
plot(theta.prime.all,q, type="l", col="purple", lty=6) #New q using theta.4
points(theta[4],0, pch=18, col="purple") #theta.3
lines(theta.prime.all,g) #g
points(theta.prime,0, pch=16, col="orange") #Drawn theta.prime

#Computing acceptance probability for new theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta[4]^2)+      #theta.4, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[4]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[4]+3)^2)

acceptance.prob = min(1, g.theta.prime/g.theta)
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta[5]<-theta.prime #theta.prime becomes theta.5
} else{
  theta[5] <- theta[4] #No move 
}
theta

#"Trace Plot":
plot(c(2,theta), type="l")
#Histogram
hist(c(2,theta))


#Continuing...... n=6:1000
#and then to n=20000
#and then to n=50000
q<-exp(-(1/2)*(theta[5]-theta.prime.all)^2) #new q using theta.5
for(n in 6:1000){
#for(n in 1001:20000){
#for(n in 20001:50000){
  theta.prime<-sample(theta.prime.all, 1, prob=q) 

  g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
    0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
    0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

  g.theta <- 0.7*exp((-1/2)*theta[n-1]^2)+      
   0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[n-1]-3)^2)+
   0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta[n-1]+3)^2)

  acceptance.prob = min(1, g.theta.prime/g.theta)

  u<- runif(1,0,1)
  if(u<acceptance.prob){
    theta[n]<-theta.prime 
  } else{
    theta[n] <- theta[n-1] 
  }
  
  q<-exp(-(1/2)*(theta[n]-theta.prime.all)^2) 
}

#"Trace Plot":
plot(c(2,theta), type="l")
#Histogram
hist(c(2,theta))

hist(c(2,theta), freq=FALSE, ylim=c(0,0.7))
lines(theta.all,g) #Unscaled


############################################
#MCMC Metroplis Hastings 
#Independent candidate density = N(0,2.5^2) 
#Shape of q(theta'):
#q<-exp(-(1/(2*2.5^2))*theta.prime^2) 

#Posterior is proportional to mixture of: 
#N(0,1), N(3, 0.5^2), N(-3, 0.5^2)
#unscaled posterior g(theta | y):
#g <- 0.7*exp((-1/2)*theta^2)+
#  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta-3)^2)+
#  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta+3)^2)


#Graph of initial densities - show the N(0,2.5^2) covers the target:
theta.prime.all<-seq(-8,8, by = 0.001)
 
q<-exp(-(1/(2*2.5^2))*(theta.prime.all)^2)
plot(theta.prime.all,q, type="l", col="red", lty=2) #shape of candidate density 

g <- 0.7*exp((-1/2)*theta.prime.all^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime.all-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime.all+3)^2)
lines(theta.prime.all,g) #adding shape of posterior mixture density to plot 

points(2,0, pch=15, col="hotpink") #Add initial theta = 2

#Starting MCMC Algorithm!
#For n=1:5
theta2<-numeric(0)

#n=1 ~ Finding theta.1
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q 
theta.prime
points(theta.prime,0, pch=16, col="green") #Add drawn theta.prime to plot


#Computing acceptance probability for theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*2^2)+      #Initial theta = 2
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(2-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(2+3)^2)


q.theta.prime<-exp(-(1/(2*2.5^2))*theta.prime^2)

q.theta<-exp(-(1/(2*2.5^2))*(2)^2) #Initial theta = 2

acceptance.prob = min(1, (g.theta.prime*q.theta)/(g.theta*q.theta.prime))
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta2[1]<-theta.prime #theta.prime becomes theta.1
} else{
  theta2[1] <- 2 #No move ~ theta.1 = initial theta (theta.0) = 2
}
theta2

#n=2 ~ Finding theta.2
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q {no change in q} 
theta.prime
points(theta.prime,0, pch=17, col="blue") #Add drawn theta.prime to plot


#Computing acceptance probability for theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta2[1]^2)+      #theta.1, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[1]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[1]+3)^2)

q.theta.prime<-exp(-(1/(2*2.5^2))*theta.prime^2) 

q.theta<-exp(-(1/(2*2.5^2))*theta2[1]^2) #theta.1, most recent theta

acceptance.prob = min(1, (g.theta.prime*q.theta)/(g.theta*q.theta.prime))
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta2[2]<-theta.prime #theta.prime becomes theta.2
} else{
  theta2[2]<-theta2[1] #No move ~ theta.2 = theta.1
}
theta2
 

#n=3 ~ Finding theta.3
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q {no change in q} 
theta.prime
points(theta.prime,0, pch=18, col="red") #Add drawn theta.prime to plot


#Computing acceptance probability for theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta2[2]^2)+      #theta.2, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[2]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[2]+3)^2)

q.theta.prime<-exp(-(1/(2*2.5^2))*theta.prime^2) 

q.theta<-exp(-(1/(2*2.5^2))*theta2[2]^2) #theta.2, most recent theta

acceptance.prob = min(1, (g.theta.prime*q.theta)/(g.theta*q.theta.prime))
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta2[3]<-theta.prime #theta.prime becomes theta.3
} else{
  theta2[3]<-theta2[2] #No move ~ theta.3 = theta.2
}
theta2


#n=4 ~ Finding theta.4
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q {no change in q} 
theta.prime
points(theta.prime,0, pch=19, col="purple") #Add drawn theta.prime to plot


#Computing acceptance probability for theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta2[3]^2)+      #theta.3, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[3]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[3]+3)^2)

q.theta.prime<-exp(-(1/(2*2.5^2))*theta.prime^2) 

q.theta<-exp(-(1/(2*2.5^2))*theta2[3]^2) #theta.3, most recent theta

acceptance.prob = min(1, (g.theta.prime*q.theta)/(g.theta*q.theta.prime))
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta2[4]<-theta.prime #theta.prime becomes theta.4
} else{
  theta2[4]<-theta2[3] #No move ~ theta.4 = theta.3
}
theta2


#n=5 ~ Finding theta.5
theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q {no change in q} 
theta.prime
points(theta.prime,0, pch=8, col="salmon") #Add drawn theta.prime to plot


#Computing acceptance probability for theta.prime
g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)

g.theta <- 0.7*exp((-1/2)*theta2[4]^2)+      #theta.4, most recent theta
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[4]-3)^2)+
  0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[4]+3)^2)

q.theta.prime<-exp(-(1/(2*2.5^2))*theta.prime^2) 

q.theta<-exp(-(1/(2*2.5^2))*theta2[4]^2) #theta.4, most recent theta

acceptance.prob = min(1, (g.theta.prime*q.theta)/(g.theta*q.theta.prime))
acceptance.prob

#Accept/Reject Step
u<- runif(1,0,1)
if(u<acceptance.prob){
  theta2[5]<-theta.prime #theta.prime becomes theta.5
} else{
  theta2[5]<-theta2[4] #No move ~ theta.5 = theta.4
}
theta2

#"Trace Plot":
plot(c(2,theta2), type="l")
#Histogram
hist(c(2,theta2))


#Continuing...... n=6:2500
#and then to n=10000
#and then to n=30000
for(n in 6:2500){
#for(n in 2501:10000){
#for(n in 10001:30000){
  theta.prime<-sample(theta.prime.all, 1, prob=q) #Draw theta.prime from q {no change in q} 

  g.theta.prime <- 0.7*exp((-1/2)*theta.prime^2)+
    0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime-3)^2)+
    0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta.prime+3)^2)
  
  g.theta <- 0.7*exp((-1/2)*theta2[n-1]^2)+      
    0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[n-1]-3)^2)+
    0.15*(1/0.50)*exp((-1/(2*0.5^2))*(theta2[n-1]+3)^2)
  
  q.theta.prime<-exp(-(1/(2*2.5^2))*theta.prime^2) 
  
  q.theta<-exp(-(1/(2*2.5^2))*theta2[n-1]^2) 
  
  acceptance.prob = min(1, (g.theta.prime*q.theta)/(g.theta*q.theta.prime))

  u<- runif(1,0,1)
  if(u<acceptance.prob){
    theta2[n]<-theta.prime 
  } else{
    theta2[n]<-theta2[n-1] #No move 
  }
}

#"Trace Plot":
plot(c(2,theta2), type="l")
#Histogram
hist(c(2,theta2))

hist(c(2,theta2), freq=FALSE, ylim=c(0,0.7))
lines(theta.prime.all,g) #Unscaled


