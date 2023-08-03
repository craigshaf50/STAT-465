
# initial value of theta
theta<-0

# Step 1: vertical slice - drawing sample aux variable g(y-axis)
c<-exp(-(theta^2)/2)
c
g.1<-runif(1,0,c)
g.1

# Step 2: horizontal slice - drawing sample theta(x-axis)
gi<-sqrt(-2*log(g.1))
gi
a <- -gi
b <- gi
theta.1<-runif(1,a,b)
theta.1

#continue
c<-exp(-(theta.1^2)/2)
c
g.2<-runif(1,0,c)
g.2

gi<-sqrt(-2*log(g.2))
gi
a <- -gi
b <- gi
theta.2<-runif(1,a,b)
theta.2
