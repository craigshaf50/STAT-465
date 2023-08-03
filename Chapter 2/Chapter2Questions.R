#install.packages('Bolstad')
library(Bolstad)
library(tidyverse)
help(sscsample.data)


#2.1
#a
boxplot(income~ethnicity, data= sscsample.data)
sscsample.data %>% 
  group_by(ethnicity) %>% summarize(AvgIncome=mean(income))

#b
mySamples = list(simple=NULL, strat=NULL, cluster=NULL)
mySamples$simple = sscsample(20,200)

colMeans(mySamples$simple$s.strata)

mean(mySamples$simple$means)
mean(sscsample.data$income)
#9.040897- 8.994273 = 0.046624

#c
mySamples$strat = sscsample(20,200,"stratified")
mySamples$strat

colMeans(mySamples$strat$s.strata)

mean(mySamples$strat$means)
mean(sscsample.data$income)
#9.031885-8.994273 = 0.037612

#d
mySamples$cluster = sscsample(20,200,"cluster")
mySamples$cluster

colMeans(mySamples$cluster$s.strata)

mean(mySamples$cluster$means)
mean(sscsample.data$income)
#9.11658-8.994273=0.122307

#e
sapply(mySamples,function(x)sd(x$means))

sapply(mySamples,function(x)IQR(x$means))

#2.2
#a
xdesign(corr = 0.8)




