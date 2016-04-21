library(ape)
library(geiger)
library(laser)
library(TreeSim)
library(TreePar)

#these three lines open the pomp tree and calculate bd parameters for simulations, etc
pomptree<-read.tree("pompilinaechronstrictoNEW2014.tre.txt")#reads latest pompilinae tree with all terminals
btimes<-branching.times(pomptree)#calculates empirical speciation times of Pompilinae tree
treebd <- make.bd(pomptree, sampling.f = 1.0) #these three lines fit a bd model and calculates mu and lambda


#the next steps are same as previous but instead of TreeSim using CorSim and divides tree into 2 clades to simulate missing branching times at appropriate ages
#this loop simulates missing branching times up to the total exoected for 2 separate clades of Pompilinae 1000 times and creates a single data frame with all the results called CorsimList 

clade1<-read.tree("clade1pompilinae.tre")#reads clade of Agenioideus+poecilo etc
btimesclade1<-branching.times(clade1)#calculates empirical speciation times of clade 1 tree

CorsimListclade1 = matrix(NA, nrow=1000, ncol=232)
for(i in 1:1000) {
  corsimclade1<-corsim(x = btimesclade1, lambda=9.111433e-02, mu=3.242355e-07, missing=214,told=20,tyoung=0)#simulates missing speciation events in the incomplete phylogeny
  corsimclade1=as.data.frame(corsimclade1)
  CorsimListclade1[i,] = corsimclade1[,1]
}

clade2<-read.tree("clade2pompilinae.tre")#reads clade of Agenioideus+poecilo etc
btimesclade2<-branching.times(clade2)#calculates empirical speciation times of clade 1 tree

CorsimListclade2 = matrix(NA, nrow=1000, ncol=648)
for(i in 1:1000) {
  corsimclade2<-corsim(x = btimesclade2, lambda=9.111433e-02, mu=3.242355e-07, missing=598,told=10,tyoung=0)#simulates missing speciation events in the incomplete phylogeny
  corsimclade2=as.data.frame(corsimclade2)
  CorsimListclade2[i,] = corsimclade2[,1]
}
#next lines can resample from corsimlist. Not necessary because it is not random so it belongs to dataset
#transposedcorsim=t(CorsimList)#transposes matrix
#prunedcorsim = transposedcorsim[sample(ncol(transposedcorsim), size=71, replace = FALSE), ]  #prunes branching times from corsimlist
#finalcorsimpruned=t(prunedcorsim)#transposes the pruned matrix and returns data to use

#the following lines run fitdAICrc on all corsimclade1 and corsimclade2 trees manually because the batch function is not working (see above)
#this loop returns a matrix with columns with AIC values for the different parameters

fitdAICrcmancorsimclade1 = matrix(NA, nrow=1000, ncol=6)
fitdAICrcdeltacorsimclade1 = matrix(NA, nrow=1000, ncol=6)
colnames(fitdAICrcmancorsimclade1) = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate")
for(i in 1:1000) {
  fitAICcorsimclade1 <-fitdAICrc(CorsimListclade1[i,], modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate"), ints = 500)
  fitdAICrcmancorsimclade1[i,] = fitAICcorsimclade1$AIC
  fitdAICrcdeltacorsimclade1[i,] = fitAICcorsimclade1$dAIC 
}

mean(fitdAICrcmancorsimclade1[,6]) #this calculates the mean AIC for each model
sd(fitdAICrcmancorsimclade1[,6]) #this calculates the standard deviation of AIC for each model
hist(fitdAICrcmancorsimclade1[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix

mean(fitdAICrcdeltacorsimclade1[,6]) #this calculates mean of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
sd(fitdAICrcdeltacorsimclade1[,6]) #this calculares sd of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
hist(fitdAICrcdeltacorsimclade1[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel

fitdAICrcmancorsimclade2 = matrix(NA, nrow=1000, ncol=6)
fitdAICrcdeltacorsimclade2 = matrix(NA, nrow=1000, ncol=6)
colnames(fitdAICrcmancorsimclade2) = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate")
for(i in 1:1000) {
  fitAICcorsimclade2 <-fitdAICrc(CorsimListclade2[i,], modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate"), ints = 500)
  fitdAICrcmancorsimclade2[i,] = fitAICcorsimclade2$AIC
  fitdAICrcdeltacorsimclade2[i,] = fitAICcorsimclade2$dAIC 
}

mean(fitdAICrcmancorsimclade2[,6]) #this calculates the mean AIC for each model
sd(fitdAICrcmancorsimclade2[,6]) #this calculares sd of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
hist(fitdAICrcmancorsimclade2[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix

mean(fitdAICrcdeltacorsimclade2[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
sd(fitdAICrcdeltacorsimclade2[,6]) #this calculares sd of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
hist(fitdAICrcdeltacorsimclade2[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel


#the following lines find the distribution of r1 r2 r3 and st1 st2 for yule3rate which was the best model in both cases

parametersclade1 = matrix(NA, nrow=1000, ncol=7)
colnames(parametersclade1) = c("LH","r1", "r2", "r3", "st1", "st2","AIC")
for(i in 1:1000) {
  parametersyule3rate <-yule3rate(CorsimListclade1[i,], ints = 500)
  parametersclade1[i,] = parametersyule3rate
}

parametersclade2 = matrix(NA, nrow=1000, ncol=7)
colnames(parametersclade2) = c("LH","r1", "r2", "r3", "st1", "st2","AIC")
for(i in 1:1000) {
  parametersyule3rate2 <-yule3rate(CorsimListclade2[i,], ints = 500)
  parametersclade2[i,] = parametersyule3rate2
}


mean(parametersclade1[,6]) #this calculates the mean for any column
hist(parametersclade1[,2]) #this creates a histogram of column n(change to calculate for different values) from the matrix

mean(parametersclade2[,6]) #this calculates the mean for any column
hist(parametersclade2[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix

#the following lines simulate random constant bd trees with the total taxa from clade 1 and clade 2 separately and run dAIC analyses on these as null distribution

treesimclade1<-sim.bd.taxa.age(n=232, numbsim=1000, lambda=9.111433e-02, mu=3.242355e-07, frac = 1, age=31, mrca = FALSE)#simulates 1000 trees with diversification rates calculated by bd(btimes)

treesimclade2<-sim.bd.taxa.age(n=648, numbsim=1000, lambda=9.111433e-02, mu=3.242355e-07, frac = 1, age=31, mrca = FALSE)#simulates 1000 trees with diversification rates calculated by bd(btimes)

for(i in 1:length(treesimclade1)){
  write.tree(treesimclade1[[i]],file="clade1simtree.tre", append=TRUE)
}#writes simulated trees to a file
Btimessimtreesclade1=getBtimes.batch(fname="clade1simtree.tre")#creates a vector of branching times for all simulated trees of clade 1

for(i in 1:length(treesimclade2)){
  write.tree(treesimclade2[[i]],file="clade2simtree.tre", append=TRUE)
}#writes simulated trees to a file
Btimessimtreesclade2=getBtimes.batch(fname="clade2simtree.tre")#creates a vector of branching times for all simulated trees of clade 2

fitdAICrcmantreesimclade1 = matrix(NA, nrow=1000, ncol=6)
fitdAICrcdeltatreesimclade1 = matrix(NA, nrow=1000, ncol=6)
colnames(fitdAICrcmantreesimclade1) = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate")
for(i in 1:1000) {
  fitAICtreesimclade1 <-fitdAICrc(Btimessimtreesclade1[i,], modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate"), ints = 500)
  fitdAICrcmantreesimclade1[i,] = fitAICtreesimclade1$AIC
  fitdAICrcdeltatreesimclade1[i,] = fitAICtreesimclade1$dAIC 
}

mean(fitdAICrcmantreesimclade1[,6]) #this calculates the mean AIC for each model
sd(fitdAICrcmantreesimclade1[,6]) #this calculates the standard deviation of AIC for each model
hist(fitdAICrcmantreesimclade1[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix

mean(fitdAICrcdeltatreesimclade1[,6]) #this calculates mean of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
sd(fitdAICrcdeltatreesimclade1[,6]) #this calculares sd of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
hist(fitdAICrcdeltatreesimclade1[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel

fitdAICrcmantreesimclade2 = matrix(NA, nrow=1000, ncol=6)
fitdAICrcdeltatreesimclade2 = matrix(NA, nrow=1000, ncol=6)
colnames(fitdAICrcmantreesimclade2) = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate")
for(i in 1:1000) {
  fitAICtreesimclade2 <-fitdAICrc(Btimessimtreesclade2[i,], modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate","yule3rate"), ints = 500)
  fitdAICrcmantreesimclade2[i,] = fitAICtreesimclade2$AIC
  fitdAICrcdeltatreesimclade2[i,] = fitAICtreesimclade2$dAIC 
}

mean(fitdAICrcmantreesimclade2[,6]) #this calculates the mean AIC for each model
sd(fitdAICrcmantreesimclade2[,6]) #this calculates the standard deviation of AIC for each model
hist(fitdAICrcmantreesimclade2[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix

mean(fitdAICrcdeltatreesimclade2[,6]) #this calculates mean of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
sd(fitdAICrcdeltatreesimclade2[,6]) #this calculares sd of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel
hist(fitdAICrcdeltatreesimclade2[,6]) #this creates a histogram of column n(change to calculate for different values) from the matrix in this case is the deltaAIC of the yule2ratemodel

#the next lines run a t test on the dAICrc of corsim and treesim trees

x=fitdAICrcdeltatreesimclade1[,1]
y=fitdAICrcdeltacorsimclade1[,1]
t.test(x,y,mu=0,conf.level=0.95)
