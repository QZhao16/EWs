library(survival)

atelopus = read.csv("atelopus.csv") 
attach(atelopus)

m1 = coxph(Surv(time1,time2,censor)~range.size*culturemortprob+range.size*log(AVMD)+range.size*tempchange*X40yr.meantemp
                      +meantemp*range.size+logaltitude+precip+cluster(speciesID)) 
m1
