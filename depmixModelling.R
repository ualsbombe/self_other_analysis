rm(list=ls())
graphics.off()

library(depmixS4)
setwd('~/Dropbox/JoenssonLouM.Fl/dataToAnalyse/')
source('analysisSelfOther.R') ## to get dt
graphics.off()
rm(list=setdiff(ls(),'dt'))


m <- depmix(response = contrast ~ 1, data=dt,nstates=3,trstart=runif(9))

fm <- fit(m)

m.cov <- depmix(contrast ~ 1,data=dt,nstates=2,transition = ~ correct + givenDopa*cond)

fm.cov <- fit(m.cov)

m.cov.add <- depmix(contrast ~ 1,data=dt,nstates=2,transition = ~ correct + givenDopa+cond)

fm.cov.add <- fit(m.cov.add)