library('depmixS4')
data("speed")
set.seed(1)
modNull <- depmix(response = rt ~ 1 ,data=speed,nstates =1, trstart=runif(1))

mod <- depmix(response = rt ~ 1,data=speed,nstates =2, trstart=runif(4))

modPacc <- depmix(response = rt ~ 1 + scale(Pacc),data=speed,nstates =2, trstart=runif(4))

fm <- fit(mod)

mod.cov <- depmix(rt ~ 1, data = speed, nstates = 2, family=gaussian(), transition = ~ scale(Pacc), instart=runif(2))
fm.cov <- fit(mod.cov,verbose=FALSE)

summary(fm.cov, which = 'transition')

mod.multi <- depmix(list(rt ~ 1,corr ~1),data= speed,nstates=2, family=list(gaussian(),multinomial('identity')),transition= ~scale(Pacc),instart=runif(2))

fm.multi <- fit(mod.multi)

summary(fm.multi)