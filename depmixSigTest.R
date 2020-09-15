#create a 2 state model with one continuous and one binary response
	data(speed)
	mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()))
	# print the model, formulae and parameter values (ie the starting values)
	
	
mod2 <- depmix(list(rt~1,corr~Pacc),data=speed,nstates=2,family=list(gaussian(),multinomial()))

fm <- fit(mod)
fm2 <- fit(mod2)

logdiff <- as.vector(abs(2*logLik(fm2) - 2*logLik(fm)))
df = freepars(fm2) - freepars(fm)
print(dchisq(logdiff,df))