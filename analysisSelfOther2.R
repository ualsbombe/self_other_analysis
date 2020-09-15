selfA <- tapply(dt$resp_awa_self,interaction(dt$id,dt$givenDopa),mean,na.rm=T)

selfA <- as.vector(selfA)[!is.na(selfA)]

otherA <- tapply(dt$resp_awa_other,interaction(dt$id,dt$givenDopa),mean,na.rm=T)
otherA <- as.vector(otherA)[!is.na(otherA)]

len <- length(selfA)
len2 <- length(otherA)

diffSelf <- selfA[(len/2+1):len] - selfA[1:(len/2)]
diffOther <- otherA[(len/2+1):len] - otherA[1:(len/2)]

data <- data.frame('awa'=c(selfA,otherA),'dopa'=c(rep(c(0,1),each=len/2),rep(c(0,1),each=len2/2)),'cond'=c(rep('self',len),rep('other',len2)))

m <- lm(awa ~ dopa*cond,data=data)
