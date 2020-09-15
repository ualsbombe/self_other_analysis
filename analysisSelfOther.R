rm(list=ls())
graphics.off()
setwd('~/Dropbox/JoenssonLouM.Fl/dataToAnalyse/')
library(matrixStats)
library(MASS)

dt <- read.csv('merged_results_new.csv',sep=';')
min <- read.csv('mincontrast.csv',sep=';')

self <- dt[which(dt$cond=='self'),]
other <- dt[which(dt$cond=='other'),]

diff <- self$contrast_self - other$contrast_other

ns = length(unique(dt$id))
maxTrials = max(dt$self_trialcnt,na.rm=TRUE)

means = matrix(nrow=ns,ncol=maxTrials)
mSelf = means
mOther = means
for(i in 1:ns)
{
	maxTrialThisN = max(dt$self_trialcnt[which(dt$id==i)],na.rm=TRUE)
	
	means[i,1:maxTrialThisN] = diff[(min(which(self$id==i))):(max(which(self$id==i)))]
	
	mSelf[i,1:maxTrialThisN] = self$contrast_self[self$id==i]
	
	mOther[i,1:maxTrialThisN] = other$contrast_other[other$id==i]
}

cMeans <- colMeans(means,na.rm=TRUE)
cMSelf <- colMeans(mSelf,na.rm=TRUE)
cMOther <- colMeans(mOther,na.rm=TRUE)

quartz()
plot(cMSelf,type='l',ylim=c(-15,60),col='red')
lines(cMOther,col='blue')
lines(cMeans,col='green')
lines(1:length(cMeans),rep(0,length(cMeans)))

### divide into self and other

## do self first
nsDopa = unique(dt$id[which(dt$givenDopa==1)])
maxTrials = max(dt$self_trialcnt,na.rm=TRUE)

meansDopa = matrix(nrow=length(nsDopa),ncol=maxTrials)
mDopaSelf = meansDopa
mDopaOther = meansDopa

j = 1
for(i in nsDopa)
{
	maxTrialThisN = max(dt$self_trialcnt[which(dt$id==i)],na.rm=TRUE)
	
	mDopaSelf[j,1:maxTrialThisN] = self$contrast_self[self$id==i]
	
	mDopaOther[j,1:maxTrialThisN] = other$contrast_other[other$id==i]
	j = j + 1
}

## then do other (ugly as hell, this script ain't it?)
nsPlac = unique(dt$id[which(dt$givenDopa==0)])
maxTrials = max(dt$self_trialcnt,na.rm=TRUE)

meansPlac = matrix(nrow=length(nsPlac),ncol=maxTrials)
mPlacSelf = meansPlac
mPlacOther = meansPlac
j = 1
for(i in nsPlac)
{
	maxTrialThisN = max(dt$self_trialcnt[which(dt$id==i)],na.rm=TRUE)
	
	mPlacSelf[j,1:maxTrialThisN] = self$contrast_self[self$id==i]
	
	mPlacOther[j,1:maxTrialThisN] = other$contrast_other[other$id==i]
	j = j + 1
}

trialsOfInterest = 28:75
cMDopaSelf <- colMeans(mDopaSelf,na.rm=TRUE)
cMDopaOther <- colMeans(mDopaOther,na.rm=TRUE)
dopaSelfCONF <- colSds(mDopaSelf,na.rm=TRUE)[trialsOfInterest]/sqrt(36)*qt(0.975,36)
dopaOtherCONF <- colSds(mDopaOther,na.rm=TRUE)[trialsOfInterest]/sqrt(38)*qt(0.975,38)
cMDopaSelf <- cMDopaSelf[trialsOfInterest]
cMDopaOther <- cMDopaOther[trialsOfInterest]

cMPlacSelf <- colMeans(mPlacSelf,na.rm=TRUE)
cMPlacOther <- colMeans(mPlacOther,na.rm=TRUE)
placSelfCONF <- colSds(mPlacSelf,na.rm=TRUE)[trialsOfInterest]/sqrt(36)*qt(0.975,36)
placOtherCONF <- colSds(mPlacOther,na.rm=TRUE)[trialsOfInterest]/sqrt(38)*qt(0.975,38)





cMPlacSelf <- cMPlacSelf[trialsOfInterest]
cMPlacOther <- cMPlacOther[trialsOfInterest]

diff1 <- cMDopaSelf - cMPlacSelf
diff2 <- cMDopaOther - cMPlacOther
x=28:75
quartz()
par(mfrow=c(1,2))
plot(x=x,cMDopaSelf,type='l',ylim=c(-15,60),col='red',lwd=3)
lines(x=x,cMPlacSelf,col='blue',lwd=3)
lines(x=x,diff1,col='green')
lines(x=x,cMDopaSelf+dopaSelfCONF,lty=2,col='red')
lines(x=x,cMDopaSelf-dopaSelfCONF,lty=2,col='red')
lines(x=x,cMPlacSelf+placSelfCONF,lty=2,col='blue')
lines(x=x,cMPlacSelf-placSelfCONF,lty=2,col='blue')
legend('topright',c('dopaSelf','placSelf','diff'),lwd=1,col=c('red','blue','green'))
lines(x=x,rep(0,length(cMPlacSelf)),lty=2)
plot(x=x,cMDopaOther,type='l',ylim=c(-15,60),col='red',lwd=3)
lines(x=x,cMPlacOther,col='blue',lwd=3)
lines(x=x,diff2,col='green')
lines(x=x,cMDopaOther+dopaOtherCONF,lty=2,col='red')
lines(x=x,cMDopaOther-dopaOtherCONF,lty=2,col='red')
lines(x=x,cMPlacOther+placOtherCONF,lty=2,col='blue')
lines(x=x,cMPlacOther-placOtherCONF,lty=2,col='blue')
legend('topright',c('dopaOther','placOther','diff'),lwd=1,col=c('red','blue','green'))
lines(x=x,rep(0,length(cMPlacSelf)),lty=2)

## we'll start from trial 28, this includes 35/74 ≈ 50 % of the trials it's found by : "sum(table(1dbStart)[1:10])"
db1StartSelf <- tapply(dt$self_trialcnt[dt$col_stepsize_self==1],dt$id[dt$col_stepsize_self==1],min,na.rm=TRUE)
db1StartOther <- tapply(dt$other_trialcnt[dt$col_stepsize_other==1],dt$id[dt$col_stepsize_other==1],min,na.rm=TRUE)
quartz()
par(mfrow=c(1,2))
hist(db1StartSelf)
hist(db1StartOther)
trialStart = 28
trialEnd = 75 ## got this from Morten, ~ 50 % of participants have this many trials


## get ready for modelling
# self <- self[which(self$self_trialcnt<=max(trialsOfInterest)),]
# self <- self[self$self_trialcnt>=trialStart,]
# other <- other[which(other$other_trialcnt<=trialEnd),]
# other <- other[other$other_trialcnt>=trialStart,]

dt <- rbind(self,other)
dt$contrast <- c(dt$contrast_self[!is.na(dt$contrast_self)],dt$contrast_other[!is.na(dt$contrast_other)])
dt$correct <- c(dt$resp_dis_self[!is.na(dt$resp_dis_self)],dt$resp_dis_other[!is.na(dt$resp_dis_other)])
dt$pas <- c(dt$resp_awa_self[!is.na(dt$resp_awa_self)],dt$resp_awa_other[!is.na(dt$resp_awa_other)])

## factorise
dt$givenDopa <- factor(dt$givenDopa)
dt$correct <- factor(dt$correct)

m.logis <- glm(correct ~ contrast + givenDopa*cond, data=dt,family=binomial)

## create data set for linear modelling

meansDopaSelf <- as.vector(tapply(dt$contrast[dt$givenDopa==1&dt$cond=='self'],dt$id[dt$givenDopa==1&dt$cond=='self'],mean))
meansDopaOther <- as.vector(tapply(dt$contrast[dt$givenDopa==1&dt$cond=='other'],dt$id[dt$givenDopa==1&dt$cond=='other'],mean))

meansPlacSelf <- as.vector(tapply(dt$contrast[dt$givenDopa==0&dt$cond=='self'],dt$id[dt$givenDopa==0&dt$cond=='self'],mean))
meansPlacOther <- as.vector(tapply(dt$contrast[dt$givenDopa==0&dt$cond=='other'],dt$id[dt$givenDopa==0&dt$cond=='other'],mean))

relTimeDopaSelf <- as.vector(tapply(dt$reltime[dt$givenDopa==1&dt$cond=='self'],dt$id[dt$givenDopa==1&dt$cond=='self'],mean))

relTimeDopaOther <- as.vector(tapply(dt$reltime[dt$givenDopa==1&dt$cond=='other'],dt$id[dt$givenDopa==1&dt$cond=='other'],mean))

relTimePlacSelf <- as.vector(tapply(dt$reltime[dt$givenDopa==0&dt$cond=='self'],dt$id[dt$givenDopa==0&dt$cond=='self'],mean))

relTimePlacOther <- as.vector(tapply(dt$reltime[dt$givenDopa==0&dt$cond=='other'],dt$id[dt$givenDopa==0&dt$cond=='other'],mean))



firstContrastsSelfDopa <- dt$contrast_self[which(dt$self_trialcnt==28&dt$givenDopa==1)]

firstContrastsOtherDopa <- dt$contrast_other[which(dt$other_trialcnt==28&dt$givenDopa==1)]

firstContrastsSelfPlac <- dt$contrast_self[which(dt$self_trialcnt==28&dt$givenDopa==0)]

firstContrastsOtherPlac <- dt$contrast_other[which(dt$other_trialcnt==28&dt$givenDopa==0)]


anovaDt <- data.frame('meanContrasts'=c(meansDopaSelf,meansDopaOther,meansPlacSelf,meansPlacOther),'givenDopa'=factor(c(rep(1,2*length(meansDopaSelf)),rep(0,2*length(meansPlacSelf)))),'cond'=factor(c(rep('self',length(meansDopaSelf)),rep('other',length(meansDopaOther)),rep('self',length(meansPlacSelf)),rep('other',length(meansPlacOther)))),'startingContrast'=c(firstContrastsSelfDopa,firstContrastsOtherDopa,firstContrastsSelfPlac,firstContrastsOtherPlac),'id'=factor(c(rep(unique(dt$id[which(dt$givenDopa==1)]),2),rep(unique(dt$id[which(dt$givenDopa==0)]),2))),'relTime' = c(relTimeDopaSelf,relTimeDopaOther,relTimePlacSelf,relTimePlacOther))

anovaDt$minContrast <- c(min$min_contrast[1:36],min$min_contrast[75:110],min$min_contrast[37:74],min$min_contrast[111:148]) ## due to differing orderings between me and mOrten


m.full <- lm(meanContrasts ~ givenDopa * cond + startingContrast + relTime,data = anovaDt)
m.add <- lm(meanContrasts ~ givenDopa + cond,data=anovaDt)

#m.error <- aov(meanContrasts ~ givenDopa*cond + Error(id/startingContrast),data=anovaDt)

library(lme4)

mm.full <- lmer(meanContrasts ~ givenDopa*cond + relTime +  + (1|id) + (1|startingContrast) ,data=anovaDt,REML=FALSE)
#mm.add <- lmer(meanContrasts ~ givenDopa+cond + (id|startingContrast)+ (1|id),data=anovaDt,REML=FALSE)

## carry on on model where all trials are included (dt)

dt$trialN <- c(dt$self_trialcnt[which(dt$cond=='self')],dt$other_trialcnt[which(dt$cond=='other')])

dt$trialN <- factor(dt$trialN)
dt$id <- factor(dt$id)

m.rep.cov.full <- lm(log(contrast) ~ trialN+cond*givenDopa + reltime + likeability,data=dt)

m.rep.cov.add <- lm(log(contrast) ~ trialN+cond+givenDopa + reltime + likeability,data=dt)


mm.rep.cov <- lmer(contrast ~ trialN+cond*givenDopa + reltime + (1|id) ,data=dt)

mm.rep.cov.noInt <- lmer(contrast ~ trialN+cond+givenDopa + reltime + (1|id) ,data=dt)

### do on minContrast


m.full <- lm(minContrast ~ givenDopa * cond,data = anovaDt)
mm.full <- lmer(minContrast ~ givenDopa*cond + (1|startingContrast:id) ,data=anovaDt,REML=FALSE)

means <- as.vector(tapply(anovaDt$minContrast,interaction(anovaDt$givenDopa,anovaDt$cond),mean))
sds <- as.vector(tapply(anovaDt$minContrast,interaction(anovaDt$givenDopa,anovaDt$cond),sd))
ns <- as.vector(tapply(anovaDt$minContrast,interaction(anovaDt$givenDopa,anovaDt$cond),length))

conf <- sds/sqrt(ns)*qt(0.975,ns)
col = c('red','blue','green','purple')
plot(means,ylim=c(20,40),col=col)
#errbar(x=1:4,y=means,yplus=means+conf,yminus=means-conf,col=col)

mm.trial <- lmer(contrast ~ cond * givenDopa + (1|id) + (1|trialN),data=dt,REML=FALSE)
mm.trial.noInt <- lmer(contrast ~ cond + givenDopa + (1|id) + (1|trialN),data=dt,REML=FALSE)

m.pas <- lm(pas ~  cond * givenDopa + likeability,data=dt)
dt$pas <- factor(dt$pas)
m.olr  <- polr(pas ~ cond * givenDopa,data=dt,Hess=TRUE)
m.olr2 <- polr(pas ~ cond + givenDopa,data=dt,Hess=TRUE)

quartz()
dt$contrast <- log(dt$contrast)
par(mfrow=c(1,2))
h <- hist(dt$contrast,xlim=c(1,6),main='hist of log(contrast)')
xfit <- seq(min(dt$contrast),max(dt$contrast),length=40)
		yfit <- dnorm(xfit,mean=mean(dt$contrast),sd=sd(dt$contrast))
		yfit <- yfit*diff(h$mids[1:2])*length(dt$contrast)
		lines(xfit,yfit,lwd=2)

qqnorm(dt$contrast)
qqline(dt$contrast)


dt$contrastNew <- c(dt$contrast[2:length(dt$contrast)],NA)
indexForRemoval <- length(dt$contrast[dt$id==1&dt$cond=='other'])

a <- seq(48,7000,48)
b <- 1:dim(dt)[1]
c <- setdiff(b,a)

dtNew <- dt[c,]

dtNewid1 <- subset(dtNew,dtNew$id==1)