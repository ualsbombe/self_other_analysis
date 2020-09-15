## clear and load data

rm(list=ls())
graphics.off()
setwd('~/Dropbox/JoenssonLouM.Fl/dataToAnalyse/')
dt <- read.csv('merged_results_new.csv',sep=';')

## get dt ordered nicely and without NAs
self <- dt[which(dt$cond=='self'),]
other <- dt[which(dt$cond=='other'),]
dt <- rbind(self,other)
dt$contrast <- c(dt$contrast_self[!is.na(dt$contrast_self)],dt$contrast_other[!is.na(dt$contrast_other)])
dt$correct <- factor(c(dt$resp_dis_self[!is.na(dt$resp_dis_self)],dt$resp_dis_other[!is.na(dt$resp_dis_other)]))
dt$pas <- factor(c(dt$resp_awa_self[!is.na(dt$resp_awa_self)],dt$resp_awa_other[!is.na(dt$resp_awa_other)]))

## we only analyse db1Trials
db1Trials <- dt$col_stepsize_self==1|dt$col_stepsize_other==1
dt <- dt[which(db1Trials),] ## only analyse db1trials
dt$givenDopa <- factor(dt$givenDopa)


## do Ordinal regression
print('ORDINAL REGRESSION')
#load packages
library(ordinal)

## ordinal regression (mixed models)
mm.olr  <- clmm(pas ~ cond * givenDopa + likeability +  (1|id),data=dt,control=list()) ## takes some time

print(s.olr <- summary(mm.olr))


print('LOGISTIC REGRESSION')
library(lme4)
## dichotomous PAS
dt$dichPas <- rep(1,length(dt$pas))
for(i in 1:length(dt$dichPas))
{
	if(as.numeric(dt$pas[i])>1) dt$dichPas[i] <- 2
}
dt$dichPas <- factor(dt$dichPas)

mm.logit <- glmer(dichPas ~ cond * givenDopa + likeability + (1|id),data=dt,family=binomial())

print(s.logit <- summary(mm.logit))

## do Markov modelling
print('MARKOV CHAINS')
library(markovchain)

## long loop for setting up steps
#graphics.off()
markovMeans <- numeric(2*length(unique(dt$id)))
normalMeans <- markovMeans
startingContrast <- markovMeans
absorbed <- markovMeans
k=0
quartz()
par(mfrow=c(2,4))
for(j in c('self','other'))
{
	if(j=='other'){k=length(markovMeans)/2}
	for(i in as.numeric(unique(dt$id)))
	{
		if(j=='self') db1Trial <- min(which(dt$col_stepsize_self[!is.na(dt$col_stepsize_self)&dt$id==i]==1))
		if(j=='other') db1Trial <- min(which(dt$col_stepsize_other[!is.na(dt$col_stepsize_other)&dt$id==i]==1))
		
		self1 <- subset(dt,dt$id==i & dt$cond==j)
		
		#self1$contrast <- exp(self1$contrast)
		#self1 <- self1[db1Trial:length(self1$contrast),]
		
		sequence <- self1$contrast
		#sequence <- dt$contrast
		sequenceMatr<-createSequenceMatrix(sequence,sanitize=FALSE)
		mcFitMLE<-markovchainFit(data=sequence)
		#mcFitBSP<-markovchainFit(data=sequence,method="bootstrap",nboot=5, name="Bootstrap Mc")
		
		ss <- steadyStates(mcFitMLE[[1]])
		#print(ss)
		as <- absorbingStates(mcFitMLE[[1]])
		if(length(as) == 0) as <- NA
		
		markovMeans[i+k] <- sum(as.numeric(colnames(ss))*ss)
		normalMeans[i+k] <- mean(self1$contrast)
		absorbed[i+k] <- as
		startingContrast[i+k] <- self1$contrast[1]
	  	if(!is.na(absorbed[i+k])) plot(self1$contrast)
#	  	print(sum(which(!is.na(absorbed)))==1)
#	  	print(which(!is.na(absorbed)==i))
	}
}

quartz()
diff <- markovMeans-normalMeans
plot(diff)
#graphics.off()
dopaSubjects <- as.numeric(unique(dt$id[dt$givenDopa==1]))
markovDt <- data.frame('id'=1:k,'contrast'=markovMeans,'cond'=c(rep('self',k),rep('other',k)),'givenDopa'=0,'startContrast'=startingContrast)
markovDt$givenDopa[c(dopaSubjects,dopaSubjects+k)] <- 1

temp <- which(!is.na(absorbed))
for(i in 1:length(temp))
{
	if(temp[i] > length(unique(dt$id))) temp[i] <- temp[i] -74
} 

subsToBeKept <- setdiff(as.numeric(row.names(markovDt)),temp)
markovDt <- markovDt[subsToBeKept,]
#markovDt <- markovDt[is.na(absorbed),]


m <- lm(contrast ~ cond*givenDopa,data=markovDt)
markovDt$id <- factor(markovDt$id)
markovDt$givenDopa <- factor(markovDt$givenDopa)
mm.full <- lmer(contrast ~ cond *givenDopa + (1|id),data=markovDt,REML=FALSE)
mm.add <- lmer(contrast ~ cond +givenDopa  + (1|id),data=markovDt,REML=FALSE)
an <- anova(mm.add,mm.full)

print(s.full <- summary(mm.full))
print(an)