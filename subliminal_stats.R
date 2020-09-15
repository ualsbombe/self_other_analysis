## clear and load data

rm(list=ls())
graphics.off()
#setwd('C:/Users/morten/Dropbox/PhD/projects/subliminal/experiment/')
setwd('/Users/Lau/Dropbox/JoenssonLouM.Fl/dataToAnalyse/')
dt <- read.csv('merged_results_new.csv',sep=';')
logit.inv = function(x){exp(x)/(1+exp(x))} ## logit inverse


## set functions
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

## get dt ordered nicely and without NAs
self <- dt[which(dt$cond=='self'),]
other <- dt[which(dt$cond=='other'),]
dt <- rbind(self,other)
dt$contrast <- c(dt$contrast_self[!is.na(dt$contrast_self)],dt$contrast_other[!is.na(dt$contrast_other)])
dt$correct <- factor(c(dt$resp_dis_self[!is.na(dt$resp_dis_self)],dt$resp_dis_other[!is.na(dt$resp_dis_other)]))
dt$pas <- factor(c(dt$resp_awa_self[!is.na(dt$resp_awa_self)],dt$resp_awa_other[!is.na(dt$resp_awa_other)]))
dt$correct_val = as.numeric(as.character(dt$correct));
dt$cond = factor(dt$cond, levels = c("self","other"), ordered=FALSE) # force order

## we only analyse db1Trials
db1Trials <- dt$col_stepsize_self==1|dt$col_stepsize_other==1
dt <- dt[which(db1Trials),] ## only analyse db1trials
dt$givenDopa <- factor(dt$givenDopa)

## extract information about accuracies
subjmean  = aggregate(correct_val ~ givenDopa * cond * id, data=dt, FUN=mean)
totalmean = aggregate(correct_val ~ givenDopa * cond, data=subjmean, FUN=mean)
totalsd   = aggregate(correct_val ~ givenDopa * cond, data=subjmean, FUN=sd)

# cond-effect? (paired thus ranksum)
x = subset(subjmean$correct_val, subjmean$cond == "self" );
y = subset(subjmean$correct_val, subjmean$cond == "other" );
wilcox.test(x,y, paired=TRUE)
#
# data:  x and y
# W = 3033, p-value = 0.2587

# group effect?
x = subset(subjmean$correct_val, subjmean$givenDopa == 1 );
y = subset(subjmean$correct_val, subjmean$givenDopa == 0 );
wilcox.test(x,y, paired=FALSE)
# data:  x and y
# W = 3197, p-value = 0.07724
# alternative hypothesis: true location shift is not equal to 0

subjstats = aov(correct_val~cond*givenDopa, data=subjmean)
summary(subjstats)

## create barplot of contrast based on subject average MJ
subjmean  = aggregate(contrast ~ givenDopa * cond * id, data=dt, FUN=mean)
totalmean = aggregate(contrast ~ givenDopa * cond, data=subjmean, FUN=mean)
totalsd   = aggregate(contrast ~ givenDopa * cond, data=subjmean, FUN=sd)
totallen  = aggregate(contrast ~ givenDopa * cond, data=subjmean, FUN=length)
se = totalsd$contrast / sqrt(totallen$contrast)

y.means = totalmean$contrast

ymax <- max(totalmean$contrast+se)
barx = barplot(y.means, ylim=c(0,ymax), names.arg=c('self placebo', 'self Sinemet','other placebo','other Sinemet'), main="Contrast vs. group x condition", ylab="Contrast")
error.bar(barx, y.means, se)

## do stats on the subject level

subjstats = aov(contrast~cond*givenDopa, data=subjmean)
summary(subjstats)
#                Df Sum Sq Mean Sq F value Pr(>F)
#cond             1      2    1.89   0.009  0.924
#givenDopa        1    148  148.39   0.716  0.399
#cond:givenDopa   1     32   31.96   0.154  0.695
#Residuals      144  29824  207.11 

# subject level analysis of awareness
# split pas in two groups low conf. (1-2) and high conf. (3-4) 
dt$dichPas_mj <- rep(0,length(dt$pas))
for(i in 1:length(dt$dichPas_mj))
{
  if(as.numeric(dt$pas[i])>1) dt$dichPas_mj[i] <- 1  #1 vs 2,3,4
}

subjmean  = aggregate(dichPas_mj ~ givenDopa * cond * id, data=dt, FUN=mean)
totalmean = aggregate(dichPas_mj ~ givenDopa * cond, data=subjmean, FUN=mean)
totalsd   = aggregate(dichPas_mj ~ givenDopa * cond, data=subjmean, FUN=sd)

kruskal.test(dichPas_mj~cond, data=subjmean)
kruskal.test(dichPas_mj~givenDopa, data=subjmean)

## do Ordinal regression (awareness)
print('ORDINAL REGRESSION')
#load packages
library(ordinal)

## ordinal regression (mixed models)
#mm.olr  <- clmm(pas ~ cond * givenDopa + likeability +  (1|id),data=dt,control=list()) ## takes some time
mm.olr  <- clmm(pas ~ cond * givenDopa + (1|id),data=dt,control=list()) ## takes some time


print(s.olr <- summary(mm.olr))

print('LOGISTIC REGRESSION')
library(lme4)
## dichotomous PAS
dt$dichPas <- rep(1,length(dt$pas))
for(i in 1:length(dt$dichPas))
{
	if(as.numeric(dt$pas[i])>1) dt$dichPas[i] <- 2  #1 vs 2,3,4
}
dt$dichPas <- factor(dt$dichPas)

##mm.logit <- glmer(dichPas ~ cond * givenDopa + likeability + (1|id),data=dt,family=binomial())
mm.logit <- glmer(dichPas ~ cond * givenDopa + (1|id),data=dt,family=binomial())

print(s.logit <- summary(mm.logit))

## MJ metacognition (logistic regression)
dt$correct_val = as.numeric(dt$correct)-1

subjmean  = aggregate(correct_val ~ dichPas*givenDopa*cond*id, data=dt, FUN=mean)
totalmean = aggregate(correct_val ~ dichPas*givenDopa*cond, data=subjmean, FUN=mean)


dt$correct_val <- factor(dt$correct_val)
dt$cond <- factor(dt$cond,ordered=FALSE)
metacog = glmer(correct_val ~ dichPas*givenDopa + dichPas*cond + (1|id), data=dt, family=binomial())    
print(s.metacog <- summary(metacog))
metacogFull <- update(metacog,~.+dichPas:cond:givenDopa,data=dt)
metacog.noInt <- update(metacog,~.-dichPas:givenDopa,data=dt)
print(an.metacog <- anova(metacog.noInt,metacog,metacogFull))

coef <- s.metacog$coefficients
selfDopaPas2 <- logit.inv(coef[1,1]+coef[2,1]+coef[3,1]+coef[5,1])
selfDopaPas1 <- logit.inv(coef[1,1]+coef[3,1])

otherDopaPas2 <- logit.inv(coef[1,1]+coef[2,1]+coef[3,1]+coef[4,1]+coef[5,1]+coef[6,1])
otherDopaPas1 <- logit.inv(coef[1,1]+coef[3,1]+coef[4,1])

selfPlacPas2 <- logit.inv(coef[1,1]+coef[2,1])
selfPlacPas1 <- logit.inv(coef[1,1])

otherPlacPas2 <- logit.inv(coef[1,1]+coef[2,1]+coef[4,1]+coef[6,1])
otherPlacPas1 <- logit.inv(coef[1,1]+coef[4,1])



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

markovDt$id <- factor(markovDt$id)
markovDt$givenDopa <- factor(markovDt$givenDopa)
mm.full <- lmer(contrast ~ cond *givenDopa + (1|id),data=markovDt,REML=FALSE)
mm.add <- lmer(contrast ~ cond +givenDopa  + (1|id),data=markovDt,REML=FALSE)
an <- anova(mm.add,mm.full)

print(s.full <- summary(mm.full))
print(an)





