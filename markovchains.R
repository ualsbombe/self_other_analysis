## you need to run analysisSelfOther to get "dt"
library(markovchain)
library(lme4)
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
		
		self1$contrast <- exp(self1$contrast)
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


m <- lm(contrast ~ cond * givenDopa ,data=markovDt)
mm.full <- lmer(contrast ~ cond *givenDopa + (1|id),data=markovDt,REML=FALSE)
mm.add <- lmer(contrast ~ cond +givenDopa  + (1|id),data=markovDt,REML=FALSE)
mm.cond <- lmer(contrast ~ cond  + (1|id),data=markovDt,REML=FALSE)
mm.dopa <- lmer(contrast ~ givenDopa  + (1|id),data=markovDt,REML=FALSE)
mm.null<- lmer(contrast ~ 1  + (1|id),data=markovDt,REML=FALSE)



an <- anova(mm.add,mm.full)
quartz()
par(mfrow=c(1,2))
plot(sort(markovDt$contrast[markovDt$cond=='self'&markovDt$givenDopa==1]),type='p',col='red',ylim=c(0,100))
points(sort(markovDt$contrast[markovDt$cond=='self'&markovDt$givenDopa==0]),col='blue')

plot(sort(markovDt$contrast[markovDt$cond=='other'&markovDt$givenDopa==1]),type='p',col='red',ylim=c(0,100))
points(sort(markovDt$contrast[markovDt$cond=='other'&markovDt$givenDopa==0]),col='blue')

markovDt$logContrast <- log(markovDt$contrast)
quartz()
par(mfrow=c(2,2))
xlim=c(2,5)
hist(markovDt$logContrast[markovDt$cond=='self'&markovDt$givenDopa==1],xlim=xlim)
hist(markovDt$logContrast[markovDt$cond=='self'&markovDt$givenDopa==0],xlim=xlim)
hist(markovDt$logContrast[markovDt$cond=='other'&markovDt$givenDopa==1],xlim=xlim)
hist(markovDt$logContrast[markovDt$cond=='other'&markovDt$givenDopa==0],xlim=xlim)

m.log <- lm(logContrast ~ cond * givenDopa+startContrast,data=markovDt)
mm.full.log <- lmer(logContrast ~ cond *givenDopa +startContrast + (1|id),data=markovDt,REML=FALSE)
mm.add.log <- lmer(logContrast ~ cond +givenDopa + startContrast + (1|id),data=markovDt,REML=FALSE)