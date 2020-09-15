
nSteps <- character(dim(dt)[1])
n = 1 ## a counter


for(i in unique(dt$id))
{
	for(j in c('self','other'))
	{
		foo <- subset(dt,dt$id==i & dt$cond == j)
		sortedContrasts <- sort(unique(foo$contrast))
		startingContrast <- foo$contrast[1]
		reference <- which(sortedContrasts==startingContrast)
		
		for(k in 1:dim(foo)[1])
		{
			current <- reference - which(sortedContrasts==foo$contrast[k])
			affix <- abs(current)
			if(current==0)
			{
				nSteps[n] <- 'n'	
			}
			else
			{
				if(current>0)
				{
					nSteps[n] <- paste('n-',affix,sep='')
				}
				else
				{
					if(current<0)
					{
						nSteps[n] <- paste('n+',affix,sep='')
					}
				}
			
			}
			n = n + 1
		}
	}
}

dt$steps <- nSteps

## do markov chain sequences
n = 1
markovContrast <- numeric(dim(dt)[1])
for(i in c('self','other'))
{
	for(j in c(1,0))
	{
		foo <- subset(dt,dt$cond==i & dt$givenDopa == j)
		sequence <- foo$steps
		sequencMatr <- createSequenceMatrix(sequence,sanitize=FALSE)
		mcFitMLE <- markovchainFit(data=sequence)
		ss <- steadyStates(mcFitMLE[[1]])
		ps <- as.vector(ss)
		names <- colnames(ss)
		for(k in 1:dim(foo)[1])
		{
			whichN <- which(names==foo$steps[k])
			markovContrast[n] <- foo$contrast[k]*ps[whichN]
			n=n+1
		}
	}
}

dt$mContrast <- markovContrast
## create mixed models
library(lme4)
mm.markov.full <- lmer(mContrast ~ givenDopa*cond + likeability + (1|id),data=dt,REML=FALSE)
mm.markov.add <- lmer(mContrast ~ givenDopa+cond+likeability + (1|id),data=dt,REML=FALSE)