## ordinal logistic regression example
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)


dat <-  read.dta('http://www.ats.ucla.edu/stat/data/ologit.dta')
head(dat)

## one at a time, table apply, pared, and public
lapply(dat[,c('apply','pared','public')],table)
## three way cross tabs (xtabs) and flatten the table

ftable(xtabs(~public + apply + pared, data=dat))

summary(dat$gpa)
sd(dat$gpa)

ggplot(dat,aes(x=apply,y=gpa)) +
	geom_boxplot(size=.75) +
	geom_jitter(alpha = .5) +
	facet_grid(pared~public,margins=TRUE) +
	theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust=1))
	
m <- polr(apply ~ pared + public + gpa,data =dat,Hess=TRUE)