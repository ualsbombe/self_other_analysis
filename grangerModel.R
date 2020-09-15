rm(list=ls())
graphics.off()

library(lme4)
setwd('/Users/Lau/Desktop/')
dt <- read.csv(file='Granger_bpfilter_order5.csv',sep=';',dec=',')

mm <- lmer(Granger ~ dir*freq + (1|group/id/trial),data=dt,verbose=TRUE)
