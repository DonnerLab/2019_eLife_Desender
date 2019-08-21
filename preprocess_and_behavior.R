#This file uses R to load the data of Experiment 1, preprocess it, perform the behavioral analysis, and then save the data to make it hddm-ready.

setwd("...")
rm(list=ls())

#Some stuff for visuals and packages etc.
mycol3 <- c('black','grey','yellow');mycol2 <- c('black','grey');par(family="sans",pch=16);library(lattice); library(reshape); library(matrixStats);library(effects);library(lme4);library(MALDIquant);library(MASS);library(scales)
cexkl <- 1.5;cexgr <- 2;lwdgr <- 3;

#code to plot error bars #
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
#Same but for horizontal error bars
error.bar.horiz <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x+upper,y, x-lower, y, angle=90, code=3, length=length, ...)
}
##Code for VIFs with mixed models##
vif.mer <- function (fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

#Load the data, these are stored here https://osf.io/83x7c/
Data = read.table('Experiment1.txt',header=T)

#coding error for pp 21
Data$sub[Data$sub==21&Data$age==20] <- 22

#offsetdot codes WHEN in time a rt was being made
Data$rt_in_time <- Data$offsetdot

#Code some previous trial variables
Data$prevrt[2:length(Data$rt)] <- Data$rt[1:length(Data$rt)-1]
Data$prevcor[2:length(Data$cor)] <- Data$cor[1:length(Data$cor)-1]
Data$prevprevcor[3:length(Data$cor)] <- Data$cor[1:(length(Data$cor)-2)]
Data$prevconf[2:length(Data$conf)] <- Data$conf[1:length(Data$conf)-1]
Data$prevcoh[2:length(Data$coh_single)] <- Data$coh_single[1:length(Data$coh_single)-1]
Data$prevvc[2:length(Data$vc)] <- Data$vc[1:length(Data$vc)-1]
Data$follconf <- NA
Data$follconf[1:length(Data$conf)-1] <- Data$conf[2:length(Data$conf)]
Data$postcor <- NA
Data$postcor[1:length(Data$cor)-1] <- Data$cor[2:length(Data$cor)]

#create a variable coding stimulus
Data$stimulus <- 1
Data$stimulus[Data$resp==2&Data$cor==1] <- 2
Data$stimulus[Data$resp==1&Data$cor==0] <- 2

#exclude training etc.
Data <- subset(Data,include==1)
Data <- subset(Data, block>3)

#recode as numeric
Data$coh <- as.numeric(Data$coh_single)
Data$vc <- as.factor(Data$vc)
Data$t2 <- as.factor(Data$t2duration)

#Check whether all makes sense
table(Data$sub)
table(Data$block,Data$sub)
table(Data$include,Data$block)
table(Data$sub,Data$cor,Data$t2duration)
table(Data$cor)
table(Data$conf)

#Code conf so that t1 and t2 have the same interpretation
Data$conf[Data$t2==0] <- Data$conf[Data$t2==0]+3

table(Data$block,Data$include,Data$cor)

#exclude nonsense RTs < .2
Data <- subset(Data, rt > .2) 

#Check meanerr
meanerr <- tapply(Data$cor,list(Data$sub),mean)
densityplot(~meanerr) #M = 75%
meancj <- tapply(Data$conf,list(Data$sub),mean)
densityplot(~meancj)  
#overall RT level
meanrt1 <- with(Data,aggregate(rt, by=list(sub),mean,na.rm=T))
densityplot(~meanrt1[2]) # mean=1.19s sd=0.23s

#Plot raw stuff
par(mfrow=c(2,4))
for(i in 1:N){
  temp <- subset(Data,sub==i)
  temp <- subset(temp, t2 != 0)
  plot(temp$rt,frame=F,main=paste("Subject",i),ylim=c(0,3))  
  tempACC <- with(temp,aggregate(cor,by=list(block,t2duration),mean))
  plot(tempACC$x[tempACC$Group.2==1],frame=F,ylab='accuracy',ylim=c(0,1),type='b',main=paste("Subject",i));abline(h=.5,lty=2)
  lines(tempACC$x[tempACC$Group.2==2],col='blue',type='b')
  tempBIAS <- with(temp,aggregate(resp,by=list(block,t2duration),mean))
  plot(tempBIAS$x[tempBIAS$Group.2==1]-1,frame=F,ylab='response bias',ylim=c(0,1),type='b',main=paste("Subject",i));abline(h=.5,lty=2)
  lines(tempBIAS$x[tempBIAS$Group.2==2]-1,col='blue',type='b')
  cjAcc <- with(temp,aggregate(cor,by=list(conf),mean))  
  plot(cjAcc$x~cjAcc$Group.1,ylim=c(0,1),xlim=c(1,6),ylab='accuracy',xlab='confidence',frame=F,type='b')
}

#6 and 30 are out because they lack data in one the confidence cells!
Data <- subset(Data, sub!=6)
Data <- subset(Data, sub!=30)

#12, 15, 21, 29 are off in t0 (immediate confidence), but t0 data is not used here so we'll leave them in

#Some useful variables
N <- length(table(Data$sub))
subs <- with(Data,aggregate(sub,by=list(sub),mean))$x


#Check statistical confidence signatures on these data
DataCor <- subset(Data,cor==1)
DataErr <- subset(Data,cor==0)
cohrt <- with(DataCor,aggregate(rt,by=list(sub,coh),median))
names(cohrt) <- c("Subject","coh","RT")
cohrt <- cast(cohrt,Subject~coh)

cohrtERR <- with(DataErr,aggregate(rt,by=list(sub,coh),median))
names(cohrtERR) <- c("Subject","coh","RT")
cohrtERR <- cast(cohrtERR,Subject~coh)

cohcor <- with(Data,aggregate(cor,by=list(sub,coh),mean))
names(cohcor) <- c("Subject","coh","cor")
cohcor <- cast(cohcor,Subject+experiment~coh)

cohcjCor <- with(DataCor,aggregate(conf,by=list(sub,coh),mean))
names(cohcjCor) <- c("Subject","coh","cj")
cohcjCor <- cast(cohcjCor,Subject~coh)

cohcjErr <- with(DataErr,aggregate(conf,by=list(sub,coh),mean))
names(cohcjErr) <- c("Subject","coh","cj")
cohcjErr <- cast(cohcjErr,Subject~coh)

#RTs on correct vs error trials
tiff(file="RTsOnCorVsErr.tiff", res=350, width=300*(350/72), height=400*(350/72))
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- snrrt[,3:7]*1000
y <- snrrtERR[,3:7]*1000

means <- sapply(x, mean);meansy <- sapply(y, mean,na.rm=T);n<- length(x)
stripchart(x, ylim=c(200,2500), xlim=c(0,as.numeric(max(names(x)))+.1), vertical = TRUE, col="white",frame=F, xaxt='n')
mtext("Reaction times (ms)",2,at=600,line=2.5)
axis(1,at=as.numeric(names(x))+.1,labels=as.numeric(names(x)))
for(i in 1:20) abline(h=(000+(i*250)),col='lightgrey',lwd=0.8)  
mtext("Coherence",1,2.5)
stripchart(x,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)
stripchart(y,at=as.numeric(names(y))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg="grey",col='white',cex=1.5)
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=21,bg="black",col="white",cex=1.5)
lines(as.numeric(names(x))+.05,means,type='b')
for(i in 1:n) points(as.numeric(names(y))[i]+.05,meansy[i],pch=24,bg="black",col="white",cex=1.5)
lines(as.numeric(names(y))+.05,meansy,type='b')
legend("bottomleft",legend=c("Correct","Error"),pch=c(16,17),lty=1,bty = "n")
dev.off()



#Confidence on cor vs errs
tiff(file="ConfOnCorVsErr.tiff", res=350, width=300*(350/72), height=390*(350/72))
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- cohcjCor[,3:7]
xErr <- cohcjErr[,3:7]

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(1,6), xlim=c(0,max(as.numeric(names(x)))+.1), vertical = TRUE, col="white",frame=F,axes=F)
axis(2,at=1:6,labels=c('error','error','error','correct','correct','correct'))
mtext(c('sure','probably','guess','guess','probably','sure'),2,at=1:6,line=1.8)
mtext("Confidence",2,at=3.5,line=3)
axis(1,at=as.numeric(names(x))+.05,labels=names(x),cex.axis=.9)
for(i in 1:20) abline(h=(0+(i*.5)),col='lightgrey')  
mtext("Coherence",1,2.5)
stripchart(x,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg='grey',col='white',cex=1.5)
stripchart(xErr,at=as.numeric(names(xErr))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg='grey',col='white')
lines(as.numeric(names(x))+.05,means,type='b')
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=21,bg="black",col="white",cex=1.5)
means <- sapply(xErr, mean,na.rm=T);
lines(as.numeric(names(x))+.05,means,type='b',pch=2)
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=24,bg="black",col="white",cex=1.5)
legend("bottomleft",legend=c("Correct","Error"),pch=c(16,17),lty=1,bty = "n")
dev.off()


}
setwd("C:/Users/kobe/Documents/Projecten/Slowing after confidence/OldData/hDDM/varConfPostSlowing")

#SAVE DATA TO TEST POST-CONF SLOWING HYPOTHESIS (no iri=0)
hddmData <- Data
DataTzero <- subset(hddmData,t2==0)
#12, 15, 21, 29 are off in t0 <- non differs from chance
DataTzero <- subset(DataTzero,sub!=12)
DataTzero <- subset(DataTzero,sub!=15)
DataTzero <- subset(DataTzero,sub!=21)
DataTzero <- subset(DataTzero,sub!=29)
hddmData <- subset(hddmData,t2 != 0)
hddmData <- subset(hddmData,prevconf > 0)
hddmData <- subset(hddmData,follconf > 0)
hddmData <- subset(hddmData,withinblocktrial > 1) #delete first trials (no prevconf)
hddmData <- subset(hddmData,withinblocktrial < 60) #delete last trials (no postconf)
hddmData$subj_idx <- hddmData$sub
hddmData$resp <- hddmData$resp-1
hddmData$stimulus <- -1
hddmData$stimulus[hddmData$resp==1 & hddmData$cor == 1] <- 1
hddmData$stimulus[hddmData$resp==0 & hddmData$cor == 0] <- 1
hddmData$response <- hddmData$resp
Data <- hddmData
# write.table(hddmData,"varConfcombined3_hddmSLOWING.csv",quote=F,row.names=F,sep=',')

#0. Inspect Data per participants
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(3,2))
for(s in 1:N){   
i = subs[s]
temp <- subset(Data,sub==i)
plot(temp$rt,frame=F,main=paste('Subject',s),ylab="RTs",xlab='trial',ylim=c(0,3),cex.lab=2.5,cex.axis=2.5,cex.main=2.5)
tempBlock <- with(temp,aggregate(cor,by=list(block),mean,na.rm=TRUE))
plot(tempBlock$x,frame=F,main="",ylab='accuracy',xlab='block',ylim=c(.4,1),cex.lab=2.5,cex.axis=2.5,cex.main=2.5);abline(h=.5,lty=2)
#  plot(temp$conf,frame=F,main="",ylab="ratings")
}


##########
#Analysis#
library(lmerTest)
Data$coh <- as.factor(Data$coh)
Data$cor <- as.factor(Data$cor)
#RT
DataCor <- subset(Data,cor==1)
DataErr <- subset(Data,cor==0)
fit <- lmer(rt~coh + (1|sub),data=DataCor)
fit1 <- lmer(rt~coh + (coh|sub),data=DataCor)
anova(fit,fit1)
anova(fit1)

#cor
fit <- glmer(cor~coh + (1|sub),data=Data,family='binomial')
fit1 <- glmer(cor~coh + (coh|sub),data=Data,family='binomial')
anova(fit,fit1)
Anova(fit)

#conf
fit <- lmer(conf~coh*cor + (1|sub),data=Data)
fit1a <- lmer(conf~coh*cor + (coh|sub),data=Data)
fit1b <- lmer(conf~coh*cor + (cor|sub),data=Data)
anova(fit,fit1a)
anova(fit,fit1b)
fit2 <- lmer(conf~coh*cor + (coh+cor|sub),data=Data)
#fit2a <- lmer(conf~coh*cor + (coh*cor|sub),data=Data)
anova(fit2)

Data$t2 <- as.factor(Data$t2);Data$coh <- as.factor(Data$coh)
fit3 <- lmer(conf~coh*cor*t2 + (coh+cor|sub),data=Data)
fit3_b <- lmer(conf~coh*cor*t2 + (coh+cor+t2|sub),data=Data)
anova(fit3,fit3_b)
anova(fit3_b)
data.frame(effect('t2',fit3_b))

#Post-hoc linear contrast, default=0%,error trials,t1
E <- matrix(c(0,1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),1) #error trials 
C <- matrix(c(0,1,2,3,4,1,0,1,2,3,4,0,0,0,0,0,0,0,0,0),1) #correts trials 
t <- glht(fit3_b, linfct = E)
t <- glht(fit3_b, linfct = C)
summary(t)

#Revision 1: show the BF for no effect of post-decisional evidence
library(BayesFactor)
confi <- with(Data,aggregate(conf,by=list(sub,coh,cor,t2),mean))
names(confi) <- c('sub','coh','cor','t2','conf')
confi$coh <- as.factor(confi$coh);confi$t2 <- as.factor(confi$t2);confi$sub <- as.factor(confi$sub);confi$cor <- as.factor(confi$cor)
fit <- anovaBF(conf~coh*t2*cor+sub,data=confi,whichRandom="sub")
fit[18]/fit[17]


#Revision 1: 4A. HOW GOOD ARE THE DATA: ANALYSIS ON SINGLE SUBJECT LEVEL
library(car)
rt_p <- rep(NA,N);acc_p <- rep(NA,N);conf_p <- rep(NA,N);
rt_b <- rep(NA,N);acc_b <- rep(NA,N);conf_b <- rep(NA,N);
rt_halve <- rep(NA,N);acc_halve <- rep(NA,N);conf_halve <- rep(NA,N);
for(i in 1:N){
  temp <- subset(Data,sub==subs[i])
  temp$cor <- as.numeric(temp$cor)
  temp$halve <- "first";
  temp$halve[round(length(temp$rt)/2):length(temp$rt)] <- "second"
  temp$halve <- as.factor(temp$halve)
  tempCor <- subset(temp,cor==1)
  
  fit <- lm(rt~coh*halve,data=tempCor)
  rt_p[i] <- anova(fit)[1,5]
  rt_halve[i] <- anova(fit)[3,5]
  rt_b[i] <- summary(fit)$coefficients[2,1]
  fit <- glm(cor~coh*halve,data=temp)
  acc_p[i] <- Anova(fit)[1,3]
  acc_halve[i] <- Anova(fit)[3,3]
  acc_b[i] <- summary(fit)$coefficients[2,1]
  fit <- glm(conf~coh*halve,data=temp)
  conf_p[i] <- Anova(fit)[1,3]
  conf_halve[i] <- Anova(fit)[3,3]
  conf_b[i] <- summary(fit)$coefficients[2,1]
}
paste("Significant effect of coherence on RT for ", length(which(rt_p<.05))," out of ",N," subs; correct sign for ", length(which(rt_b<0)))
paste("Significant effect of coherence on accuracy for ", length(which(acc_p<.05))," out of ",N," subs; correct sign for ", length(which(acc_b>0)))
paste("Significant effect of coherence on confidence for ", length(which(conf_p<.05))," out of ",N," subs; correct sign for ", length(which(conf_b>0)))

paste("Coherence*halve for rt was p<05 for ", length(which(rt_halve<.05))," out of ",N," subs")
paste("Coherence*halve for accuracy was p<05 for ", length(which(acc_halve<.05))," out of ",N," subs")
paste("Coherence*halve for confidence was p<05 for ", length(which(conf_halve<.05))," out of ",N," subs")
#who are these people
subs[which(rt_p>.05|acc_p>.05|conf_p>.05)] #no scaling with coherence
subs[which(rt_halve<.05|acc_halve<.05|conf_halve<.05)] #interaction with halve


#Revision 1: 4A. HOW STABLE ARE THE DATA: INCLUDE BLOCK
DataCor$block <- as.factor(DataCor$block)
Data$block <- as.factor(Data$block)
fit1 <- lmer(rt~coh*block + (1|sub),data=DataCor)
fit1 <- glmer(cor~coh*block + (1|sub),data=Data,family='binomial')
Anova(fit1)
plot(effect('coh:block',fit1))
#Revision 1: 4A. HOW STABLE ARE THE DATA: check at the subject level
rt_p <- rep(NA,N);acc_p <- rep(NA,N);conf_p <- rep(NA,N);meta_p <- rep(NA,N);
for(i in 1:N){
  temp <- subset(Data,sub==subs[i])
  temp$cor <- as.numeric(temp$cor)
  temp$trial <- 1:length(temp$trial)
  
  fit <- lm(rt~trial,data=temp)
  rt_p[i] <- anova(fit)[1,5]
  fit <- glm(cor~trial,data=temp)
  acc_p[i] <- Anova(fit)[1,3]
  fit <- glm(conf~trial,data=temp)
  conf_p[i] <- Anova(fit)[1,3]
  temp$conf <- as.factor(temp$conf)
  fit <- glm(cor~conf*trial,data=temp)
  meta_p[i] <- Anova(fit)[3,3]
}
paste("Significant change in RT for ", length(which(rt_p<.05))," out of ",N," subs")
paste("Significant change in accuracy for ", length(which(acc_p<.05))," out of ",N," subs")
paste("Significant change in confidence for ", length(which(conf_p<.05))," out of ",N," subs")
paste("Significant change in cor~conf for ", length(which(meta_p<.05))," out of ",N," subs")


#post-hoc seperations
t2 <- subset(Data,t2==1)
fitFU <- lmer(conf~coh*cor + (coh+cor|sub),data=t2)
anova(fitFU)


#Autocorrelation
fit <- lmer(rt~prevrt + (1|sub),data=Data)
fita <- lmer(rt~prevrt + (prevrt|sub),data=Data)
anova(fit,fita)
anova(fita)

fit <- glmer(cor~prevcor + (1|sub),data=Data,family='binomial')
fita <- glmer(cor~prevcor + (prevcor|sub),data=Data,family='binomial')
anova(fit,fita)
Anova(fit)


fit <- lmer(conf~prevconf + (1|sub),data=Data)
fita <- lmer(conf~prevconf + (prevconf|sub),data=Data)
anova(fit,fita)
anova(fita)


#### Recode previous and following cj in three bins, guess is the reference
Data$prevcj[Data$prevconf==1] = 'C_error' #sure error
Data$prevcj[Data$prevconf==2] = 'C_error' #prob error
Data$prevcj[Data$prevconf==3] = 'B_guess' #guess error
Data$prevcj[Data$prevconf==4] = 'B_guess' #guess cor => this becomes the reference, because lowest value!!! (change it?)
Data$prevcj[Data$prevconf==5] = 'A_correct' #prob cor
Data$prevcj[Data$prevconf==6] = 'A_correct' #sure cor

#### Recode previous and following cj in three bins, guess is the reference
Data$follcj[Data$follconf==1] = 'C_error' #sure error
Data$follcj[Data$follconf==2] = 'C_error' #prob error
Data$follcj[Data$follconf==3] = 'B_guess' #guess error
Data$follcj[Data$follconf==4] = 'B_guess' #guess cor => this becomes the reference, because lowest value!!! (change it?)
Data$follcj[Data$follconf==5] = 'A_correct' #prob cor
Data$follcj[Data$follconf==6] = 'A_correct' #sure cor

#########################################
# "Simulate" bound based on actual data #
# Data$bound <- NA
# Data$bound[1] <- 1
# for(i in 2:length(Data$rt)){
#   if(Data$prevcj[i]=='A_correct'){ Data$bound[i] <- Data$bound[i-1] -.045} 
#   if(Data$prevcj[i]=='B_guess'){ Data$bound[i] <- Data$bound[i-1] +.08}
#   if(Data$prevcj[i]=='C_error'){ Data$bound[i] <- Data$bound[i-1] +.2}
# }
# plot(Data$bound)
# with(Data,aggregate(bound,by=list(prevcj),mean)) #only condition on previous confidence gives flawed result
# with(Data,aggregate(bound,by=list(prevcj),mean))$x - with(Data,aggregate(bound,by=list(follcj),mean))$x #statistical control
# with(Data,aggregate(bound,by=list(prevcj,follcj),mean)) #only selecting post-same-conf trials works
# with(Data,aggregate(bound,by=list(prevcj,postcor),mean)) #only conditio on postcor gives flawed result

#crosstabs
table(Data$sub,Data$prevcj)
table(Data$sub,Data$follcj)

table(Data$sub,Data$prevcj,Data$t2) #10 is out for t1, 3 and 15 for t2
table(Data$sub,Data$follcj,Data$t2) #10 is out for t1, 3, 13 and 15, 28 for t2

table(Data$sub,Data$prevcj,Data$prevcor) #2,15,28 are out for postcor!!
table(Data$sub,Data$follcj,Data$prevcor) #and 26?

#show trial numbers
x <- table(Data$prevcj,Data$sub)[c(3,2,1),]
ebar <- barplot(x,ylab='# trials',xlab='participants',main="Experiment 1",col=c(rgb(0,1,0,.5),rgb(1,0,0,.5),rgb(0,0,1,.5)),border='white')
for(i in 1:N) text(x=ebar[i],y=sum(x[1:2,i])+x[3,i]/2,x[3,i], srt = 90)
for(i in 1:N) text(x=ebar[i],y=x[1,i]+x[2,i]/2,x[2,i], srt = 90)  
for(i in 1:N) text(x=ebar[i],y=x[1,i]/2,x[1,i], srt = 90) 



################################################################################
###  REPLICATE PREVIOUS BEHAVIORAL ANALYSIS
################################################################################
#FIRST show that behavior scales with SNR THEN show that previous confidence adds somethin
#(seperately for exp1 and exp2)
Data$snr <- as.factor(Data$coh)
Data$cj <- as.numeric(Data$conf)
DataCor <- subset(Data,cor==1)
DataErr <- subset(Data,cor==0)

#SNR - individual differences plot
DataCor$snr <- as.factor(DataCor$coh_single)
snrrt <- with(DataCor,aggregate(rt,by=list(sub,snr),median))
names(snrrt) <- c("Subject","snr","RT")
snrrt <- cast(snrrt,Subject~snr)

snrrt_sd <- with(DataCor,aggregate(rt,by=list(sub,snr),sd))
names(snrrt_sd) <- c("Subject","snr","RT")
snrrt_sd <- cast(snrrt_sd,Subject~snr)
ntrials <- table(DataCor$sub,DataCor$snr)
snrrt_sd[,2:6] <- snrrt_sd[,2:6]/sqrt(ntrials) 

snrrtERR <- with(DataErr,aggregate(rt,by=list(sub,snr),median))
names(snrrtERR) <- c("Subject","snr","RT")
snrrtERR <- cast(snrrtERR,Subject~snr)

snrcor <- with(Data,aggregate(cor,by=list(sub,snr),mean))
names(snrcor) <- c("Subject","snr","cor")
snrcor <- cast(snrcor,Subject~snr)

DataCorTemp <- subset(DataCor, t2==2) #1=postevidence, 2=blank
snrcjCor <- with(DataCorTemp,aggregate(cj,by=list(sub,snr),mean))
names(snrcjCor) <- c("Subject","snr","cj")
snrcjCor <- cast(snrcjCor,Subject~snr)

DataErrTemp <- subset(DataErr, t2==2) #1=postevidence, 2=blank
snrcjErr <- with(DataErrTemp,aggregate(cj,by=list(sub,snr),mean))
names(snrcjErr) <- c("Subject","snr","cj")
snrcjErr <- cast(snrcjErr,Subject~snr)


#RTs and Accuracy
setwd("C:/Users/kobe/Documents/Projecten/Slowing after confidence/OldData/hDDM/varConfPostSlowing")
tiff(file="figureswcond/snr1_WITH_EBAR.tiff", res=350, width=300*(350/72), height=400*(350/72))
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- snrrt[,2:6]
y <- snrcor[,2:6]

#x_pred <- colMeans(RTpreds)*1000
#y_pred <- colMeans(ACCpreds)

library(scales);y_rescale <- data.frame(matrix(NaN,dim(y)[1],dim(y)[2]));
for(i in 1:dim(y)[2]) y_rescale[,i] <- scales::rescale(y[,i],to=c(-325,0),from=c(.5,1))
means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(-400,2200), xlim=c(0,as.numeric(max(names(x)))+.1), vertical = TRUE, col="white",frame=F, axes=F)
axis(2,at=seq(400,2200,200),labels=seq(400,2200,200))
mtext("Reaction times (ms)",2,at=1200,line=2.5)
axis(4,at=seq(-325,0,325/2),labels=seq(50,100,25))
mtext("Accuracy (%)",4,at=-162.5,line=2.5)
for(i in 1:20) abline(h=(200+(i*200)),col='lightgrey',lwd=0.8)  
for(i in 1:3) abline(h=(-325-(325/2)+(i*325/2)),col='lightgrey',lwd=0.8)  
axis(1,at=as.numeric(names(x))+.05,labels=as.numeric(names(x))*100,cex.axis=.9)
mtext("Coherence (%)",1,2.5)
#stripchart(x*1000,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)
dotPos <- x;
for(i in 1:n) dotPos[,i] <-  stripchart_kobe(x=x[,i]*1000,at=as.numeric(names(x))[i]+.05)
for(i in 1:n) x[,i] <- sort(x[,i])
for(i in 1:N) for(j in 1:n) error.bar(dotPos[i,j],x[i,j]*1000,snrrt_sd[i,j+1]*1000,col='grey')

for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i]*1000,pch=21,bg="black",col="white",cex=1.5)
stripchart(y_rescale,at=as.numeric(names(y))+.05,method = 'jitter',jitter=.01,add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)
meansy <- sapply(y_rescale, mean,na.rm=T);
for(i in 1:n) points(as.numeric(names(y))[i]+.05,meansy[i],pch=21,bg="black",col="white",cex=1.5)
lines(as.numeric(names(x))+.05,means*1000,type='b')
lines(as.numeric(names(y))+.05,meansy,type='b')

# lines(as.numeric(names(x))+.05,x_pred,type='p',pch=4,col='green',lwd=2)
# lines(as.numeric(names(y))+.05,rescale(y_pred,to=c(-325,0),from=c(.5,1)),type='p',pch=4,col='green',lwd=2)
dev.off()

#RTs on correct vs error trials
tiff(file="figureswcond/RTcorerr.tiff", res=350, width=300*(350/72), height=400*(350/72))
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- snrrt[,2:6]
y <- snrrtERR[,2:6]

means <- sapply(x, mean);meansy <- sapply(y, mean,na.rm=T);n<- length(x)
stripchart(x, ylim=c(200,2200), xlim=c(0,as.numeric(max(names(x)))+.05), vertical = TRUE, col="white",frame=F, xaxt='n')
mtext("Reaction times (ms)",2,at=600,line=2.5)
axis(1,at=as.numeric(names(x))+.05,labels=as.numeric(names(x))*100)
for(i in 1:20) abline(h=(0+(i*250)),col='lightgrey',lwd=0.8)  
mtext("Coherence",1,2.5)
stripchart(x*1000,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)
stripchart(y*1000,at=as.numeric(names(y))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg="grey",col='white',cex=1.5)
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i]*1000,pch=21,bg="black",col="white",cex=1.5)
lines(as.numeric(names(x))+.05,means*1000,type='b')
for(i in 1:n) points(as.numeric(names(y))[i]+.05,meansy[i]*1000,pch=24,bg="black",col="white",cex=1.5)
lines(as.numeric(names(y))+.05,meansy*1000,type='b')
dev.off()


#Confidence
tiff(file="figureswcond/coh_cj.tiff", res=350, width=300*(350/72), height=390*(350/72))
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- snrcjCor[,2:6]
xErr <- snrcjErr[,2:6]

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(1,6), xlim=c(0,max(as.numeric(names(x)))+.1), vertical = TRUE, col="white",frame=F,axes=F,
           main="all")
axis(2,at=1:6,labels=rep('',6))
mtext(c('error','error','error','correct','correct','correct'),2,at=1:6,line=.8)
mtext(c('sure','probably','guess','guess','probably','sure'),2,at=1:6,line=1.8)
mtext("Confidence",2,at=3.5,line=3)
axis(1,at=as.numeric(names(x))+.05,labels=as.numeric(names(x))*100, cex.axis=.9)
for(i in 1:20) abline(h=(0+(i*.5)),col='lightgrey')  
mtext("Coherence (%)",1,2.5)
stripchart(x,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg='grey',col='white',cex=1.5)
stripchart(xErr,at=as.numeric(names(xErr))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg='grey',col='white')
lines(as.numeric(names(x))+.05,means,type='b')
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=21,bg="black",col="white",cex=1.5)
means <- sapply(xErr, mean,na.rm=T);
lines(as.numeric(names(x))+.05,means,type='b',pch=2)
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=24,bg="black",col="white",cex=1.5)
legend("bottomleft",legend=c("Correct","Error"),pch=c(16,17),lty=1,bty = "n")
dev.off()

#X. RT predicts accuracy (tobi) 
{
  Data$RTquint <- NA;SubsID <- names(table(Data$sub))
  for(j in 1:length(table(Data$sub))){
    tempSub <- subset(Data, sub==SubsID[j])
    quants <- quantile(tempSub$rt,seq(0,1,.1))
    
    for(i in 1:length(quants)-1){
      Data$RTquint[Data$sub==SubsID[j] & Data$rt>=quants[i] & Data$rt <= quants[i+1]] <- i
    }
  }
  table(is.na(Data$RTquint))
  boxplot(Data$rt~Data$RTquint)
  with(Data,aggregate(rt,by=list(RTquint,sub),median))
  
  #Show accuracy per RT level
  accsnrrt1 <- with(Data,aggregate(cor,by=list(RTquint,sub),mean))
  names(accsnrrt1) <- c('rtquint','sub','acc')
  accsnrrt <- cast(accsnrrt1,sub~rtquint)
  
  rtsnrrt <- with(Data,aggregate(rt,by=list(RTquint,sub),median))
  names(rtsnrrt) <- c('rtquint','sub','mrt')
  rtsnrrt <- cast(rtsnrrt,sub~rtquint)
  
  tiff(file=paste('figureswcond/rtacc.tiff',sep=""), res=350, width=350*(350/72), height=350*(350/72))
  plot(colMeans(rtsnrrt)*1000,colMeans(accsnrrt)*100,xlim=c(500,2200),ylim=c(0,100),type='b',pch=20,frame=F,ylab='Accuracy (%)',xlab='Reaction time (ms)')
  for(i in 1:20) abline(h=(40+(i*5)),col='lightgrey')  
  error.bar(colMeans(rtsnrrt)*1000,colMeans(accsnrrt)*100,(colSds(as.matrix(accsnrrt))*100)/sqrt(dim(accsnrrt)[2]))
  dev.off()
  
  tiff(file=paste('figureswcond/rtacc_raw.tiff',sep=""), res=350, width=350*(350/72), height=350*(350/72))
  par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
  x <- accsnrrt[,11:-1:2]*100
  y <- rtsnrrt[,11:-1:2]*1000
  
  means <- sapply(x, mean);n<- length(x)
  stripchart(x, ylim=c(0,100), xlim=c(2300,500), vertical = TRUE, col="white",frame=F,xaxt='n',
             main="all")
  mtext("Accuracy (%)",2,at=50,line=3)
  axis(1)
  for(i in 1:20) abline(h=(-20+(i*20)),col='lightgrey')  
  mtext("Reaction time (ms)",1,2.5)
  #for(i in 1:n) points(y=x[,i],x=y[,i],pch=21,bg='grey',col='white',cex=1.5)
  stripchart(x,at=colMeans(y),jitter=30, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg='grey',col='white',cex=1.5)
  lines(colMeans(y),means,type='b')
  for(i in 1:n) points(colMeans(y)[i],means[i],pch=21,bg="black",col="white",cex=1.5)
  dev.off()
  
  
  fit <- lmer(acc~rtquint + (1|sub),data=accsnrrt1)
  fit1 <- lmer(acc~rtquint + (rtquint|sub),data=accsnrrt1)
  anova(fit,fit1)
  summary(fit1)
}


#X. Seperate axis for accuracy (per label) and frequency
Data$ones <- rep(1,length(Data$sub))
tempDat <- subset(Data, t2==1) #1=postevidence, 2=blank
#tempDat <- DataTzero

AccPerLab <- with(tempDat, aggregate(cor,by=list(sub,conf),mean))
names(AccPerLab) <- c("Subject","cj","meanACC")
AccPerLab$meanACC <- AccPerLab$meanACC*100
AccPerLab <- cast(AccPerLab,Subject~cj)
t.test(AccPerLab$`6`,mu=.5)
t.test(AccPerLab$`1`,mu=.5)

for(i in 1:6) print(t.test(AccPerLab1[,i+1],AccPerLab2[,i+1],paired=T))

FreqPerLab <- with(tempDat, aggregate(ones,by=list(sub,conf),sum))
names(FreqPerLab) <- c("Subject","cj","meanFreq")
FreqPerLab <- cast(FreqPerLab,Subject~cj)
TotalFreq <- with(tempDat, aggregate(ones,by=list(sub),sum))
FreqPerLab[,c(2:7)] <- (FreqPerLab[,c(2:7)]/TotalFreq$x)*100
FreqPerLab <- as.matrix(FreqPerLab[,c(2:7)],23,6)
dotSize <- rescale(FreqPerLab,to=c(.75,3.2))

tiff(file="figureswcond/resolution_t1.tiff", res=350, width=300*(375/72), height=325*(350/72))
x <- AccPerLab[,c(2:7)]
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
means <- sapply(x, mean,na.rm=T);n<- length(x)
stripchart(x, ylim=c(0,100),xlim=c(.5,6.5),vertical = TRUE, col="white",frame=F,xaxt='n', xlab='Confidence',ylab='Accuracy %',
           main=c('t1'))
axis(1,at=1:6,labels=c('sure','probably','guess','guess','probably','sure'),cex.axis=.8)
mtext(c('error','error','error','correct','correct','correct'),1,at=1:6,line=1.8,cex=.8)
for(i in 1:20) abline(h=(-20+(i*20)),col='lightgrey')  
for(i in 1:n) stripchart(x[,i],at=as.numeric(names(x))[i],jitter=0.25, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=dotSize[,i])
for(i in 1:n) points(as.numeric(names(x))[i],means[i],pch=21,bg="black",col="white",cex=colMeans(dotSize,na.rm=T)[i])
lines(as.numeric(names(x)),means,type='l')
legend(4.3,35,legend=c("0.5%","40%","80%"),pch=19,col='grey',box.lty=0,pt.cex=seq(.75,3.2,length.out=3))
dev.off()

#################################
## T ZERO STUFF #################
tiff(file="figureswcond/resolution_t0.tiff", res=350, width=300*(375/72), height=325*(350/72))
x <- AccPerLab[,c(2:4)]
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
means <- sapply(x, mean,na.rm=T);n<- length(x)
stripchart(x, ylim=c(0,100),xlim=c(.5,3.5),vertical = TRUE, col="white",frame=F,xaxt='n', xlab='Confidence',ylab='Accuracy %',
           main=c('t0'))
axis(1,at=1:3,labels=c('guess','probably','sure'),cex.axis=.8)
mtext(c('correct','correct','correct'),1,at=1:3,line=1.8,cex=.8)
for(i in 1:20) abline(h=(-20+(i*20)),col='lightgrey')  
for(i in 1:n) stripchart(x[,i],at=as.numeric(names(x))[i]-3,jitter=0.25, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=dotSize[,i])
for(i in 1:n) points(as.numeric(names(x))[i]-3,means[i],pch=21,bg="black",col="white",cex=colMeans(dotSize)[i])
lines(as.numeric(names(x))-3,means,type='l')
legend(2,35,legend=c("0.5%","40%","80%"),pch=19,col='grey',box.lty=0,pt.cex=seq(.75,3.2,length.out=3))
dev.off()
t.test(AccPerLab[,2],mu=.5)
t.test(AccPerLab[,3],mu=.5)
t.test(AccPerLab[,4],mu=.5)

#Confidence
DataTzeroCor <- subset(DataTzero,cor==1)
snrcjCor <- with(DataTzeroCor,aggregate(conf,by=list(sub,coh_single),mean))
names(snrcjCor) <- c("Subject","snr","cj")
snrcjCor <- cast(snrcjCor,Subject~snr)

DataTzeroErr <- subset(DataTzero,cor==0)
snrcjErr <- with(DataTzeroErr,aggregate(conf,by=list(sub,coh_single),mean))
names(snrcjErr) <- c("Subject","snr","cj")
snrcjErr <- cast(snrcjErr,Subject~snr)

tiff(file="figureswcond/coh_cjTZERO.tiff", res=350, width=300*(350/72), height=390*(350/72))
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- snrcjCor[,2:6]
xErr <- snrcjErr[,2:6]

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(4,6), xlim=c(0,max(as.numeric(names(x)))+.1), vertical = TRUE, col="white",frame=F,axes=F,
           main="t0")
axis(2,at=4:6,labels=rep('',3))
mtext(c('correct','correct','correct'),2,at=4:6,line=.8)
mtext(c('guess','probably','sure'),2,at=4:6,line=1.8)
mtext("Confidence",2,at=5,line=3)
axis(1,at=as.numeric(names(x))+.05,labels=as.numeric(names(x))*100, cex.axis=.9)
for(i in 1:20) abline(h=(0+(i*.5)),col='lightgrey')  
mtext("Coherence (%)",1,2.5)
stripchart(x,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg='grey',col='white',cex=1.5)
stripchart(xErr,at=as.numeric(names(xErr))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg='grey',col='white')
lines(as.numeric(names(x))+.05,means,type='b')
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=21,bg="black",col="white",cex=1.5)
means <- sapply(xErr, mean,na.rm=T);
lines(as.numeric(names(x))+.05,means,type='b',pch=2)
for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i],pch=24,bg="black",col="white",cex=1.5)
legend(x=.25,y=4.25,legend=c("Correct","Error"),pch=c(16,17),lty=1,bty = "n")
dev.off()







#X. RTs as a function of prevconf
Data$include <- 0
Data$include[Data$prevcj=='C_error'] <- 1
Data$include[Data$prevcj=='A_correct'&Data$follcj=='C_error'] <- 2
temp <- subset(Data,include>0)
ntrials <- table(temp$sub,temp$prevcj)

plot(density(Data$rt[Data$prevcj=='A_correct'&Data$cor==1]),xlab='RTs',ylab='density',main="",xlim=c(-3,3),frame=F,type='n')
polygon(density(Data$rt[Data$prevcj=='A_correct'&Data$cor==1]),col=rgb(0,0,1,.5))
polygon(density(Data$rt[Data$prevcj=='B_guess'&Data$cor==1]),col=rgb(1,0,0,.5))
polygon(density(Data$rt[Data$prevcj=='C_error'&Data$cor==1]),col=rgb(0,1,0,.5))
polygon(density(-Data$rt[Data$prevcj=='A_correct'&Data$cor==0]),col=rgb(0,0,1,.5))
polygon(density(-Data$rt[Data$prevcj=='B_guess'&Data$cor==0]),col=rgb(1,0,0,.5))
polygon(density(-Data$rt[Data$prevcj=='C_error'&Data$cor==0]),col=rgb(0,1,0,.5))

DataCor <- subset(Data,cor==1)
prevcj <- with(DataCor,aggregate(rt,by=list(sub,prevcj),median))
names(prevcj) <- c("Subject","prevcj","RT")
prevcjrt <- cast(prevcj,Subject~prevcj)[,2:4]

follcj <- with(DataCor,aggregate(rt,by=list(sub,follcj),median))
names(follcj) <- c("Subject","follcj","RT")
postcjrt <- cast(follcj,Subject~follcj)[,2:4]

prevcjacc1 <- with(Data,aggregate(cor,by=list(sub,prevcj),mean))
names(prevcjacc1) <- c("Subject","prevcj","acc")
prevcjacc <- cast(prevcjacc1,Subject~prevcj)[,2:4]

follcjacc1 <- with(Data,aggregate(cor,by=list(sub,follcj),mean))
names(follcjacc1) <- c("Subject","follcj","acc")
postcjacc <- cast(follcjacc1,Subject~follcj)[,2:4]

ntrialsCor <- with(DataCor,aggregate(ones,by=list(sub,prevcj),sum))
names(ntrialsCor) <- c("Subject","prevcj","ntrials")
ntrialsCor <- cast(ntrialsCor,Subject~prevcj)[,2:4]
ntrialsAll <- with(Data,aggregate(ones,by=list(sub,prevcj),sum))
names(ntrialsAll) <- c("Subject","prevcj","ntrials")
ntrialsAll <- cast(ntrialsAll,Subject~prevcj)[,2:4]
for(i in 1:3) ntrialsAll[,i] <- rescale(ntrialsAll[,i],to=c(.75,3.5),from=range(ntrialsAll))
for(i in 1:3) ntrialsCor[,i] <- rescale(ntrialsCor[,i],to=c(.75,3.5),from=range(ntrialsCor))

#Query predictions
{
  # setwd("C:/Users/kobe/Documents/Projecten/Variance confidence coupling/model")
  # source('rw_hddmbased.R')
  # setwd("C:/Users/kobe/Documents/Projecten/Slowing after confidence/OldData/hDDM/varConfPostSlowing")
  # 
  # #contruct matrices with v and a for each level of confidence prev and foll
  # a <- c(group$a_intercept,group$a_intercept+group$a_prevguess,group$a_intercept+group$a_preverror,
  #        group$a_intercept,group$a_intercept+group$a_follguess,group$a_intercept+group$a_follerror)
  # v <- c(group$v_intercept,group$v_intercept+group$v_prevguess,group$v_intercept+group$v_preverror,
  #        group$v_intercept,group$v_intercept+group$v_follguess,group$v_intercept+group$v_follerror)
  # t <- group$t
  # labels <- c("A1intercept",'A2prevguess','A3preverror','B1intercept','B2follguess','B3follerror')
  # 
  # #v <- abs(v) + abs(.1*group$v_coh)
  # 
  # nsim = 5000 #5000 for simuls
  # preds <- data.frame(matrix(NA,nrow=nsim*6))
  # for(row in 1:length(v)){
  #   sigma = 1 
  #   mu = v[row]   
  #   bound = a[row]
  #   ter = t
  #   
  #   params <- c(mu,bound,ter,0);names(params) <- c('v','a','t','z')
  #   output <- RW_hddmbased(params,samples=nsim,intra_sv=sigma,t2time=0,
  #                          evidence_bias = 'no_bias')
  #   
  #   variance = 1
  #   index = ((1*(nsim*(row-1))+1)+(variance-1)*(nsim*length(v))):((1*(nsim*(row)))+(variance-1)*(nsim*length(v)))
  #   preds$acc[index] <- output[,2]
  #   preds$acc[preds$acc==-1] <- 0
  #   preds$rt[index] <- output[,1]
  #   preds$condition[index] <- labels[row]  
  # }
  
  # predsCor <- subset(preds,acc==1)
  # RTpreds <- with(predsCor,aggregate(rt,by=list(condition),median))
  # names(RTpreds) <- c("condition","RT")
  # RTpreds <- cast(RTpreds,~condition)
  # 
  # ACCpreds <- with(preds,aggregate(acc,by=list(condition),mean))
  # names(ACCpreds) <- c("condition","acc")
  # ACCpreds <- cast(ACCpreds,~condition)
  
  #Or load HDDM predictions
  preds <- read.table('simuls_varConfSlowing.csv',sep=',',header=T)
  preds$cor <- 1;preds$cor[preds$rt_sampled>0] <- 0
  preds$rt <- abs(preds$rt_sampled)
  predsCor <- subset(preds,cor==1)

  # meanerr_dat <- with(Data,aggregate(cor,by=list(sub),mean))
  # meanerr_pred <- with(preds,aggregate(cor,by=list(sub),mean))
  # plot(meanerr_dat$x,meanerr_pred$x);cor.test(meanerr_dat$x,meanerr_pred$x)
  # meanrt_dat <- with(Data_RT,aggregate(rt,by=list(sub),median))
  # meanrt_pred <- with(predsCor,aggregate(rt,by=list(sub),median))
  # plot(meanrt_dat$x,meanrt_pred$x);cor.test(meanrt_dat$x,meanrt_pred$x)

  prevcjrt_sim <- with(predsCor,aggregate(rt,by=list(prevconf,node),median))
  prevcjrt_sim <- cast(prevcjrt_sim,Group.2~Group.1)
  postcjrt_sim <- with(predsCor,aggregate(rt,by=list(follconf,node),median))
  postcjrt_sim <- cast(postcjrt_sim,Group.2~Group.1)
  
  prevcjacc_sim <- with(preds,aggregate(cor,by=list(prevconf,node),mean))
  prevcjacc_sim <- cast(prevcjacc_sim,Group.2~Group.1)
  postcjacc_sim <- with(preds,aggregate(cor,by=list(follconf,node),mean))
  postcjacc_sim <- cast(postcjacc_sim,Group.2~Group.1)
}


#COMPARE DISTRIBUTIONS
#1. Distribution
tiff(file="figureswcond/distribution_preds.tiff", res=500, width=500*(350/72), height=500*(350/72))
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
tempC <- hist(Data$rt[Data$cor==1]*1000,breaks=seq(100,3000,100),xlim=c(0,3000),prob=F,col=rgb(0,1,0,.5),border="white",ylab="",xlab="Reaction times (ms)",cex.lab=2, cex.main=3, cex.axis=1.5,main="",yaxt='n');axis(2,labels=F);mtext("Frequency",side=2,line=1,cex=2)
tempE <- hist(Data$rt[Data$cor==0]*1000,breaks=seq(100,3000,100),prob=F,add=T,col=rgb(1,0,0,.5),border='white')
Cors <- hist(preds$rt[preds$cor==1],breaks=seq(.1,19,.1),plot=F)
Errs <- hist(abs(preds$rt[preds$cor==0]),breaks=seq(.15,21,.1),plot=F)
xNew <- Cors$mids*1000;lines(rescale(Cors$counts,to=c(0,max(tempC$counts)))~xNew,type='l',col='green',lwd=3)
xNew <- Errs$mids*1000;lines(rescale(Errs$counts,to=c(0,max(tempE$counts)))~xNew,type='l',col='red',lwd=3)
legend("topright",fill=c("white","white","green","red"),border=F,legend=c("Fits (corrects)","Fits (errors)","Data (corrects)","Data (errors)"),col=rep(c("Green","Red"),2),bty='n',lwd=c(3,3,-1,-1),cex=1)
dev.off()

#RTs & ACCs
par(mfrow=c(1,1))
#Previous Cj
tiff(file="figureswcond/behav_prevcjRT.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- prevcjrt*1000
x_pred <- prevcjrt_sim[2:4]*1000

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(400,1800), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab=substitute(paste("Reaction time on trial ",italic('n+1'))),
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(200+(i*200)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n'))),1,2.5,cex=2)
stripchart(x,at=c(1:3)+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
stripchart(x_pred,at=c(1:3)-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=3,lwd=2)
for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
# lines(means,type='b',lty=2,lwd=3)
# lines(c(1:3)-.2,colMeans(x_pred),type='p',pch=19,cex=2,lwd=3,lty=2,col='blue')
# for(i in 1:3){ temp <- t.test(x_pred[,i]); error.bar(i-.2,mean(x_pred[,i]),temp$conf.int[2]-mean(x_pred[,i]),col='blue',lwd=3)}
dev.off()

tiff(file="figureswcond/behav_prevcjACC.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- prevcjacc*100
x_pred <- prevcjacc_sim[2:4]*100

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(50,100), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab=substitute(paste("Accuracy on trial ",italic('n+1'))),
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(40+(i*10)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n'))),1,2.5,cex=2)
stripchart(x,at=c(1:3)+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
stripchart(x_pred,at=c(1:3)-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=3,lwd=2)
for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
# lines(means,type='b',lty=2,lwd=3)
# lines(c(1:3)-.2,colMeans(x_pred),type='p',pch=19,cex=2,lwd=3,lty=2,col='blue')
# for(i in 1:3){ temp <- t.test(x_pred[,i]); error.bar(i-.2,mean(x_pred[,i]),temp$conf.int[2]-mean(x_pred[,i]),col='blue',lwd=3)}
dev.off()


#Following Cj
tiff(file="figureswcond/behav_postcjRT.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- postcjrt*1000
x_pred <- postcjrt_sim[2:4]*1000

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(400,1800), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab=substitute(paste("Reaction time on trial ",italic('n+1'))),
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(200+(i*200)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n+2'))),1,2.5,cex=2)
stripchart(x,at=c(1:3)+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
stripchart(x_pred,at=c(1:3)-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=3,lwd=2)
for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
# lines(means,type='b',lty=2,lwd=3)
# lines(c(1:3)-.2,colMeans(x_pred),type='p',pch=19,cex=2,lwd=3,lty=2,col='blue')
# for(i in 1:3){ temp <- t.test(x_pred[,i]); error.bar(i-.2,mean(x_pred[,i]),temp$conf.int[2]-mean(x_pred[,i]),col='blue',lwd=3)}
dev.off()

tiff(file="figureswcond/behav_postcjACC.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- postcjacc*100
x_pred <- postcjacc_sim[2:4]*100

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(50,100), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab=substitute(paste("Accuracy on trial ",italic('n+1'))),
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(40+(i*10)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n+2'))),1,2.5,cex=2)
stripchart(x,at=c(1:3)+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
stripchart(x_pred,at=c(1:3)-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=3,lwd=2)
for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
# lines(means,type='b',lty=2,lwd=3)
# lines(c(1:3)-.2,colMeans(x_pred),type='p',pch=19,cex=2,lwd=3,lty=2,col='blue')
# for(i in 1:3){ temp <- t.test(x_pred[,i]); error.bar(i-.2,mean(x_pred[,i]),temp$conf.int[2]-mean(x_pred[,i]),col='blue',lwd=3)}
dev.off()


#DIFFERENCE SCORES
tiff(file="figureswcond/behav_diffRT.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- (prevcjrt-postcjrt)*1000
x_pred <- (prevcjrt_sim[2:4]-postcjrt_sim[2:4])*1000

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(-400,1000), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab='Subsequent reaction time (ms)',
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(-600+(i*200)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence")),1,2.5,cex=2)
stripchart(x,at=c(1:3)+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
stripchart(x_pred,at=c(1:3)-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=3,lwd=2)
for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
#lines(1:3+.15,means,type='b',lty=2,lwd=3)
legend(.75,900,legend=c("Data","Model fits"),col=c("grey","darkgrey"),box.lty=0,pch=19,cex=2)
dev.off()

tiff(file="figureswcond/behav_diffACC.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- (prevcjacc-postcjacc)*100
x_pred <- (prevcjacc_sim[,2:4]-postcjacc_sim[,2:4])*100

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(-20,40), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab='Subsequent accuracy (%)',
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(-30+(i*10)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence")),1,2.5,cex=2)
stripchart(x,at=c(1:3)+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
stripchart(x_pred,at=c(1:3)-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=3,lwd=2)
for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
#for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
# lines(c(1:3)-.2,colMeans(x_pred),type='p',pch=19,cex=2,lwd=3,lty=2,col='blue')
# for(i in 1:3){ temp <- t.test(x_pred[,i]); error.bar(i-.2,mean(x_pred[,i]),temp$conf.int[2]-mean(x_pred[,i]),col='blue',lwd=3)}
dev.off()

tiff(file="figureswcond/behav_diffRT_and_diffAcc.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 4.1, 4.1, 4.1))
x <- (prevcjrt-postcjrt)*1000
x_pred <- (prevcjrt_sim[,2:4]-postcjrt_sim[,2:4])*1000
means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(-800,800), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,axes=F,ylab='',
           cex.axis=1.5,cex.lab=2)
axis(2,at=seq(-400,1200,200),cex.axis=1.5);mtext("Subsequent RT (ms)",2,line=2.5,cex=2,at=50)
axis(4,at=seq(-800,-550,length.out=3),labels=c(-30,0,30),cex.axis=1.5);mtext("Subsequent acc. (%)",4,line=2.5,cex=2,at=-625)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(-600+(i*200)),col='lightgrey',lwd=0.8)  
for(i in 1:3) abline(h=(-800-125+(i*(125))),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence")),1,2.5,cex=2)
#add rts
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=ntrialsCor[,i],lwd=2)
lines(means,type='b',lty=2,lwd=3,cex=3,pch=19)
#for(i in 1:n) stripchart(x_pred[,i],at=i-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=ntrialsCor[,i],lwd=2)
#for(i in 1:n) lines(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2,type='b')
for(i in 1:n) lines(i,sapply(x_pred, mean)[i],type='p',pch=4,col='green',lwd=3,cex=2.5)
#add accuracy
x <- prevcjacc*100-postcjacc*100;x <- rescale(x,to=c(-800,-550),from=c(-30,30))
x_pred <- (prevcjacc_sim[,2:4]-postcjacc_sim[,2:4])*100;x_pred <- rescale(x_pred,to=c(-800,-550),from=c(-30,30))
means <- sapply(x, mean);n<- length(x)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=ntrialsCor[,i],lwd=2)
lines(means,type='b',lty=2,lwd=3,cex=3,pch=19)
#for(i in 1:n) stripchart(x_pred[,i],at=i-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=ntrialsCor[,i],lwd=2)
#for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
for(i in 1:n) lines(i,sapply(x_pred, mean)[i],type='p',pch=4,col='green',lwd=3,cex=2.5)
dev.off()


tiff(file="figureswcond/behav_diffRTbyACC.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
prevIES <- (prevcjrt*1000)*(prevcjacc)
postIES <- (postcjrt*1000)*(postcjacc)
x <- prevIES-postIES
prevIES_sim <- (prevcjrt_sim[,2:4]*1000)*(prevcjacc_sim[,2:4])
postIES_sim <- (postcjrt_sim[,2:4]*1000)*(postcjacc_sim[,2:4])
x_pred <- prevIES_sim-postIES_sim

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(-500,600), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab='Subsequent RT * accuracy',
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(-1200+(i*200)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence")),1,2.5,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=ntrialsCor[,i],lwd=2)
lines(means,type='b',lty=2,lwd=3,cex=3,pch=19)
for(i in 1:n) lines(i,sapply(x_pred, mean)[i],type='p',pch=4,col='green',lwd=3,cex=2.5)
#for(i in 1:n) stripchart(x_pred[,i],at=i-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="lightgreen",col='white',cex=ntrialsCor[,i],lwd=2)
#for(i in 1:n) points(i-.15,sapply(x_pred, mean)[i],pch=21,bg="darkgreen",col="white",cex=3.5,lwd=2)
#for(i in 1:n) points(i+.15,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
#lines(1:3+.15,means,type='b',lty=2,lwd=3)
#legend(.75,900,legend=c("Data","Model fits"),col=c("grey","darkgrey"),box.lty=0,pch=19,cex=2)
dev.off()

tiff(file="figureswcond/behav_RTbyACC_postcj.tiff", res=500, width=800*(300/72), height=700*(300/72))
par(mar=c(5.1, 5.1, 4.1, 4.1))
prevIES <- (prevcjrt*1000)*(prevcjacc)
postIES <- (postcjrt*1000)*(postcjacc)
#x <- prevIES
x <- postIES
prevIES_sim <- (prevcjrt_sim[,2:4]*1000)*(prevcjacc_sim[,2:4])
postIES_sim <- (postcjrt_sim[,2:4]*1000)*(postcjacc_sim[,2:4])
#x_pred <- prevIES_sim
x_pred <- postIES_sim

means <- sapply(x, mean);n<- length(x)
stripchart(x, ylim=c(300,1500), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3,xaxt='n',ylab=substitute(paste("RT * accuracy on trial ",italic(' n+1'))),
           cex.axis=1.5,cex.lab=2)
axis(1,at=1:3,labels=c("High","Low","Perceived error"), cex.axis=1.5)
for(i in 1:20) abline(h=(-1200+(i*200)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n+2'))),1,2.5,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=ntrialsCor[,i],lwd=2)
lines(means,type='b',lty=2,lwd=3,cex=3,pch=19)
for(i in 1:n) lines(i,sapply(x_pred, mean)[i],type='p',pch=4,col='green',lwd=3,cex=2.5)
dev.off()

#analysis 
#RTs
library(multcomp)
prevcj[,3] <- prevcj[,3]-follcj[,3]
prevcj$prevcj <- as.factor(prevcj$prevcj)
fit <- lmer(RT~prevcj+(1|Subject),data=prevcj)
anova(fit)
plot(effect('prevcj',fit))

contrast.matrix <- rbind("cor:guess"  = c(0,1,0),
                         "cor:error"  = c(0,0,1),
                         "guess:error"  = c(0,1,-1)
)
summary(glht(fit, linfct=contrast.matrix), test=adjusted("none"))

#ACC
prevcjacc1[,3] <- prevcjacc1[,3]-follcjacc1[,3]
prevcjacc1$prevcj <- as.factor(prevcjacc1$prevcj)
fit <- lmer(acc~prevcj+(1|Subject),data=prevcjacc1)
anova(fit)
plot(effect('prevcj',fit))

contrast.matrix <- rbind("cor:guess"  = c(0,1,0),
                         "cor:error"  = c(0,0,1),
                         "guess:error"  = c(0,1,-1)
)
summary(glht(fit, linfct=contrast.matrix), test=adjusted("none"))


#how about IES (rt/acc; this is a measure of drift)
prevIES <- prevcjrt/prevcjacc
postIES <- postcjrt/postcjacc
diffIES <- prevIES-postIES
ebar <- barplot(colMeans(diffIES)); error.bar(ebar,colMeans(diffIES),colSds(as.matrix(diffIES))/sqrt(N))
diffIES_long <- melt(diffIES)
diffIES_long$sub <- rep(1:N,3)
fit <- lmer(value~variable+(1|sub),diffIES_long)
anova(fit)


#how about a combined measure of bound (rt*cor)
prevIES <- prevcjrt*prevcjacc
postIES <- postcjrt*postcjacc
diffIES <- prevIES-postIES
ebar <- barplot(colMeans(diffIES)); error.bar(ebar,colMeans(diffIES),colSds(as.matrix(diffIES))/sqrt(N))
diffIES_long <- melt(diffIES)
diffIES_long$sub <- rep(1:N,3)
fit <- lmer(value~variable+(1|sub),diffIES_long)
anova(fit)

contrast.matrix <- rbind("high vs low"= c(0, 1, 0),
                         "high vs err"= c(0, 0, 1),
                         "low vss err"= c(0, 1, -1)
)
summary(glht(fit, linfct=contrast.matrix), test=adjusted("none"))


#how about preselection
whattotest <- "C_error"
if(whattotest=="C_error"){
  Data$include <- 0
  Data$include[Data$prevcj=='A_correct'&Data$follcj=='C_error'] <- 1
  Data$include[Data$prevcj=='C_error'] <- 2
}else{
  Data$include <- 0
  Data$include[Data$prevcj=='A_correct'&Data$follcj=='B_guess'] <- 1
  Data$include[Data$prevcj=='B_guess'] <- 2
}
table(Data$include,Data$prevcj)
temp <- subset(Data,include>0);tempCor <- subset(temp,cor==1)

prevcjrt <- with(tempCor,aggregate(rt,by=list(sub,prevcj),median))
prevcjrt <- cast(prevcjrt,Group.1~Group.2)
prevcjacc <- with(temp,aggregate(cor,by=list(sub,prevcj),mean))
prevcjacc <- cast(prevcjacc,Group.1~Group.2)
t.test(prevcjrt[,2],prevcjrt[,3],paired=T);colMeans(prevcjrt)
t.test(prevcjacc[,2],prevcjacc[,3],paired=T);colMeans(prevcjacc)
t.test(prevcjrt[,2]*prevcjacc[,2],prevcjrt[,3]*prevcjacc[,3],paired=T)



###############################################
#How do confidence judgments relate to RT/acc #
par(mfrow=c(1,1))
plot(density(Data$prevrt[Data$prevcj=='A_correct'&Data$prevcor==1]),xlab='RTs',ylab='density',main="",xlim=c(-3,3),frame=F)
polygon(density(Data$prevrt[Data$prevcj=='A_correct'&Data$prevcor==1]),col=rgb(0,0,1,.5))
polygon(density(Data$prevrt[Data$prevcj=='B_guess'&Data$prevcor==1]),col=rgb(1,0,0,.5))
polygon(density(Data$prevrt[Data$prevcj=='C_error'&Data$prevcor==1]),col=rgb(0,1,0,.5))
polygon(density(-Data$prevrt[Data$prevcj=='A_correct'&Data$prevcor==0]),col=rgb(0,0,1,.5))
polygon(density(-Data$prevrt[Data$prevcj=='B_guess'&Data$prevcor==0]),col=rgb(1,0,0,.5))
polygon(density(-Data$prevrt[Data$prevcj=='C_error'&Data$prevcor==0]),col=rgb(0,1,0,.5))




###########################################################???
### Resision 1: REMOVE LOW-FREQUENCY OSCILLATIONS
{
  library(signal)
  Data$rt_filtered <- NA
  for(i in 1:N){
    temp <- subset(Data,sub==subs[i])
    # bf <- butter(5, 1/10, type="high")
    # Data$rt_filtered[Data$sub==subs[i]] <- filter(bf, temp$rt)
  }
  plot(Data$rt,Data$rt_filtered)

  hist(Data$rt,breaks=50,col='darkgrey',xlim=c(-1,2))
  hist(Data$rt_filtered,breaks=50,col='lightgrey',add=T)
  legend("topright",inset=.2, legend=c("Raw","Filtered"),col=c('darkgrey','lightgrey'),lwd=3,cex=1.5,box.lty=0)

  DataCor <- subset(Data,cor==1)
  DataCor$prevcj <- as.factor(DataCor$prevcj)
  fit <- lmer(rt~prevcj+(1|sub),data=DataCor)
  fit <- lmer(rt_filtered~prevcj+(1|sub),data=DataCor)
  anova(fit)
  plot(effect('prevcj',fit))



  ##############################################################
  #Revision 1: INSPECT FREQUENCY SPECTRA OF DATA
  #Check blocks on RT and accuracy
  library(bspec);
  power_over_subs <- matrix(NA,nrow=N,ncol=26);random_power_over_subs <- matrix(NA,nrow=N,ncol=26);residuals_power_over_subs <- matrix(NA,N,26)
  #pdf(file="slowdrift/slowdrift exp1.pdf")
  for(i in 1:N){
    temp <- subset(Data,sub==subs[i])
    Bs <- unique(temp$block)

    #1. Show raw data
    # par(mfrow=c(3,2))
    # for(b in 1:length(Bs)){
    #   tempB <- subset(temp,block==Bs[b])
    #   tempB$trial <- 1:length(tempB$sub)
    #   #tempB$rt <- residuals(lm(tempB$rt~tempB$coh))+lm(tempB$rt~tempB$coh)$coefficients[1]
    #   #tempB$cor <- round(residuals(lm(tempB$cor~tempB$coh))+lm(tempB$cor~tempB$coh)$coefficients[1])
    # 
    #   if(b==1){lab=paste("sub",i)} else{lab=""}
    #   plot(tempB$rt,pch=19,frame=F,ylab='RT',ylim=c(-1.4,3.3),main=lab,yaxt='n');axis(2,at=0:3);axis(2,at=c(-.2,-1.2),labels=c('cor','err'),las=2)
    #   abline(lm(tempB$rt~tempB$trial),lty=2,col="blue")
    #   points(tempB$cor-1.2,pch=4)
    #   abline(lm(tempB$cor-1.2~tempB$trial),lty=2,col="blue")
    #   text(30,3.2,paste("RT: p = ", round(anova(lm(tempB$rt~tempB$trial))[1,5],3)),col='blue',cex=.6)
    #   text(30,3,paste("ACC: p = ", round(anova(lm(tempB$cor~tempB$trial))[1,5],3)),col='blue',cex=.6)
    # }

    #2. RT Distributions
    # par(mfrow=c(2,2))
    # tempC <- hist(temp$rt[temp$cor==1],breaks=seq(0,3,.2),col=rgb(0,1,0,.5),xlim=c(0,3),main=paste('sub',i))
    # tempE <- hist(temp$rt[temp$cor==0],breaks=seq(0,3,.2),col=rgb(1,0,0,.5),add=T)
    # Cors <- hist(temp$rt[temp$cor==1],breaks=seq(0,3,.2),plot=F);Errs <- hist(temp$rt[temp$cor==0],breaks=seq(0,3,.2),plot=F)
    # lines(rescale(Cors$counts,to=c(0,max(tempC$counts)))~Cors$mids,type='l',col='green',lwd=3)
    # lines(rescale(Errs$counts,to=c(0,max(tempE$counts)))~Errs$mids,type='l',col='red',lwd=3)

    #3. Power spectral density estimation using Welch's
    #3.1. Trial as unit of measurement
    #3.1.2 Ignore blocks
    outp = welchPSD(x=as.ts(temp$conf),seglength=50)
    plot(outp$frequency[-1],outp$power[-1],log="xy",xlab="Frequency (Hz)",ylab="Power",frame=F,type='b', main="Unit = trial",lwd=2,pch=19)
    power_over_subs[i,] = outp$power
    random_power_over_subs[i,] = welchPSD(x=as.ts(sample(temp$rt)),seglength=50)$power
    #residuals_power_over_subs[i,] = welchPSD(x=as.ts(residuals(lm(temp$rt~temp$follcj))),seglength=50)$power

    #3.1.2. Per block and then average
    # freqS <- matrix(NA,ncol=length(Bs),nrow=13); Powers <- matrix(NA,ncol=length(Bs),nrow=13)
    # for(b in 1:length(Bs)){
    #   tempB <- subset(temp,block==Bs[b])
    #   outp = welchPSD(x=as.ts(tempB$rt),seglength=25)
    #   Powers[,b] <- outp$power
    # }
    # plot(outp$frequency[-1],rowMeans(Powers)[-1],log="xy",xlab="Frequency (Hz)",ylab="Power",frame=F,main="Unit = trial")
    # for(b in 1:length(Bs)){lines(outp$frequency[-1],Powers[-1,b],col="grey",type='b')};lines(outp$frequency,rowMeans(Powers),type='b',lwd=2,pch=19)
    # legend("bottom",inset=.2, legend="blocks",col='grey',lwd=3,cex=1.5,box.lty=0)
  }
  #dev.off()

  #Group Average
  #jpeg(file="slowdrift/groupAverage exp1.jpeg")
  par(mfrow=c(1,1))
  plot(outp$frequency[-1],colMeans(power_over_subs)[-1],log="xy",xlab="Frequency",ylab="Power",frame=F,main="Measurement unit = trial")
  for(b in 1:N){lines(outp$frequency[-1],power_over_subs[b,-1],col="grey",type='b')};lines(outp$frequency[-1],colMeans(power_over_subs)[-1],lwd=2,type='b',pch=19)
  legend("bottom",inset=.2, legend=c("participants","mean","regress out post-conf"),col=c('grey','black','red'),lwd=3,cex=1.5,box.lty=0)
  #lines(outp$frequency[-1],colMeans(random_power_over_subs)[-1],lwd=2,type='b',pch=19,col='green')
  lines(outp$frequency[-1],colMeans(residuals_power_over_subs)[-1],lwd=2,type='b',pch=19,col='red')
  #dev.off()

  #Linear fit to log-log spectra
  plot(log(outp$frequency[-1]),log(colMeans(power_over_subs)[-1]))
  fit <- lm(log(colMeans(power_over_subs)[-1])~log(outp$frequency[-1]))
  summary(fit)
  abline(fit)
  
  #stripchart differences
  stripchart(power_over_subs[,]-random_power_over_subs[,],ylim=c(-2,2),xlim=c(0,.5), vertical = TRUE,frame=F,col='white')
  abline(h=0,lty=2)
  for(i in 2:26){ stripchart(power_over_subs[,i]-random_power_over_subs[,i],at=outp$frequency[i],method = 'jitter',add=T,jitter=.005,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)}
  lines(outp$frequency[-1],colMeans(power_over_subs)[-1]-colMeans(random_power_over_subs)[-1],pch=19,type='b',lwd=3)
  for(i in 1:26){
    if((t.test(power_over_subs[,i]-random_power_over_subs[,i])$p.value)<.05){
      points(y=-1.5,x=outp$frequency[i],pch=4)
    }
  }


  #FREQ SPECTRA TAKING SPECIFIC MOMENT INTO ACCOUNT
  #Assume a sampling rate of 5Hz (5 measures per second)
  library(bspec);
  power_over_subs <- matrix(NA,nrow=N,ncol=26)
  subs = with(Data,aggregate(sub,by=list(sub),mean))$x
  for(i in 1:N){
    temp <- subset(Data,sub==subs[i])
    temp$rt_in_time <- (temp$rt_in_time-temp$rt_in_time[1])+1 #start at 1

    rt_timeseries <- rep(NA,ceiling(max(temp$rt_in_time))*5) #in steps of 200ms
    index = 1;whichIDs <- NA
    for(t in 1:length(rt_timeseries)){
      if( round(temp$rt_in_time[index],1) >= t/5 & round(temp$rt_in_time[index],1) < (t+1)/5){
        rt_timeseries[t] <- temp$rt[index]
        whichIDs[index] <- t
        index=index+1
        if(index>dim(temp)[1]){break}
      }
    }

    #linear interpolation of missing values
    for(t in 1:(length(whichIDs)-1)){
      rt_timeseries[(whichIDs[t]+1):(whichIDs[t+1]-1)] <- seq(rt_timeseries[whichIDs[t]], rt_timeseries[whichIDs[t+1]],length.out=length(rt_timeseries[(whichIDs[t]+1):(whichIDs[t+1]-1)]))
    }
    rt_timeseries = rt_timeseries[whichIDs[1]:whichIDs[length(whichIDs)]] #cut off the edges

    #Power spectral density estimation using Welch's
    outp = welchPSD(x=as.ts(rt_timeseries),seglength=50)
    plot(outp$frequency[-1],outp$power[-1],log="xy",xlab="Frequency",ylab="Power",frame=F,type='b', main="Unit = trial",lwd=2,pch=19)
    power_over_subs[i,] = outp$power
  }

  #Group Average
  #jpeg(file="slowdrift/groupAverage exp1.jpeg")
  par(mfrow=c(1,1))
  plot(outp$frequency[-1],colMeans(power_over_subs)[-1],log="xy",xlab="Frequency",ylab="Power",frame=F,main="Sampling rate = 5Hz")
  for(b in 1:N){lines(outp$frequency[-1],power_over_subs[b,-1],col="grey",type='b')};lines(outp$frequency[-1],colMeans(power_over_subs)[-1],lwd=2,type='b',pch=19)
  #dev.off()
}

