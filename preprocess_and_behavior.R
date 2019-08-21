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
Data$t2 <- as.factor(Data$t2duration) #t2=0=immediate confidence (excluded here), only t2=1 and t2=2 are relevant 

#exclude nonsense RTs < .2
Data <- subset(Data, rt > .2) 

#Plot raw stuff
N <- length(table(Data$sub))
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


#############################
#SAVE THESE DATA FOR THE DDM
hddmData <- Data
hddmData <- subset(hddmData,t2 != 0) #drop the immediate confidence conditions
hddmData <- subset(hddmData,prevconf > 0) #drop nonsense trials
hddmData <- subset(hddmData,follconf > 0) #drop nonsense trials
hddmData <- subset(hddmData,withinblocktrial > 1) #delete first trials (no prevconf)
hddmData <- subset(hddmData,withinblocktrial < 60) #delete last trials (no postconf)
hddmData$subj_idx <- hddmData$sub #some labels that hddm needs
hddmData$resp <- hddmData$resp-1
hddmData$stimulus <- -1
hddmData$stimulus[hddmData$resp==1 & hddmData$cor == 1] <- 1
hddmData$stimulus[hddmData$resp==0 & hddmData$cor == 0] <- 1
hddmData$response <- hddmData$resp
Data <- hddmData
write.table(hddmData,"Experiment1_hddm.csv",quote=F,row.names=F,sep=',')


##############################################################
### Behavioral analysis
##############################################################
#start with some recoding
Data$coh <- as.factor(Data$coh)
DataCor <- subset(Data,cor==1)
DataErr <- subset(Data,cor==0)

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

#Show that behavior scales with coherence 
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

#RTs and Accuracy; Figure 2B
par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- snrrt[,2:6]
y <- snrcor[,2:6]

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
stripchart(x*1000,at=as.numeric(names(x))+.05,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)

for(i in 1:n) points(as.numeric(names(x))[i]+.05,means[i]*1000,pch=21,bg="black",col="white",cex=1.5)
stripchart(y_rescale,at=as.numeric(names(y))+.05,method = 'jitter',jitter=.01,add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)
meansy <- sapply(y_rescale, mean,na.rm=T);
for(i in 1:n) points(as.numeric(names(y))[i]+.05,meansy[i],pch=21,bg="black",col="white",cex=1.5)
lines(as.numeric(names(x))+.05,means*1000,type='b')
lines(as.numeric(names(y))+.05,meansy,type='b')


#Confidence; Figure 2C
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

#RT predicts accuracy; 
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

par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- accsnrrt[,11:-1:2]*100
y <- rtsnrrt[,11:-1:2]*1000

#Figure 2E
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

#analysis
fit <- lmer(acc~rtquint + (1|sub),data=accsnrrt1)
fit1 <- lmer(acc~rtquint + (rtquint|sub),data=accsnrrt1)
anova(fit,fit1)
summary(fit1)


#Seperate axis for accuracy (per label) and frequency
Data$ones <- rep(1,length(Data$sub))
tempDat <- subset(Data, t2!=0) #To get more specific, change this to 1=postevidence, or 2=blank

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

#Figure 2D
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


#RTs as a function of prevconf
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

#Query predictions from the hddm to overlay on the figure
preds <- read.table('simulations_Experiment1.csv',sep=',',header=T) #these are saved in python
preds$cor <- 1;preds$cor[preds$rt_sampled>0] <- 0
preds$rt <- abs(preds$rt_sampled)
predsCor <- subset(preds,cor==1)

prevcjrt_sim <- with(predsCor,aggregate(rt,by=list(prevconf,node),median))
prevcjrt_sim <- cast(prevcjrt_sim,Group.2~Group.1)
postcjrt_sim <- with(predsCor,aggregate(rt,by=list(follconf,node),median))
postcjrt_sim <- cast(postcjrt_sim,Group.2~Group.1)

prevcjacc_sim <- with(preds,aggregate(cor,by=list(prevconf,node),mean))
prevcjacc_sim <- cast(prevcjacc_sim,Group.2~Group.1)
postcjacc_sim <- with(preds,aggregate(cor,by=list(follconf,node),mean))
postcjacc_sim <- cast(postcjacc_sim,Group.2~Group.1)


#COMPARE DISTRIBUTIONS
#inset Figure 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
tempC <- hist(Data$rt[Data$cor==1]*1000,breaks=seq(100,3000,100),xlim=c(0,3000),prob=F,col=rgb(0,1,0,.5),border="white",ylab="",xlab="Reaction times (ms)",cex.lab=2, cex.main=3, cex.axis=1.5,main="",yaxt='n');axis(2,labels=F);mtext("Frequency",side=2,line=1,cex=2)
tempE <- hist(Data$rt[Data$cor==0]*1000,breaks=seq(100,3000,100),prob=F,add=T,col=rgb(1,0,0,.5),border='white')
Cors <- hist(preds$rt[preds$cor==1],breaks=seq(.1,19,.1),plot=F)
Errs <- hist(abs(preds$rt[preds$cor==0]),breaks=seq(.15,21,.1),plot=F)
xNew <- Cors$mids*1000;lines(rescale(Cors$counts,to=c(0,max(tempC$counts)))~xNew,type='l',col='green',lwd=3)
xNew <- Errs$mids*1000;lines(rescale(Errs$counts,to=c(0,max(tempE$counts)))~xNew,type='l',col='red',lwd=3)
legend("topright",fill=c("white","white","green","red"),border=F,legend=c("Fits (corrects)","Fits (errors)","Data (corrects)","Data (errors)"),col=rep(c("Green","Red"),2),bty='n',lwd=c(3,3,-1,-1),cex=1)

#Figure 3A, left
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
for(i in 1:n) lines(i,sapply(x_pred, mean)[i],type='p',pch=4,col='green',lwd=3,cex=2.5)
#add accuracy
x <- prevcjacc*100-postcjacc*100;x <- rescale(x,to=c(-800,-550),from=c(-30,30))
x_pred <- (prevcjacc_sim[,2:4]-postcjacc_sim[,2:4])*100;x_pred <- rescale(x_pred,to=c(-800,-550),from=c(-30,30))
means <- sapply(x, mean);n<- length(x)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=ntrialsCor[,i],lwd=2)
lines(means,type='b',lty=2,lwd=3,cex=3,pch=19)
for(i in 1:n) lines(i,sapply(x_pred, mean)[i],type='p',pch=4,col='green',lwd=3,cex=2.5)

#Figure 3A, right
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


#analysis 
library(multcomp)
prevcj[,3] <- prevcj[,3]-follcj[,3] #so from now on, prevcj is controlled for follcj
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
