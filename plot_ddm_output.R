#This file uses R to load in the results of the hddm fit and creates nice figures for them
#
setwd("...")
#so visuals to be used later on
rm(list=ls());col1 <- rgb(106/255,106/255,106/255,1)
cols=c('grey',rgb(1,0,0,.2),rgb(0,0,1,.2))
library(scales)

#Load the hddm output
allData = read.csv("Experiment1_HDDMestimates.csv",header=T,sep=',')

#make clearer dataframes
bound3 <- data.frame(matrix(0,28,1));names(bound3) <- c("baseline")
drift3 <- data.frame(matrix(0,28,1));names(drift3) <- c("baseline")
bound3$prevcj_guess <- allData$mean[grep("a_prevconfT.B_guess_subj.",allData$X)]
bound3$prevcj_error <- allData$mean[grep("a_prevconfT.C_error_subj.",allData$X)]
bound3$follcj_guess <- allData$mean[grep("a_follconfT.B_guess_subj.",allData$X)]
bound3$follcj_error <- allData$mean[grep("a_follconfT.C_error_subj.",allData$X)]
drift3$prevcj_guess <- allData$mean[grep("v_prevconfT.B_guess_subj.",allData$X)]
drift3$prevcj_error <- allData$mean[grep("v_prevconfT.C_error_subj.",allData$X)]
drift3$follcj_guess <- allData$mean[grep("v_follconfT.B_guess_subj.",allData$X)]
drift3$follcj_error <- allData$mean[grep("v_follconfT.C_error_subj.",allData$X)]
#control for following cj
bound3$guess_diff <- bound3$prevcj_guess-bound3$follcj_guess
bound3$error_diff <- bound3$prevcj_error-bound3$follcj_error
drift3$guess_diff <- drift3$prevcj_guess-drift3$follcj_guess
drift3$error_diff <- drift3$prevcj_error-drift3$follcj_error

#extract group level stats
group <- data.frame(matrix(0,1,1));
group$t <- allData$mean[1]
group$a_intercept <- allData$mean[31]
group$a_prevguess <- allData$mean[61]
group$a_preverror <- allData$mean[91]
group$a_follguess <- allData$mean[121]
group$a_follerror <- allData$mean[151]
group$v_intercept <- allData$mean[211]
group$v_prevguess <- allData$mean[241]
group$v_preverror <- allData$mean[271]
group$v_follguess <- allData$mean[301]
group$v_follerror <- allData$mean[331]
group$v_coh <- allData$mean[391]

#load the individual traces
a_guess_prevcj <- read.csv("hddm/a_guess_prevcj.csv",header=F)$V1
a_error_prevcj <- read.csv("hddm/a_error_prevcj.csv",header=F)$V1
v_guess_prevcj <- read.csv("hddm/v_guess_prevcj.csv",header=F)$V1
v_error_prevcj <- read.csv("hddm/v_error_prevcj.csv",header=F)$V1
a_guess_postcj <- read.csv("hddm/a_guess_postcj.csv",header=F)$V1
a_error_postcj <- read.csv("hddm/a_error_postcj.csv",header=F)$V1
v_guess_postcj <- read.csv("hddm/v_guess_postcj.csv",header=F)$V1
v_error_postcj <- read.csv("hddm/v_error_postcj.csv",header=F)$V1
#control for following cj
a_guess_diff <- a_guess_prevcj-a_guess_postcj
a_error_diff <- a_error_prevcj-a_error_postcj
v_guess_diff <- v_guess_prevcj-v_guess_postcj
v_error_diff <- v_error_prevcj-v_error_postcj


#1.Plot bound results
#raw effect previous cj; Figure 3 - Figure Supplement 3
cols=c('grey',rgb(1,0,0,.2),rgb(0,0,1,.2))
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- bound3[,c("baseline","prevcj_guess","prevcj_error")]
names(x) <- c("High","Low","Perceived error")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.2,.6),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.axis=1.5,cex.lab=2,ylab=substitute(paste("Decision bound on trial ",italic('n+1'))),
           main="")
for(i in 1:20) abline(h=(-0.6+(i*.1)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n'))),1,2.8,cex=2)
#stripchart(x,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey,col='white',cex=3,lwd=2)
polygon(rescale(density(a_guess_prevcj)$y,to=c(-.1,.8),from=c(0,10)),density(a_guess_prevcj)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(a_error_prevcj)$y,to=c(-.1,.8),from=c(0,10)),density(a_error_prevcj)$x,col=rgb(0,0,1,.2),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='b',lty=2,lwd=3)

#raw effect following cj; Figure 3 - Figure Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- bound3[,c("baseline","follcj_guess","follcj_error")]
names(x) <- c("High","Low","Perceived error")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.2,.6), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.axis=1.5,cex.lab=2,ylab=substitute(paste("Decision bound on trial ",italic('n+1'))),
           main="")
for(i in 1:40) abline(h=(-0.6+(i*.1)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n+2'))),1,2.8,cex=2)
polygon(rescale(density(a_guess_postcj)$y,to=c(-.1,.8),from=c(0,16)),density(a_guess_postcj)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(a_error_postcj)$y,to=c(-.1,.8),from=c(0,16)),density(a_error_postcj)$x,col=rgb(0,0,1,.2),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='b',lty=2,lwd=3)

#key finding; Figure 3B, left
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- bound3[,c("baseline","guess_diff","error_diff")]
names(x) <- c("High","Low","Perceived error")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.25,.6), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.axis=1.5,cex.lab=2,ylab="Subsequent decision bound",
           main="")
for(i in 1:20) abline(h=(-0.6+(i*.1)),col='lightgrey',lwd=0.8)  
mtext("Confidence (correct trials)",1,2.8,cex=2,at=2)
#lines(x=c(-2,.8),y=c(0,0),lty=2,lwd=2)
polygon(rescale(density(a_guess_diff)$y,to=c(-.1,.8),from=c(0,8)),density(a_guess_diff)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(a_error_diff)$y,to=c(-.1,.8),from=c(0,8)),density(a_error_diff)$x,col=rgb(0,0,1,.2),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='b',lty=2,lwd=3)


#Plot DRIFT findings 
#raw effect previous cj; Figure 3 - Figure Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- drift3[,c("baseline","prevcj_guess","prevcj_error")]
names(x) <- c("High","Low","Perceived error")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.2,.3), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.axis=1.5,cex.lab=2,ylab=substitute(paste("Drift rate on trial",italic('n+1'))),
           main="")
for(i in 1:12) abline(h=(-0.8+(i*.1)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n'))),1,2.8,cex=2)
polygon(rescale(density(v_guess_prevcj)$y,to=c(-.1,.8),from=c(0,13)),density(v_guess_prevcj)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v_error_prevcj)$y,to=c(-.1,.8),from=c(0,13)),density(v_error_prevcj)$x,col=rgb(0,0,1,.2),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='b',lty=2,lwd=3)

#raw effect following cj; Figure 3 - Figure Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- drift3[,c("baseline","follcj_guess","follcj_error")]
names(x) <- c("High","Low","Perceived error")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.2,.3), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.axis=1.5,cex.lab=2,ylab=substitute(paste("Drift rate on trial",italic('n+1'))),
           main="")
for(i in 1:16) abline(h=(-1+(i*.1)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Confidence on trial ",italic('n+2'))),1,2.8,cex=2)
polygon(rescale(density(v_guess_postcj)$y,to=c(-.1,.8),from=c(0,12)),density(v_guess_postcj)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v_error_postcj)$y,to=c(-.1,.8),from=c(0,12)),density(v_error_postcj)$x,col=rgb(0,0,1,.2),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='b',lty=2,lwd=3)

#main finding, Figure 3B
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- drift3[,c("baseline","guess_diff","error_diff")]
names(x) <- c("High","Low","Perceived error")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.2,.3), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.axis=1.5,cex.lab=2,ylab="Subsequent drift rate",
           main="")
for(i in 1:30) abline(h=(-.5+(i*.1)),col='lightgrey',lwd=0.8)  
mtext("Confidence (correct trials)",1,2.8,cex=2,at=2)
polygon(rescale(density(v_guess_diff)$y,to=c(-.1,.8),from=c(0,9)),density(v_guess_diff)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v_error_diff)$y,to=c(-.1,.8),from=c(0,9)),density(v_error_diff)$x,col=rgb(0,0,1,.2),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='b',lty=2,lwd=3)



##########################################################
#1.1 Model with only coherence (cf. Figure 2B)
SNR1 = read.csv("varConf_CohOnly.csv",header=T,sep=','); snr1 <- data.frame(matrix(NA,11,1))

snr1 <- data.frame(NA,1:28)
snr1$sub <- SNR1[c(37:64),1]
snr1$v0 <- SNR1[c(37:64),2]
snr1$v5 <- SNR1[c(65:92),2]
snr1$v10 <- SNR1[c(93:120),2]
snr1$v20 <- SNR1[c(121:148),2]
snr1$v40 <- SNR1[c(149:176),2]
snr1$a <- SNR1[c(3:30),2]
snr1$t <- SNR1[c(179:206),2]

#Query predictions from the DDM, using a custom made function rw_hddmbased.r
source('rw_hddmbased.R')

#contruct matrices with v and a for each level of confidence prev and foll
a <- snr1$a
v <- snr1[,4:8]
t <- snr1$t
labels=c(0,5,10,20,40)

nsim = 5000 #5000 for simuls
preds <- data.frame(matrix(NA,nrow=nsim*length(v)*length(table(Data$sub))))
for(sub in 1:length(table(Data$sub))){
  for(row in 1:length(v)){
    sigma = 1 
    mu = v[sub,row]   
    bound = a[sub]
    ter = t[sub]
    
    params <- c(mu,bound,ter,0);names(params) <- c('v','a','t','z')
    output <- RW_hddmbased(params,samples=nsim,intra_sv=sigma,t2time=0,
                           evidence_bias = 'no_bias')
    
    variance = sub
    index = ((1*(nsim*(row-1))+1)+(variance-1)*(nsim*length(v))):((1*(nsim*(row)))+(variance-1)*(nsim*length(v)))
    preds$acc[index] <- output[,2]
    preds$rt[index] <- output[,1]
    preds$condition[index] <- labels[row]
    preds$sub[index] <- sub
  }
  print(paste('running simulations for subject',sub,'out of ', length(table(Data$sub))))
}

preds$acc[preds$acc==-1] <- 0
predsCor <- subset(preds,acc==1)
RTpreds <- with(predsCor,aggregate(rt,by=list(condition,sub),median))
names(RTpreds) <- c("condition","sub","RT")
RTpreds <- cast(RTpreds,sub~condition)

ACCpreds <- with(preds,aggregate(acc,by=list(condition,sub),mean))
names(ACCpreds) <- c("condition","sub","acc")
ACCpreds <- cast(ACCpreds,sub~condition)

#cf figure 2B
par(mar=c(5.1, 5.1, 4.1, 4.1))
x <- data.frame(cbind(snr1$v0,snr1$v5,snr1$v10,snr1$v20,snr1$v40))
names(x) <- c(0,5,10,20,40)
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x, ylim=c(-.4,4),ylab="Drift rate",xlim=c(-5,45), vertical = TRUE, col="white",frame=F, xaxt='n',main="Experiment 1",yaxt='n')
axis(1,at=as.numeric(names(x)),labels=names(x))
axis(2,at=c(0:4))
for(i in 1:20) abline(h=(-1+(i*.5)),col='lightgrey')  
mtext("Signal-to-noise ratio",1,2.8)
stripchart(x,at=as.numeric(names(x)), method = 'jitter',jitter = 1.5,add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5)
for(i in 1:n) points(names(x)[i],means[i],pch=21,bg="black",col="white",cex=1.5)
lines(means~as.numeric(names(x)),type='b')



