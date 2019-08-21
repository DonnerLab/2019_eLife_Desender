setwd("...")
rm(list=ls())
col1 <- rgb(106/255,106/255,106/255,1)
cols=c('grey',rgb(1,0,0,.2),rgb(0,0,1,.2))
col5 <- colorRampPalette(c(rgb(1,0,0,.5),rgb(0,1,0,.5)))
col5 <- col5(5)
par(cex.lab=1.5,lwd=3,cex=1.5)








###############################
# EEGbinned Pe   
###############################
#Load the data
allPE = read.csv("eeg_binnedPe_HDDMestimates.csv",header=T,sep=',')

bound5 <- data.frame(matrix(0,16,1));names(bound5) <- c("baseline")
drift5 <- data.frame(matrix(0,16,1));names(drift5) <- c("baseline")
bound5$prevpe2 <- allPE[c(39:54),2]
bound5$prevpe3 <- allPE[c(57:72),2]
bound5$prevpe4 <- allPE[c(75:90),2]
bound5$prevpe5 <- allPE[c(93:108),2]
bound5$postpe2 <- allPE[c(111:126),2]
bound5$postpe3 <- allPE[c(129:144),2]
bound5$postpe4 <- allPE[c(147:162),2]
bound5$postpe5 <- allPE[c(165:180),2]

drift5$prevpe2 <- allPE[c((219):(234)),2]
drift5$prevpe3 <- allPE[c((237):(252)),2]
drift5$prevpe4 <- allPE[c((255):(270)),2]
drift5$prevpe5 <- allPE[c((273):(288)),2]
drift5$postpe2 <- allPE[c((291):(306)),2]
drift5$postpe3 <- allPE[c((309):(324)),2]
drift5$postpe4 <- allPE[c((327):(342)),2]
drift5$postpe5 <- allPE[c((345):(360)),2]

#control for following pe bin
bound5$diff2 <- bound5$prevpe2-bound5$postpe2
bound5$diff3 <- bound5$prevpe3-bound5$postpe3
bound5$diff4 <- bound5$prevpe4-bound5$postpe4
bound5$diff5 <- bound5$prevpe5-bound5$postpe5
drift5$diff2 <- drift5$prevpe2-drift5$postpe2
drift5$diff3 <- drift5$prevpe3-drift5$postpe3
drift5$diff4 <- drift5$prevpe4-drift5$postpe4
drift5$diff5 <- drift5$prevpe5-drift5$postpe5

#3.1.1. PE
a_prevcj_2 <- read.csv("PEonly_noern_a_prevcj_bin2.csv",header=F)$V1
a_prevcj_3 <- read.csv("PEonly_noern_a_prevcj_bin3.csv",header=F)$V1
a_prevcj_4 <- read.csv("PEonly_noern_a_prevcj_bin4.csv",header=F)$V1
a_prevcj_5 <- read.csv("PEonly_noern_a_prevcj_bin5.csv",header=F)$V1
v_prevcj_2 <- read.csv("PEonly_noern_v_prevcj_bin2.csv",header=F)$V1
v_prevcj_3 <- read.csv("PEonly_noern_v_prevcj_bin3.csv",header=F)$V1
v_prevcj_4 <- read.csv("PEonly_noern_v_prevcj_bin4.csv",header=F)$V1
v_prevcj_5 <- read.csv("PEonly_noern_v_prevcj_bin5.csv",header=F)$V1

a_postcj_2 <- read.csv("PEonly_noern_a_postcj_bin2.csv",header=F)$V1
a_postcj_3 <- read.csv("PEonly_noern_a_postcj_bin3.csv",header=F)$V1
a_postcj_4 <- read.csv("PEonly_noern_a_postcj_bin4.csv",header=F)$V1
a_postcj_5 <- read.csv("PEonly_noern_a_postcj_bin5.csv",header=F)$V1
v_postcj_2 <- read.csv("PEonly_noern_v_postcj_bin2.csv",header=F)$V1
v_postcj_3 <- read.csv("PEonly_noern_v_postcj_bin3.csv",header=F)$V1
v_postcj_4 <- read.csv("PEonly_noern_v_postcj_bin4.csv",header=F)$V1
v_postcj_5 <- read.csv("PEonly_noern_v_postcj_bin5.csv",header=F)$V1

a_diff_2 <- a_prevcj_2-a_postcj_2
a_diff_3 <- a_prevcj_3-a_postcj_3
a_diff_4 <- a_prevcj_4-a_postcj_4
a_diff_5 <- a_prevcj_5-a_postcj_5
v_diff_2 <- v_prevcj_2-v_postcj_2
v_diff_3 <- v_prevcj_3-v_postcj_3
v_diff_4 <- v_prevcj_4-v_postcj_4
v_diff_5 <- v_prevcj_5-v_postcj_5



#PLOT THIS STUFF
cols = c(rgb(.7,.7,.7,.5),rgb(.5,.6,.5,.5),rgb(.3,.3,.3,.5),rgb(.1,.1,.1,.5),'black')
cols = cols[5:1]
#BOUND
#Previous pe; Figure 8 - Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- bound5[,c("baseline","prevpe2","prevpe3","prevpe4","prevpe5")]
names(x) <- c("Lowest","Low","Medium","High","Highest")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.1,.1), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.lab=2,cex.axis=1.5,ylab=substitute(paste("Decision bound on trial ",italic('n+1'))))
for(i in 1:10) abline(h=(-0.2+(i*.05)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Pe amplitude on trial ",italic('n'))),1,2.8,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
lines(means,type='l',lty=2,lwd=3)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
polygon(rescale(density(a_prevcj_2)$y,to=c(-.1,.8),from=c(0,26)),density(a_prevcj_2)$x,col=cols[2],border=F)
polygon(rescale(density(a_prevcj_3)$y,to=c(-.1,.8),from=c(0,26)),density(a_prevcj_3)$x,col=cols[3],border=F)
polygon(rescale(density(a_prevcj_4)$y,to=c(-.1,.8),from=c(0,26)),density(a_prevcj_4)$x,col=cols[4],border=F)
polygon(rescale(density(a_prevcj_5)$y,to=c(-.1,.8),from=c(0,26)),density(a_prevcj_5)$x,col=cols[5],border=F)

#Following cj,; Figure 8 - Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- bound5[,c("baseline","postpe2","postpe3","postpe4","postpe5")]
names(x) <- c("Lowest","Low","Medium","High","Highest")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.7,.1), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.lab=2, cex.axis=1.5,ylab=substitute(paste("Decision bound on trial ",italic('n+1'))))
for(i in 1:10) abline(h=(-0.8+(i*.2)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Pe amplitude on trial ",italic('n+2'))),1,2.8,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
lines(means,type='l',lty=2,lwd=3)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
polygon(rescale(density(a_postcj_2)$y,to=c(-.1,.8),from=c(0,19)),density(a_postcj_2)$x,col=cols[2],border=F)
polygon(rescale(density(a_postcj_3)$y,to=c(-.1,.8),from=c(0,19)),density(a_postcj_3)$x,col=cols[3],border=F)
polygon(rescale(density(a_postcj_4)$y,to=c(-.1,.8),from=c(0,19)),density(a_postcj_4)$x,col=cols[4],border=F)
polygon(rescale(density(a_postcj_5)$y,to=c(-.1,.8),from=c(0,19)),density(a_postcj_5)$x,col=cols[5],border=F)

#Figure 8B
par(mfrow=c(1,1),mar=c(6.1, 5.1, 4.1, 4.1))
x <- bound5[,c("baseline","diff2","diff3","diff4","diff5")]
names(x) <- c("Lowest","Low","Medium","High","Highest")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
ebar <- stripchart(x, ylim=c(0,.6), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.lab=2, cex.axis=1.5,ylab="Subsequent decision bound")
for(i in 1:10) abline(h=(-0.20+(i*.1)),col='lightgrey',lwd=0.8)  
mtext("Error positivity (Pe) amplitude",1,2.8,cex=2,at=3)
#stripchart(x,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg="grey",col='white',cex=3,lwd=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
lines(means,type='l',lty=2,lwd=3)
polygon(rescale(density(a_diff_2)$y,to=c(-.1,.8),from=c(0,15)),density(a_diff_2)$x,col=cols[2],border=F)
polygon(rescale(density(a_diff_3)$y,to=c(-.1,.8),from=c(0,15)),density(a_diff_3)$x,col=cols[3],border=F)
polygon(rescale(density(a_diff_4)$y,to=c(-.1,.8),from=c(0,15)),density(a_diff_4)$x,col=cols[4],border=F)
polygon(rescale(density(a_diff_5)$y,to=c(-.1,.8),from=c(0,15)),density(a_diff_5)$x,col=cols[5],border=F)


#DRIFT
#Previous pe; Figure 8 - Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- drift5[,c("baseline","prevpe2","prevpe3","prevpe4","prevpe5")]
names(x) <- c("Lowest","Low","Medium","High","Highest")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.3,.2), xlim=c(0.5,n+.5), vertical = TRUE, col="white",frame=F, cex.lab=2, cex.axis=1.5,ylab=substitute(paste("Drift rate on trial ",italic('n+1'))))
for(i in 1:10) abline(h=(-0.6+(i*.1)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Pe amplitude on trial ",italic('n'))),1,2.8,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
lines(means,type='l',lty=2,lwd=3)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
polygon(rescale(density(v_prevcj_2)$y,to=c(-.1,.8),from=c(0,7)),density(v_prevcj_2)$x,col=cols[2],border=F)
polygon(rescale(density(v_prevcj_3)$y,to=c(-.1,.8),from=c(0,7)),density(v_prevcj_3)$x,col=cols[3],border=F)
polygon(rescale(density(v_prevcj_4)$y,to=c(-.1,.8),from=c(0,7)),density(v_prevcj_4)$x,col=cols[4],border=F)
polygon(rescale(density(v_prevcj_5)$y,to=c(-.1,.8),from=c(0,7)),density(v_prevcj_5)$x,col=cols[5],border=F)

#Following pe; Figure 8 - Supplement 3
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
x <- drift5[,c("baseline","postpe2","postpe3","postpe4","postpe5")]
names(x) <- c("Lowest","Low","Medium","High","Highest")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.3,2.2), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.lab=2, cex.axis=1.5,ylab=substitute(paste("Drift rate on trial ",italic('n+1'))))
for(i in 1:10) abline(h=(-1+(i*.5)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Pe amplitude on trial ",italic('n+2'))),1,2.8,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
lines(means,type='l',lty=2,lwd=3)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
polygon(rescale(density(v_postcj_2)$y,to=c(-.1,.8),from=c(0,5)),density(v_postcj_2)$x,col=cols[2],border=F)
polygon(rescale(density(v_postcj_3)$y,to=c(-.1,.8),from=c(0,5)),density(v_postcj_3)$x,col=cols[3],border=F)
polygon(rescale(density(v_postcj_4)$y,to=c(-.1,.8),from=c(0,5)),density(v_postcj_4)$x,col=cols[4],border=F)
polygon(rescale(density(v_postcj_5)$y,to=c(-.1,.8),from=c(0,5)),density(v_postcj_5)$x,col=cols[5],border=F)

#3Figure 8B
par(mfrow=c(1,1),mar=c(6.1, 5.1, 4.1, 4.1))
x <- drift5[,c("baseline","diff2","diff3","diff4","diff5")]
names(x) <- c("Lowest","Low","Medium","High","Highest")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x, ylim=c(-2.5,.2), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.lab=2, cex.axis=1.5,ylab="Subsequent drift rate")
for(i in 1:10) abline(h=(-3.5+(i*.5)),col='lightgrey',lwd=0.8)  
mtext("Error positivity (Pe) amplitude",1,2.8,cex=2)
mtext("(error trials)",1,4.5,cex=2)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=3,lwd=2)
lines(means,type='l',lty=2,lwd=3)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=3.5,lwd=2)
polygon(rescale(density(v_diff_2)$y,to=c(-.1,.8),from=c(0,4)),density(v_diff_2)$x,col=cols[2],border=F)
polygon(rescale(density(v_diff_3)$y,to=c(-.1,.8),from=c(0,4)),density(v_diff_3)$x,col=cols[3],border=F)
polygon(rescale(density(v_diff_4)$y,to=c(-.1,.8),from=c(0,4)),density(v_diff_4)$x,col=cols[4],border=F)
polygon(rescale(density(v_diff_5)$y,to=c(-.1,.8),from=c(0,4)),density(v_diff_5)$x,col=cols[5],border=F)

#STATS
#pairwise comparing adjacent bins bound
length(which(a_prevcj_2<0))/length(a_prevcj_2) 
length(which(a_prevcj_3<a_prevcj_2))/length(a_prevcj_2) 
length(which(a_prevcj_4>a_prevcj_3))/length(a_prevcj_3) 
length(which(a_prevcj_5<a_prevcj_4))/length(a_prevcj_4) 

length(which(a_postcj_2>0))/length(a_postcj_2) #
length(which(a_postcj_3>a_postcj_2))/length(a_postcj_2) 
length(which(a_postcj_4>a_postcj_3))/length(a_postcj_3) 
length(which(a_postcj_5>a_postcj_4))/length(a_postcj_4) 

length(which(a_diff_2<0))/length(a_diff_2) #
length(which(a_diff_3<a_diff_2))/length(a_diff_2) #
length(which(a_diff_4<a_diff_3))/length(a_diff_3) #
length(which(a_diff_5<a_diff_4))/length(a_diff_4) #

#pairwise comparing adjacent bins drift
length(which(v_prevcj_2>0))/length(v_prevcj_2) #
length(which(v_prevcj_3>v_prevcj_2))/length(v_prevcj_2) #
length(which(v_prevcj_4<v_prevcj_3))/length(v_prevcj_3) #
length(which(v_prevcj_5>v_prevcj_4))/length(v_prevcj_4) #

length(which(v_postcj_2<0))/length(v_postcj_2) #
length(which(v_postcj_3<v_postcj_2))/length(v_postcj_2) 
length(which(v_postcj_4<v_postcj_3))/length(v_postcj_3) 
length(which(v_postcj_5<v_postcj_4))/length(v_postcj_4) 

length(which(v_diff_2>0))/length(v_diff_2) #
length(which(v_diff_3<v_diff_2))/length(v_diff_2) #
length(which(v_diff_4>v_diff_3))/length(v_diff_3) #
length(which(v_diff_5>v_diff_4))/length(v_diff_4) #




##############################################################################"
#7. ERPs As continous variables
#Load the data
#ERPs esimates
allData = read.csv("eegRegression_HDDMestimates.csv",header=T,sep=',')

bound3 <- data.frame(matrix(0,16,1));names(bound3) <- c("baseline")
drift3 <- data.frame(matrix(0,16,1));names(drift3) <- c("baseline")
bound3$prevpe <- allData$mean[grep("a_prevpe_bin_subj.",allData$X)]
bound3$prevern <- allData$mean[grep("a_prevernfcz_bin_subj.",allData$X)]
bound3$postpe <- allData$mean[grep("a_postpe_bin_subj.",allData$X)]
bound3$postern <- allData$mean[grep("a_posternfcz_bin_subj.",allData$X)]
drift3$prevpe <- allData$mean[grep("v_prevpe_bin_subj.",allData$X)]
drift3$prevern <- allData$mean[grep("v_prevernfcz_bin_subj.",allData$X)]
drift3$postpe <- allData$mean[grep("v_postpe_bin_subj.",allData$X)]
drift3$postern <- allData$mean[grep("v_posternfcz_bin_subj.",allData$X)]
bound3$pe_diff <- bound3$prevpe-bound3$postpe
bound3$ern_diff <- bound3$prevern-bound3$postern
drift3$pe_diff <- drift3$prevpe-drift3$postpe
drift3$ern_diff <- drift3$prevern-drift3$postern

a_prevpe <- read.csv("traces/a_prevpe_binnedpertrial.csv",header=F)$V1
a_prevern <- read.csv("traces/a_prevern_binnedpertrial.csv",header=F)$V1
v_prevpe <- read.csv("traces/v_prevpe_binnedpertrial.csv",header=F)$V1
v_prevern <- read.csv("traces/v_prevern_binnedpertrial.csv",header=F)$V1
a_postpe <- read.csv("traces/a_postpe_binnedpertrial.csv",header=F)$V1
a_postern <- read.csv("traces/a_postern_binnedpertrial.csv",header=F)$V1
v_postpe <- read.csv("traces/v_postpe_binnedpertrial.csv",header=F)$V1
v_postern <- read.csv("traces/v_postern_binnedpertrial.csv",header=F)$V1
a_pe_diff <- a_prevpe-a_postpe
a_ern_diff <- a_prevern-a_postern
v_pe_diff <- v_prevpe-v_postpe
v_ern_diff <- v_prevern-v_postern

#PLOTS
#Figure 8 - Supplement 1
par(mfrow=c(1,1),mgp=c(2,1,0))
x <- bound3[,c("prevern","prevpe")]
names(x) <- c("Error-related","Error positivity")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.0003,.0003),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F,ylab=substitute(paste("Decision bound on trial ",italic('n+1'))),main="")
mtext('Effect of EEG signature on',2,line=3)
mtext(c("(Pe)","negativity"),1,at=2:1,line=2)
mtext("(ERN)",1,at=1,line=3)
mtext(substitute(paste("Trial ",italic('n'))),1,at=1.5,line=4)
for(i in 1:10) abline(h=(-.0004+(i*.0001)),col='lightgrey',lwd=0.8)  
polygon(rescale(density(a_prevpe)$y,to=c(-.1,.8),from=c(0,19200)),density(a_prevpe)$x,col=rgb(0,1,0,.3),border=F)
polygon(rescale(density(a_prevern)$y,to=c(-.1,.8),from=c(0,19200)),density(a_prevern)$x,col=rgb(0,.5,0,.3),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)

#Figure 8 - Supplement 1
par(mfrow=c(1,1),mgp=c(2,1,0))
x <- bound3[,c("postern","postpe")]
names(x) <- c("Error-related","Error positivity")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.001,.0002),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F,ylab=substitute(paste("Decision bound on trial ",italic('n+1'))),main="")
mtext('Effect of EEG signature on',2,line=3)
mtext(c("(Pe)","negativity"),1,at=2:1,line=2)
mtext("(ERN)",1,at=1,line=3)
mtext(substitute(paste("Trial ",italic('n+2'))),1,at=1.5,line=4)
for(i in 1:20) abline(h=(-.002+(i*.0002)),col='lightgrey',lwd=0.8)  
polygon(rescale(density(a_postpe)$y,to=c(-.1,.8),from=c(0,15100)),density(a_postpe)$x,col=rgb(0,1,0,.3),border=F)
polygon(rescale(density(a_postern)$y,to=c(-.1,.8),from=c(0,15100)),density(a_postern)$x,col=rgb(0,.5,0,.3),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)

#Figure 8A
cols = c(rgb(0,.5,0,.3),rgb(0,1,0,.3))
par(mfrow=c(1,1),mgp=c(2,1,0))
x <- bound3[,c("ern_diff","pe_diff")]
names(x) <- c("Error-related","Error positivity")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.0005,.001),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F,ylab="subsequent decision bound",main="")
mtext('Effect of EEG signature on',2,line=3)
mtext(c("(Pe)","negativity"),1,at=2:1,line=2)
mtext("(ERN)",1,at=1,line=3)
for(i in 1:10) abline(h=(-.002+(i*.0005)),col='lightgrey',lwd=0.8)  
abline(h=0,lty=2)
polygon(rescale(density(a_pe_diff)$y,to=c(-.1,.8),from=c(0,11800)),density(a_pe_diff)$x,col=rgb(0,1,0,.3),border=F)
polygon(rescale(density(a_ern_diff)$y,to=c(-.1,.8),from=c(0,11800)),density(a_ern_diff)$x,col=rgb(0,.5,0,.3),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)

#DRIFT
#Figure 8 - Supplement 1
par(mfrow=c(1,1),mgp=c(2,1,0))
x <- drift3[,c("prevern","prevpe")]
names(x) <- c("Error-related","Error positivity")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.001,.00025),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F,ylab=substitute(paste("Drift rate on trial ",italic('n+1'))),main="")
mtext('Effect of EEG signature on',2,line=3)
mtext(c("negativity","(Pe)"),1,at=1:2,line=2)
mtext("(ERN)",1,at=1,line=3)
mtext(substitute(paste("Trial ",italic('n'))),1,at=1.5,line=4)
for(i in 1:20) abline(h=(-.002+(i*.00025)),col='lightgrey',lwd=0.8)  
polygon(rescale(density(v_prevpe)$y,to=c(-.1,.8),from=c(0,4340)),density(v_prevpe)$x,col=rgb(0,1,0,.3),border=F)
polygon(rescale(density(v_prevern)$y,to=c(-.1,.8),from=c(0,4340)),density(v_prevern)$x,col=rgb(0,.5,0,.3),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)

#Figure 8 - Supplement 1
par(mfrow=c(1,1),mgp=c(2,1,0))
x <- drift3[,c("postern","postpe")]
names(x) <- c("Error-related","Error positivity")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.001,.004),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F,ylab=substitute(paste("Drift rate on trial ",italic('n+1'))),main="")
mtext('Effect of EEG signature on',2,line=3)
mtext(c("negativity","(Pe)"),1,at=1:2,line=2)
mtext("(ERN)",1,at=1,line=3)
mtext(substitute(paste("Trial ",italic('n+2'))),1,at=1.5,line=4)
for(i in 1:20) abline(h=(-.005+(i*.001)),col='lightgrey',lwd=0.8)  
polygon(rescale(density(v_postpe)$y,to=c(-.1,.8),from=c(0,4000)),density(v_postpe)$x,col=rgb(0,1,0,.3),border=F)
polygon(rescale(density(v_postern)$y,to=c(-.1,.8),from=c(0,4000)),density(v_postern)$x,col=rgb(0,.5,0,.3),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)

#Figure 8A
par(mfrow=c(1,1),mgp=c(2,1,0))
x <- drift3[,c("ern_diff","pe_diff")]
names(x) <- c("Error-related","Error positivity")
means <- sapply(x, mean);n<- length(x);min(x);max(x)
stripchart(x,ylim=c(-.004,.001),xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F,ylab="Subsequent drift rate",main="")
mtext('Effect of EEG signature on',2,line=3)
mtext(c("negativity","(Pe)"),1,at=1:2,line=2)
mtext("(ERN)",1,at=1,line=3)
for(i in 1:20) abline(h=(-.005+(i*.001)),col='lightgrey',lwd=0.8)  
abline(h=0,lty=2)
polygon(rescale(density(v_pe_diff)$y,to=c(-.1,.8),from=c(0,2550)),density(v_pe_diff)$x,col=rgb(0,1,0,.3),border=F)
polygon(rescale(density(v_ern_diff)$y,to=c(-.1,.8),from=c(0,2550)),density(v_ern_diff)$x,col=rgb(0,.5,0,.3),border=F)
for(i in 1:n) stripchart(x[,i],at=i,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=cols[i],col='white',cex=2)
for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)

#Statistics
#bound
length(which(a_prevpe>0))/length(a_prevpe) #
length(which(a_prevern>0))/length(a_prevern) #
length(which(a_postpe>0))/length(a_postpe) #
length(which(a_postern>0))/length(a_postern) #
length(which(a_pe_diff<0))/length(a_pe_diff) #
length(which(a_ern_diff>0))/length(a_ern_diff) #

length(which(a_pe_diff<a_ern_diff))/length(a_pe_diff) #

#drift
length(which(v_prevpe<0))/length(v_prevpe) #
length(which(v_prevern>0))/length(v_prevern) #
length(which(v_postpe<0))/length(v_postpe) #
length(which(v_postern>0))/length(v_postern) #
length(which(v_pe_diff>0))/length(v_pe_diff) #
length(which(v_ern_diff>0))/length(v_ern_diff) #

length(which(v_pe_diff>v_ern_diff))/length(v_pe_diff) #
