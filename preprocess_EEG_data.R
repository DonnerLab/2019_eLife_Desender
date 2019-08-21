setwd("...")
rm(list=ls())

#Due to privacy issues, only averaged EEG data are available on GH. 
#Please contact Annika Boldt directly at a.boldt@ucl.ac.uk to get single trial data
#Then, filter the data at 10hz and read these data in below:
Data <- read.table('....txt',header=T) 
Data$ACC <-  1;
Data$ACC[Data$cor==1] <-  0;
Data$cj <- Data$conf

Data$stimulus <- NA #code as 1 and -1
Data$stimulus[Data$resp==1&Data$ACC==1] <- 1
Data$stimulus[Data$resp==0&Data$ACC==0] <- 1
Data$stimulus[Data$resp==0&Data$ACC==1] <- -1
Data$stimulus[Data$resp==1&Data$ACC==0] <- -1
table(Data$stimulus,Data$resp,Data$ACC)
#
Data$prevcor[2:length(Data$ACC)] <- Data$ACC[1:length(Data$ACC)-1]
Data$prevprevcor[3:length(Data$ACC)] <- Data$ACC[1:(length(Data$ACC)-2)]
Data$postcor <- NA
Data$postcor[1:length(Data$ACC)-1] <- Data$ACC[2:length(Data$ACC)]
Data$prevcj[2:length(Data$cj)] <- Data$cj[1:length(Data$cj)-1]
Data$postcj <- NA
Data$postcj[1:length(Data$cj)-1] <- Data$cj[2:length(Data$cj)]
Data$postsub <- NA
Data$postsub[1:length(Data$sub)-1] <- Data$sub[2:length(Data$sub)]
Data$postsub[dim(Data)[1]] <- 99

#Regress ERN from Pe and Pe from ERN
Data$ERNfcz_nope <- NA;Data$PEpz_noern <- NA;pe_ern_cor <- NA; pe_ern_cor_filt <- NA;
for(s in 1:16){
  temp <- subset(Data,sub==s)
  pe_ern_cor[s] <- cor(temp$PEpz,temp$ERNfcz) #quantify pe/ern on trial level
  
  tempPe <- lm(PEpz~ERNfcz,data=temp)
  tempErn <- lm(ERNfcz~PEpz,data=temp)
  Data$PEpz_noern[Data$sub==s] <- tempPe$residuals
  Data$ERNfcz_nope[Data$sub==s] <-tempErn$residuals 
  
  temp <- subset(Data,sub==s)
  pe_ern_cor_filt[s] <- cor(temp$PEpz_noern,temp$ERNfcz_nope) #quantify pe/ern on trial level
}  
mean(pe_ern_cor);sd(pe_ern_cor)  
mean(pe_ern_cor_filt);sd(pe_ern_cor_filt)  

#Replicate Boldt & Yeung, 2015
par(mfrow=c(2,2))
Pe <- with(Data,aggregate(PEpz_noern,by=list(cj,sub),mean))
Pe <- cast(Pe, Group.2~Group.1)
plot(colMeans(Pe))
ERN <- with(Data,aggregate(ERNfcz_nope,by=list(cj,sub),mean))
ERN <- cast(ERN, Group.2~Group.1)
plot(colMeans(ERN))
ERN <- with(Data,aggregate(ERNfcz_nope,by=list(cor,sub),mean))
ERN <- cast(ERN, Group.2~Group.1)
t.test(ERN[,2],ERN[,3],paired=T)

#### Recode ERPs into bins & standardize neural data!
Data$Pe_bin <- NA; Data$ERNfcz_bin <- NA; Data$ERNfcz_bin_nope <- NA;
Data$Pe_z <- NA; Data$ERNfcz_z <- NA; Data$PE_bin_noern
for(s in 1:16){
  Data$trialnr[Data$sub==s] <- 1:length(Data$sub[Data$sub==s])
}
for(s in 1:16){
  temp <- subset(Data,sub==s)
  quantPe <- quantile(temp$PEpz,probs=seq(0,1,.2)); #.2 before, 5 bins
  quantERNfcz <- quantile(temp$ERNfcz,probs=seq(0,1,.2));
  quantPe_noern <- quantile(temp$PEpz_noern,probs=seq(0,1,.2)); #.2 before, 5 bins
  quantERNfcz_nope <- quantile(temp$ERNfcz_nope,probs=seq(0,1,.2));
  for(i in 1:5){
    Data$Pe_bin[Data$sub==s&Data$PEpz >= quantPe[i]] <- i
    Data$ERNfcz_bin[Data$sub==s&Data$ERNfcz >= quantERNfcz[i]] <- i
    Data$Pe_bin_noern[Data$sub==s&Data$PEpz_noern >= quantPe_noern[i]] <- i
    Data$ERNfcz_bin_nope[Data$sub==s&Data$ERNfcz_nope >= quantERNfcz_nope[i]] <- i
  }
  Data$Pe_z[Data$sub==s] <- scale(Data$PEpz[Data$sub==s])
  Data$ERNfcz_z[Data$sub==s] <- scale(Data$ERNfcz[Data$sub==s])
}

#Order all trials according to specific variable, and then number these
Data <- Data[order(Data$PEpz),]
for(s in 1:16){
  Data$Pe_bin[Data$sub==s] <- 1:length(Data$sub[Data$sub==s])
}
Data <- Data[order(Data$ERNfcz),]
for(s in 1:16){
  Data$ERNfcz_bin[Data$sub==s] <- 1:length(Data$sub[Data$sub==s])
}
Data <- Data[order(Data$sub,Data$trialnr),]

#
Data$prevconf[2:length(Data$conf)] <- Data$conf[1:length(Data$conf)-1]
Data$prevpe[2:length(Data$PEpz)] <- Data$PEpz[1:length(Data$PEpz)-1]
Data$postpe <- NA
Data$postpe[1:length(Data$PEpz)-1] <- Data$PEpz[2:length(Data$PEpz)]
Data$prevernfcz[2:length(Data$ERNfcz)] <- Data$ERNfcz[1:length(Data$ERNfcz)-1]
Data$posternfcz <- NA
Data$posternfcz[1:length(Data$ERNfcz)-1] <- Data$ERNfcz[2:length(Data$ERNfcz)]
Data$prevernpz[2:length(Data$ERNpz)] <- Data$ERNpz[1:length(Data$ERNpz)-1]
Data$posternpz <- NA
Data$posternpz[1:length(Data$ERNpz)-1] <- Data$ERNpz[2:length(Data$ERNpz)]
#same for z scored data
Data$prevpe_z[2:length(Data$Pe_z)] <- Data$Pe_z[1:length(Data$Pe_z)-1]
Data$postpe_z <- NA
Data$postpe_z[1:length(Data$Pe_z)-1] <- Data$Pe_z[2:length(Data$Pe_z)]
Data$prevernfcz_z[2:length(Data$ERNfcz_z)] <- Data$ERNfcz_z[1:length(Data$ERNfcz_z)-1]
Data$posternfcz_z <- NA
Data$posternfcz_z[1:length(Data$ERNfcz_z)-1] <- Data$ERNfcz_z[2:length(Data$ERNfcz_z)]
#same for binned data
Data$prevpe_bin[2:length(Data$Pe_bin)] <- Data$Pe_bin[1:length(Data$Pe_bin)-1]
Data$postpe_bin <- NA
Data$postpe_bin[1:length(Data$PE_bin)-1] <- Data$Pe_bin[2:length(Data$Pe_bin)]
Data$prevernfcz_bin[2:length(Data$ERNfcz_bin)] <- Data$ERNfcz_bin[1:length(Data$ERNfcz_bin)-1]
Data$posternfcz_bin <- NA
Data$posternfcz_bin[1:length(Data$ERNfcz_bin)-1] <- Data$ERNfcz_bin[2:length(Data$ERNfcz_bin)]
#same for binned data after regressing the other out
Data$prevpe_bin_noern[2:length(Data$Pe_bin_noern)] <- Data$Pe_bin_noern[1:length(Data$Pe_bin_noern)-1]
Data$postpe_bin_noern <- NA
Data$postpe_bin_noern[1:length(Data$PE_bin_noern)-1] <- Data$Pe_bin_noern[2:length(Data$Pe_bin_noern)]
Data$prevernfcz_bin_nope[2:length(Data$ERNfcz_bin_nope)] <- Data$ERNfcz_bin_nope[1:length(Data$ERNfcz_bin_nope)-1]
Data$posternfcz_bin_nope <- NA
Data$posternfcz_bin_nope[1:length(Data$ERNfcz_bin_nope)-1] <- Data$ERNfcz_bin_nope[2:length(Data$ERNfcz_bin_nope)]
#check overlap
table(Data$prevpe_bin,Data$prevpe_bin_noern)
#
Data$prevpe_filt[2:length(Data$PEpz_filt)] <- Data$PEpz_filt[1:length(Data$PEpz_filt)-1]
Data$postpe_filt <- NA
Data$postpe_filt[1:length(Data$PEpz_filt)-1] <- Data$PEpz_filt[2:length(Data$PEpz_filt)]
Data$prevernfcz_filt[2:length(Data$ERNfcz_filt)] <- Data$ERNfcz_filt[1:length(Data$ERNfcz_filt)-1]
Data$posternfcz_filt <- NA
Data$posternfcz_filt[1:length(Data$ERNfcz_filt)-1] <- Data$ERNfcz_filt[2:length(Data$ERNfcz_filt)]
Data$prevernpz_filt[2:length(Data$ERNpz_filt)] <- Data$ERNpz_filt[1:length(Data$ERNpz_filt)-1]
Data$posternpz_filt <- NA
Data$posternpz_filt[1:length(Data$ERNpz_filt)-1] <- Data$ERNpz_filt[2:length(Data$ERNpz_filt)]

#outliers?
plot(density(Data$prevpe_z,na.rm=T),xlim=c(-5,5))
table(abs(Data$prevpe_z) > 3)
plot(density(Data$prevernfcz_z,na.rm=T),xlim=c(-5,5))
table(abs(Data$prevernfcz_z) > 3)


#Exclude all first trial from the experiment
Data$outlier <- 0
Data$outlier[Data$trial == 1] <- 1 #first trial
Data$outlier[Data$sub != Data$postsub] <- 2 #first trial
Data$outlier[abs(Data$prevpe_z) > 3] <- 3 #outliers prev pe
Data$outlier[abs(Data$prevernfcz_z) > 3] <- 4 #outliers prev ern
table(Data$outlier)
Data <- subset(Data, outlier==0)

#Replicate Boldt & Yeung, 2015 again
par(mfrow=c(2,2))
library(reshape)
Pe <- with(Data,aggregate(PEpz,by=list(cj,sub),mean))
Pe <- cast(Pe, Group.2~Group.1)
plot(colMeans(Pe))
ERN <- with(Data,aggregate(ERNfcz,by=list(cj,sub),mean))
ERN <- cast(ERN, Group.2~Group.1)
plot(colMeans(ERN))
ERN <- with(Data,aggregate(ERNfcz,by=list(cor,sub),mean))
ERN <- cast(ERN, Group.2~Group.1)
t.test(ERN[,2],ERN[,3],paired=T)

#Save data to fit to hddm (code response (0 1))
Data$subj_idx <- Data$sub
Data$response <- Data$resp
Data$rt <- Data$rt/1000
write.table(Data,"EEGdata_hddm.csv",quote=F,row.names=F,sep=',')
