##################################################
###############qSIP宏基因组代码###################
##################################################

data <- read.table("./data.txt", header=TRUE,sep = "\t")
data$Yijk <- data$Rel_Abun * data$DNA_content
##WAD.func
WAD.func <- function(y, x){
  WAD <- sum(x*(y/sum(y)))
  WAD
}
#sample.vec ##Create a version of 'sample' that avoids the problem when length(x)==1
sample.vec <- function(x, ...){
  x[sample(length(x), ...)]
}

#boot.TUBE.func
boot.TUBE.func <- function(X,  vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000){ 
  # Create a dataframe of only x, y, and rep: 
  test.data <- data.frame(x=X[,vars[1]], y=X[,vars[2]], rep=factor(X[,vars[3]]))
  
  # Calculate observed weighted average density (WAD), total 16S copies, and mass of soil for each rep:
  obs.wads <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=2))
  names(obs.wads) <- c("wad",  "rep")
  for (r in 1:length(levels(test.data$rep))){
    obs.wads$rep[r] <- levels(test.data$rep)[r]
    obs.wads$wad[r] <- WAD.func(y=test.data$y[test.data$rep == levels(test.data$rep)[r]], x=test.data$x[test.data$rep == levels(test.data$rep)[r]])
  }
  obs.wads$rep <- factor(obs.wads$rep)
  
  # Bootstrapping: Calculate a bootstrap vector of mean WADs across reps along with total 16S copies & mass of soil:
  boot.wads <- numeric()
  for (i in 1:draws){
    indices <- sample.vec(1:dim(obs.wads)[1], dim(obs.wads)[1], replace=TRUE)
    boot.wads[i] <- mean(obs.wads$wad[indices], na.rm=T)
  }
  boot.wads.CI <- quantile(boot.wads, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
  
  reps.NAs <- obs.wads$rep[is.na(obs.wads$wad)]
  if (length(reps.NAs) == 0){
    message <- "none"
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")
  
  return(list(boot.wads=boot.wads, 
              obs.wads=obs.wads, 
              obs.wad.mean=mean(obs.wads$wad, na.rm=T), 
              boot.wads.mean=mean(boot.wads, na.rm=T), 
              boot.wads.median=median(boot.wads, na.rm=T), 
              boot.wads.CI=boot.wads.CI, 
              message=message))
}

library(tidyr)
data=tidyr::unite(data, "MAG_Rep",MAG_ID,Rep,remove = FALSE)
output = boot.TUBE.func(data, vars=c("Density", "Yijk", "MAG_Rep"),  CI=0.90, draws=1000)

WLABi=output$obs.wads
WLABi=separate(WLABi, rep,c("MAG_ID","Rep"), sep = "_")
colnames(WLABi)=c("WLABi","MAG_ID","Rep")

data=merge(data,WLABi,by=c("MAG_ID","Rep"),all=T)

data=data[,c("MAG_ID","Rep","Sample_ID","GC_content" , "WLABi" )]
data=data[!duplicated(data), ]

data$WLighti=data$GC_content*0.088+1.6689
data$Zi=data$WLABi-data$WLighti
data$MLIGHTi=0.496*data$GC_content+307.691
data$MHEAVYMAXi=-0.4987282*data$GC_content+9.974564+data$MLIGHTi
data$MLABi=((data$Zi/data$WLighti)+1)*data$MLIGHTi
data$EAF=((data$MLABi-data$MLIGHTi)/(data$MHEAVYMAXi-data$MLIGHTi))*(1-0.01111233)

func <- function(x)(c(n = length(x),mean=mean(x,na.rm = T),sd=sd(x,na.rm = T),cv=sd(x,na.rm = T)/mean(x,na.rm = T)*100,se=sd(x,na.rm = T)/sqrt(length(x))))
mean=aggregate(EAF ~ MAG_ID, data=data,func)##
result <- as.data.frame(do.call(cbind, mean))
write.csv(result,'EAF.csv')  
