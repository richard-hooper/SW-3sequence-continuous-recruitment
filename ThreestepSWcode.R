######################################################
## R code to accompany "Efficient designs for three-
## sequence stepped wedge trials with continuous recruitment"
## by R Hooper, O Quintin, J Kasza
##
## This R file includes functions for:
## 1. simulating data from the a three-sequence
## stepped wedge trial with given properties;
## 2. analysing that data;
## 3. replicating the simulation and analysis;
## 4. the calculation of theoretical power to detect an effect
######################################################


require(glmmTMB)
require(MASS)

#A wrapper function to replicate data generation and analysis nsims times
ThreeStepSW_wrap <- function(nrep,rho,tau,totvar,mfirst, mmid,Kfirst, Kmid,theta){
  #Wrapper function forthe ThreeStepSWdata and ThreeStepSWanalysis functions. 
  #nrep = number of replications in the simulation study
  # ---- design parameters 
  #rho = intracluster correlation (ICC)
  #tau = cluster autocorrelation (CAC)
  #totvar = total variance
  #mfirst = number of participants per cluster in first (and last) period
  #mmid = number of participants per cluster in each middle period
  #Kfirst = number of clusters assigned to first (and last) sequence
  #Kmid = number of clusters assigned to middle sequence
  #theta = treatment effect of interest
  
  output <- replicate(nrep, ThreeStepSWanalysis(rho,tau,totvar,mfirst, mmid,Kfirst, Kmid,theta))
  
  
  #Here we are primary interested in validating large-sample power calculations
  #- bias in treatment effect estimation
  #- empirical standard error 
  #- MSE
  #- Average Model SE
  #- Rejection percentage
  
  Rxeff <- theta
  
  myresults <- NULL
  #Number of converged
  nconv <- nrep - sum(is.na(output[1,]))
  
  #Bias
  bias <- mean(output[1,], na.rm=TRUE)-Rxeff
  #Empirical SE
  empSE <- sqrt(sum((output[1,]-mean(output[1,], na.rm=TRUE))^2, na.rm=TRUE)/(nconv-1))
  #MSE
  mse <- (sum(output[1,]-Rxeff, na.rm=TRUE)^2)/(nconv)
  #Average model SE
  avemodelSE <- sqrt(sum((output[2,])^2, na.rm=TRUE)/(nconv))
  #Rejection percentage
  rejPercent <-  sum(abs(output[1,])/output[2,]>1.96, na.rm=TRUE)/nconv
  #Coverage
  temp1 <- (output[1,] -1.96*output[2,] < Rxeff)
  temp2 <- (output[1,] +1.96*output[2,] > Rxeff)
  coverage <- sum(temp1 & temp2, na.rm=TRUE)/nconv
  rm(temp1)
  rm(temp2)
  
  #Number of converged models
 # myresults[7] <- nrep-sum(is.na(output[1,]))

  return(c(bias, empSE, mse, avemodelSE, rejPercent, coverage, nconv))
  
}


#Simulate one dataset from a three-step SW design with continuous recruitment:
ThreeStepSWdata <- function(rho,tau,totvar,mfirst, mmid,Kfirst, Kmid,theta){
  
  # ---- design parameters 
  #rho = intracluster correlation (ICC)
  #tau = cluster autocorrelation (CAC)
  #totvar = total variance
  #mfirst = number of participants per cluster in first (and last) period
  #mmid = number of participants per cluster in each middle period
  #Kfirst = number of clusters assigned to first (and last) sequence
  #Kmid = number of clusters assigned to middle sequence
  #theta = treatment effect of interest
  
  #Total number of clusters and participants per cluster
  Ktot <- 2*Kfirst + Kmid
  mtot <- 2*mfirst + 2*mmid
  
  
  #Generate the vector of treatment indicator variables for each sequence
  if(mfirst > 0){
    Xseq1 <- c(rep(0, times=mfirst), rep(1, times= (mtot-mfirst)))
    Xseq2 <- c(rep(0, times=(mfirst + mmid)), rep(1, times=(mfirst + mmid)))
    Xseq3 <- c(rep(0, times=(mtot-mfirst)), rep(1, times=mfirst))
    
    #Generate the matrix of time effects for each cluster
    Tmat <- matrix(data = 0, nrow = mtot, ncol = 4)
    Tmat[1:mfirst,1] <- 1
    Tmat[(mfirst + 1):(mmid + mfirst),2] <- 1
    Tmat[(mmid + mfirst + 1):(2*mmid + mfirst),3] <- 1
    Tmat[(2*mmid + mfirst + 1):mtot,4] <- 1
  }
  else {
    Xseq1 <- rep(0, mtot)
    Xseq2 <- c(rep(0, floor(mtot/2)), rep(1, mtot-floor(mtot/2)))
    Xseq3 <- rep(1, mtot)
    
    #Generate the matrix of time effects for each cluster
    Tmat <- matrix(data = 0, nrow = mtot, ncol = 4)
    Tmat[1:floor(floor(mtot/2)/2),1] <- 1
    Tmat[(floor(floor(mtot/2)/2)+1):floor(mtot/2),2] <- 1
    Tmat[(floor(mtot/2) + 1):floor(3*mtot/4),3] <- 1
    Tmat[(floor(3*mtot/4) + 1):mtot,4] <- 1
  }
  
  
  #Putting together the design matrix for all observations
  Xvecall <- c(rep(Xseq1, times=Kfirst), rep(Xseq2, times=Kmid), rep(Xseq3, times=Kfirst))
  Xmattotal <- cbind(matrix(rep(t(Tmat), times=Ktot), ncol=4, byrow=TRUE), Xvecall)

  #Cluster indicators
  Cind <- rep(seq(1:Ktot), each = mtot)
  #Period indicators
  Pind <- rep(c(rep(1, times=mfirst), rep(2, times=mmid), rep(3, times=mmid), rep(4, times=mfirst)), times=Ktot)
  #Time of each observation (equally spaced observations in each cluster)
  Tind <- rep(seq(0,1, length.out = mtot), times=Ktot)

  treatf <- as.factor(Xvecall)
  periodf <- as.factor(Pind)
  clusterf <- as.factor(Cind)
  arrivalf <- as.factor(Tind)
  
  des_factors <- data.frame(treatf, periodf, clusterf, arrivalf)
    
  #Generate the correlation matrix of a single cluster
  corr1 <- diag(1-rho,mtot)  +  rho*(tau^abs(matrix(seq(0,1, length.out = mtot),nrow=mtot, ncol=mtot, byrow=FALSE) - matrix(seq(0,1, length.out = mtot),nrow=mtot, ncol=mtot, byrow=TRUE)))
  #The covariance matrix of a single cluster is then:
  cov1 <- totvar*corr1
  
  #Generate the random part of each observation
  #Each row corresponds to observations from one cluster
  randpart <- mvrnorm(n=Ktot, mu = rep(0, mtot), Sigma = cov1)
  randpartvec <- as.vector(t(randpart))
  
  # ---- add all the fixed and random effects together to get the outcome for 
  #each observation
  #Note that here I'm just adding in the Pinds for the period effects
  YY = Xvecall*theta+ Pind + randpartvec # vector of observations
  
  # make a data frame
  mydat = cbind(data.frame(YY),des_factors) 
  return(mydat)
}

#Generate and analyse a single dataset 
ThreeStepSWanalysis <- function(rho,tau,totvar, mfirst, mmid, Kfirst, Kmid, theta){
  # ---- design parameters 
  #rho = intracluster correlation (ICC)
  #tau = cluster autocorrelation (CAC)
  #totvar = total variance
  #mfirst = number of participants per cluster in first (and last) period
  #mmid = number of participants per cluster in each middle period
  #Kfirst = number of clusters assigned to first (and last) sequence
  #Kmid = number of clusters assigned to middle sequence
  #theta = treatment effect of interest
  
  mydata <- ThreeStepSWdata(rho,tau,totvar, mfirst, mmid, Kfirst, Kmid, theta)
  
  mymodel <- glmmTMB(YY~treatf+periodf+ar1(arrivalf+0|clusterf), data=mydata, REML=T) 
  trt_est    <- summary(mymodel)$coefficients$cond["treatf1",1]
  se_trt_est <- summary(mymodel)$coefficients$cond["treatf1",2]
  
  return( c(trt_est, se_trt_est) )
  
}


#Calcualte the theoretical power of a three-step SW design with continuous recruitment:
ThreeStepSWpower <- function(rho,tau,totvar,mfirst, mmid,Kfirst, Kmid,theta, alpha){
  
  # ---- design parameters 
  #rho = intracluster correlation (ICC)
  #tau = cluster autocorrelation (CAC)
  #totvar = total variance
  #mfirst = number of participants per cluster in first (and last) period
  #mmid = number of participants per cluster in each middle period
  #Kfirst = number of clusters assigned to first (and last) sequence
  #Kmid = number of clusters assigned to middle sequence
  #theta = treatment effect of interest
  #alpha = 2-sided sigificance level
  
  #Total number of clusters and participants per cluster
  Ktot <- 2*Kfirst + Kmid
  mtot <- 2*mfirst + 2*mmid
  
  
  #Generate the vector of treatment indicator variables for each sequence
  if(mfirst > 0){
    Xseq1 <- c(rep(0, times=mfirst), rep(1, times= (mtot-mfirst)))
    Xseq2 <- c(rep(0, times=(mfirst + mmid)), rep(1, times=(mfirst + mmid)))
    Xseq3 <- c(rep(0, times=(mtot-mfirst)), rep(1, times=mfirst))
    
    #Generate the matrix of time effects for each cluster
    Tmat <- matrix(data = 0, nrow = mtot, ncol = 4)
    Tmat[1:mfirst,1] <- 1
    Tmat[(mfirst + 1):(mmid + mfirst),2] <- 1
    Tmat[(mmid + mfirst + 1):(2*mmid + mfirst),3] <- 1
    Tmat[(2*mmid + mfirst + 1):mtot,4] <- 1
  }
  else {
    Xseq1 <- rep(0, mtot)
    Xseq2 <- c(rep(0, floor(mtot/2)), rep(1, mtot-floor(mtot/2)))
    Xseq3 <- rep(1, mtot)
    
    #Generate the matrix of time effects for each cluster
    Tmat <- matrix(data = 0, nrow = mtot, ncol = 4)
    Tmat[1:floor(floor(mtot/2)/2),1] <- 1
    Tmat[(floor(floor(mtot/2)/2)+1):floor(mtot/2),2] <- 1
    Tmat[(floor(mtot/2) + 1):floor(3*mtot/4),3] <- 1
    Tmat[(floor(3*mtot/4) + 1):mtot,4] <- 1
  }
  
  
  #Putting together the design matrix for all observations
  Xvecall <- c(rep(Xseq1, times=Kfirst), rep(Xseq2, times=Kmid), rep(Xseq3, times=Kfirst))
  Xmattotal <- cbind(matrix(rep(t(Tmat), times=Ktot), ncol=4, byrow=TRUE), Xvecall)
  
  
  #Generate the correlation matrix of a single cluster
#  corr1 <- diag(1-rho,mtot)  + rho*(tau^abs(matrix(1:mtot,nrow=mtot, ncol=mtot, byrow=FALSE) - matrix(1:mtot,nrow=mtot, ncol=mtot, byrow=TRUE)))
  corr1 <- diag(1-rho,mtot)  +  rho*(tau^abs(matrix(seq(0,1, length.out = mtot),nrow=mtot, ncol=mtot, byrow=FALSE) - matrix(seq(0,1, length.out = mtot),nrow=mtot, ncol=mtot, byrow=TRUE)))
  V1 <-  totvar*corr1
  VarMat <- kronecker(diag(1,Ktot), V1)
  
  sigma2 <- (solve((t(Xmattotal)%*%solve(VarMat)%*%Xmattotal))[ncol(Xmattotal),ncol(Xmattotal)])
  
  power <- pnorm(abs(theta)/sqrt(sigma2)-qnorm(1-alpha/2))
  
  return(power)
}