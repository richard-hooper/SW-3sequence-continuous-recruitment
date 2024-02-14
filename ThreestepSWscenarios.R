######################################################
## R code to accompany "Efficient designs for three-
## sequence stepped wedge trials with continuous recruitment"
## by R Hooper, O Quintin, J Kasza
##
## This R file includes code to simulate data from a set
## of configurations given in the Scenarios.csv file.
######################################################

#The commented code is for the calculation of theoretical power
# scenarios <- read.csv("Scenarios.csv", header = TRUE)
# scenarios$s[scenarios$Standard == "y"] <- 0.25
# scenarios$s[scenarios$Standard != "y"] <- 1/12
# scenarios$w = 1/3
# 
# 
# scenarios$mfirst <- floor(scenarios$mtot*scenarios$s)
# scenarios$mmid <- floor(scenarios$mtot*(0.5 - scenarios$s))
# scenarios$Kfirst <- floor(scenarios$Ktot*scenarios$w)
# scenarios$Kmid <- floor(scenarios$Ktot*(1-2*scenarios$w))
# 
# scenarios$totvar <- 10.7^2

#Add theoretical power to the scenarios dataset
#Note: this is quite computationally intensive due to the need to 
#invert large matrices. I've calculated and saved these but these
#calculations can be checked
# for(i in 1:nrow(scenarios)){
#  scenarios$power[i] <- ThreeStepSWpower(scenarios$rho[i],scenarios$tau[i], scenarios$totvar[i],
#                                          scenarios$mfirst[i], scenarios$mmid[i],scenarios$Kfirst[i],
#                                          scenarios$Kmid[i], scenarios$theta[i], 0.05)
# 
#  }

# write.csv(scenarios, file="scenarios_power.csv")

##Code for simulation to obtain empirical power
scenarios <- read.csv("scenarios_power.csv", header = TRUE)

nrep <- 1000

set.seed(1234)
myresults <- matrix(data=NA, nrow=nrow(scenarios), ncol=7)


for(i in 1:nrow(scenarios)){
    myresults[i,] <- ThreeStepSW_wrap(nrep,scenarios$rho[i],scenarios$tau[i],
                  scenarios$totvar[i],scenarios$mfirst[i], scenarios$mmid[i], 
                  scenarios$Kfirst[i], scenarios$Kmid[i], scenarios$theta[i])
  
}

