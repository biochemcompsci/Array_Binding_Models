library(deSolve)
library(reshape2)
library(ggplot2)

source("Generate_Mixed_Distribution.R")

# Ab + Pep <-> AbPep* <-> AbPep
# dAb/dt = k1off(AbPep*) - k1on(Ab X Pep)
# dPep/dt = k1off(AbPep*) - k1on(Ab X Pep)
# dAbPep*/dt = k1on(Ab X Pep) - k1off(AbPep*) - k2on(AbPep*) + k2off(AbPep)
# dAbPep/dt = k2on(AbPep*) - k2off(AbPep)
# k1on = kdiffusion = constant
# k1off = kdesorption = constant
# k2on = interaction on-rate
# k2off = interaction off-rate

# k1offParam <- 0.5
# k1onParam <- 300
# k2offParam <- 0.000001
# k2onParam <- "Mixed Modal"
# maxAbPepConc <- 3.8E-19

# k1offParam <- 5
# k1onParam <- 500
# k2offParam <- 0.000001
# k2onParam <- "Mixed Modal"
# maxAbPepConc <- 2.2E-19
# minAbPepConc <- 0

k1offParam <- 5
k1onParam <- 500
k2offParam <- 0.000001
k2offParam <- "Mixed Modal"
k2onParam <- "Mixed Modal"
maxAbPepConc <- 6.0E-20
minAbPepConc <- 0

#concInit <- c("Ab" = 100E-12, "Pep" = 100E-15, "AbPep*" = 0, "AbPep" = 0)
concInit <- c("Ab" = 100E-12, "Pep" = 100E-15, "AbPep*" = 0, "AbPep" = 0)
timeSeries <- seq(0, 120, 0.5)
#params <- list(k1off = 1, k1on = 300, k2off = 0.000001, k2on = 2.4)
#params <- list(k1off = 1, k1on = 300, k2off = 0.000001, k2on = 2.2)

# Compare Simulated to Experimental
Chagas_OneTwenty_Min_10X_Unscaled.Combined_Specificity <- read.csv("~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Green_Channel_Chagas/120_min_incubation/10X_Wash/sample_level/Chagas_OneTwenty_Min_10X_Unscaled-Combined_Specificity.csv", stringsAsFactors=FALSE)
Chagas6041Specific <- Chagas_OneTwenty_Min_10X_Unscaled.Combined_Specificity[Chagas_OneTwenty_Min_10X_Unscaled.Combined_Specificity$Sample.ID == 'Chagas.G1.6041',1:4]
Chagas6051Specific <- Chagas_OneTwenty_Min_10X_Unscaled.Combined_Specificity[Chagas_OneTwenty_Min_10X_Unscaled.Combined_Specificity$Sample.ID == 'Chagas.G2.6051',1:4]

ExperimentalData <- Chagas6051Specific
ExperimentalDataScaled <- ExperimentalData
ExperimentalDataScaled$Raw.Mean <- ExperimentalDataScaled$Raw.Mean / max(ExperimentalDataScaled$Raw.Mean)

# Generate a Bimodal Distribution of On-Rates Based on Experimental Observation
#sampleSize <- 2000
sampleSize <- nrow(ExperimentalData)

dist1Weight <- 0.4
#dist2Weight <- 0.
dist1Mean <- 0.15
dist2Mean <- 25
dist1SD <- 1
dist2SD <- 0.1
min <- 0
max <- 35

#konSample <- generateBimodalDistribution(sampleSize, dist1Mean, dist2Mean, dist1SD, dist2SD, dist1Weight, dist2Weight, min, max)
konSample <- generateBimodalDistribution("Norm", sampleSize, dist1Mean, dist2Mean, dist1SD, dist2SD, dist1Weight, min, max)
plot(density(konSample))

# Generate a Bimodal Distribution of Off-Rates Based on Experimental Observation
dist1Weight <- 0.9999
#dist2Weight <- 0.35
dist1Mean <- 0.000005
dist2Mean <- 0.0005
dist1SD <- 2.5
dist2SD <- 3
min <- 0.000003
max <- 0.5

koffSample <- generateBimodalDistribution("LogNorm", sampleSize, dist1Mean, dist2Mean, dist1SD, dist2SD, dist1Weight, min, max)
#koffSample <- rlnorm(sampleSize, meanlog = dist1Mean, sdlog = dist1SD)
plot(density(koffSample))

# Function that Encapsulates the System of Differential Equations
kineticModel <- function(time, conc, params) {
  
  k1on <- params$k1on
  k1off <- params$k1off
  k2on <- params$k2on
  k2off <- params$k2off
  
  odeSystem <- rep(0, length(conc))
  
  odeSystem[1] <- k1off*conc["AbPep*"] - k1on*(conc["Ab"]*conc["Pep"]) # dAb/dt
  odeSystem[2] <- k1off*conc["AbPep*"] - k1on*(conc["Ab"]*conc["Pep"]) # dPep/dt
  odeSystem[3] <- k1on*(conc["Ab"]*conc["Pep"]) - k1off*conc["AbPep*"] - k2on*conc["AbPep*"] + k2off*conc["AbPep"] # dAbPep*/dt
  odeSystem[4] <- k2on*conc["AbPep*"] - k2off*conc["AbPep"] # dAbPep/dt
  
  return(list(odeSystem))
  
}

#k2onSample <- rlnorm(100, meanlog = 1.5, sdlog = 1)
k2onSample <- konSample
k2offSample <- koffSample

simulationOutput <- data.frame()
idxList <- c()
currIdx <- 1

for(currk2on in k2onSample) {
  
  #Slowest On-Rate has Fastest Off-Rate (choose opposite sides of kon and koff distribution)
  currk2off <- k2offSample[length(k2offSample) - currIdx + 1]
  paste0('kon = ', currk2on, 'koff = ', currk2off)
  
  params <- list(k1off = k1offParam, k1on = k1onParam, k2offParam = currk2off, k2on = currk2on)
  
  if(length(simulationOutput) > 0) {
    
    simulationOutput <- rbind(simulationOutput, ode(y = concInit, times = timeSeries, func = kineticModel, parms = params, method = "lsoda"))
    
  } else {
    
    simulationOutput <- ode(y = concInit, times = timeSeries, func = kineticModel, parms = params, method = "lsoda")
    
  }
  
  idxList <- append(idxList, currIdx)
  currIdx <- currIdx + 1
  
}

filteredOutput <-simulationOutput[,'AbPep']
filteredOutput <- filteredOutput[filteredOutput < maxAbPepConc]
filteredOutput <- filteredOutput[filteredOutput > minAbPepConc]
plot(density(filteredOutput))





simulationOutputDF <- as.data.frame(cbind(idxList, simulationOutput))
simulationOutputDFMelt <- melt(simulationOutputDF, measure.vars = c("AbPep"), id.vars = c("time"))
simulationOutputTimeMeans <- dcast(simulationOutputDFMelt, time ~ variable, fun.aggregate = mean, value.var ='value')
minConc <- min(simulationOutputTimeMeans$AbPep)

if(minConc < 0) {
  
  simulationOutputTimeMeans$AbPep <- simulationOutputTimeMeans$AbPep + abs(minConc) + simulationOutputTimeMeans[2,'AbPep']
  
}

#simulationOutput <- ode(y = concInit, times = timeSeries, func = kineticModel, parms = params, method = "lsoda")
#plot(x=simulationOutput[,"time"], y=simulationOutput[,"AbPep"], type="l", xlab = "Incubation Time (min)", ylab = "[AbPep]", main = "Simulated Two-State Ab-Peptide Binding vs. Incubation Time")
#plot(x=simulationOutputDF$time, y=simulationOutputDF$AbPep, type="l", xlab = "Incubation Time (min)", ylab = "[AbPep]", main = "Simulated Two-State Ab-Peptide Binding vs. Incubation Time")
plot(x=simulationOutputTimeMeans[,"time"], y=log10(simulationOutputTimeMeans[,"AbPep"]), type="l", xlab = "Incubation Time (min)", ylab = "[AbPep]", main = "Simulated Two-State Ab-Peptide Binding vs. Incubation Time")

currPlot <- ggplot(data=simulationOutputDF, aes(log10(simulationOutputDF$AbPep)), xmin=0) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab("[AbPep]") + ylab("Density") + 
  ggtitle(paste0("k1on = ", k1onParam, " k1off = ", k1offParam, " k2on = ", k2onParam, " k2off = ", k2offParam))
show(currPlot)


# Function to caculate the residuals between simulated and experimental data
# ssq=function(parms){
#   
#   # inital concentration
#   cinit=c(A=1,B=0,C=0)
#   # time points for which conc is reported
#   # include the points where data is available
#   t=c(seq(0,5,0.1),df$time)
#   t=sort(unique(t))
#   # parameters from the parameter estimation routine
#   k1=parms[1]
#   k2=parms[2]
#   # solve ODE for a given set of parameters
#   out=ode(y=cinit,times=t,func=rxnrate,parms=list(k1=k1,k2=k2))
#   
#   # Filter data that contains time points where data is available
#   outdf=data.frame(out)
#   outdf=outdf[outdf$time %in% df$time,]
#   # Evaluate predicted vs experimental residual
#   preddf=melt(outdf,id.var="time",variable.name="species",value.name="conc")
#   expdf=melt(df,id.var="time",variable.name="species",value.name="conc")
#   ssqres=preddf$conc-expdf$conc
#   
#   # return predicted vs experimental residual
#   return(ssqres)
#   
# }
# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
# parms=c(k1=0.5,k2=0.5)
# # fitting
# fitval=nls.lm(par=parms,fn=ssq)
# # Summary of fit
# summary(fitval)
# 
# Parameters:
#   Estimate  Std. Error  t value  Pr(>|t|)
# k1 2.01906  0.04867     41.49      <2e-16 ***
#   k2 0.99297  0.01779     55.82      <2e-16 ***
#   ---
#   Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0212 on 58 degrees of freedom
# Number of iterations to termination: 7
# Reason for termination: Relative error in the sum of squares is at most `ftol'.
#  
# # Estimated parameter
# parest=as.list(coef(fitval))
# parest
# $k1
# [1] 2.019065
#  
# $k2
# [1] 0.992973
#  
# # degrees of freedom: # data points - # parameters
# dof=3*nrow(df)-2
# dof
# [1] 58
# # mean error
# ms=sqrt(deviance(fitval)/dof)
# ms
# [1] 0.02119577
#  
# # variance Covariance Matrix
# S=vcov(fitval)
# S
#  
#       k1              k2
# k1 0.0023685244 -0.0003605831
# k2 -0.0003605831 0.0003164724
# plot of predicted vs experimental data

# simulated predicted profile at estimated parameter values
# cinit=c(A=1,B=0,C=0)
# t=seq(0,5,0.2)
# parms=as.list(parest)
# out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
# outdf=data.frame(out)
# names(outdf)=c("time","ca_pred","cb_pred","cc_pred")
# 
# # Overlay predicted profile with experimental data
# tmppred=melt(outdf,id.var=c("time"),variable.name="species",value.name="conc")
# tmpexp=melt(df,id.var=c("time"),variable.name="species",value.name="conc")
# p=ggplot(data=tmppred,aes(x=time,y=conc,color=species,linetype=species))+geom_line()
# p=p+geom_line(data=tmpexp,aes(x=time,y=conc,color=species,linetype=species))
# p=p+geom_point(data=tmpexp,aes(x=time,y=conc,color=species))
# p=p+scale_linetype_manual(values=c(0,1,0,1,0,1))
# p=p+scale_color_manual(values=rep(c("red","blue","green"),each=2))+theme_bw()
# print(p)



