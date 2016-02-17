library(ggplot2)
library(reshape)

outputFile <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Green_Channel_Chagas/AB1_Twelve_100kDa_Non_Frac_10X-3X_Specific_Unscaled_Off-Rate_LogTF.csv"
#outputFile <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/12_min_incubation/Twelve_3X_10X_AB1_Unscaled_Off-Rate_LogTF.csv"

# zeroAplotAxis <- "0 - 12 Minute 10X Fluorescence Intensity per Min (Unscaled)"
# ABplotAxis <- "12 - 12 Minute 10X Chagas 100kDa Fractionated Non-Fractionated"
# BCplotAxis <- "60 - 120 Minute 10X Fluorescence Intensity per Min (Unscaled)"
zeroAplotAxis <- "0X - 3X AB1 Fluorescence Intensity per Min (Unscaled) 12 Min Incubation"
ABplotAxis <- "10X - 3X AB1 100kDa Non-Fractionated 12 Min Incubation"
BCplotAxis <- "0 - 120 Minute 3X Fluorescence Intensity per Min (Unscaled)"

# timePointFileA <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Green_Channel_Chagas/12_min_incubation/10X_Wash/sample_level/Chagas_Twelve_Min_10X_Unscaled-Combined_Specificity.csv"
# timePointFileB <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Green_Channel_Chagas/60_min_incubation/10X_Wash/sample_level/Chagas_Sixty_Min_10X_Unscaled-Combined_Specificity.csv"
# timePointFileC <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Green_Channel_Chagas/120_min_incubation/10X_Wash/sample_level/Chagas_OneTwenty_Min_10X_Unscaled-Combined_Specificity.csv"
timePointFileA <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/12_min_incubation/3X_Wash/sample_level/AB1_Twelve_Min_3X_Unscaled_Specific-Combined_Specificity.csv"
timePointFileB <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/12_min_incubation/10X_Wash/sample_level/AB1_Twelve_Min_10X_Unscaled_Specific-Combined_Specificity.csv"
timePointFileC <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/120_min_incubation/3X_Wash/sample_level/AB1_OneTwenty_Min_3X_Unscaled_Specific-Combined_Specificity.csv"

numTimePoints <- 2
logTransform <- TRUE

timeAMin <- 12
timeBMin <- 12
timeCMin <- 120

timeA <- read.csv(timePointFileA, header = TRUE, stringsAsFactors = FALSE)
timeB <- read.csv(timePointFileB, header = TRUE, stringsAsFactors = FALSE)

if(numTimePoints > 2){
  timeC <- read.csv(timePointFileC, header = TRUE, stringsAsFactors = FALSE)
}

time0 <- 0
scaleToMax <- FALSE

valueCol <- 'Z.Score.Mean'
#valueCol <- 'Raw.Mean'

dataCols <- c(2,3,4)
idColumn <- 1
signalColumn <- 2
saturatedSignal <- 65000
maxFeatures <- 10000

timeA <- timeA[,dataCols]
timeB <- timeB[,dataCols] 
  
if(numTimePoints > 2){
  
  timeC <- timeC[,dataCols]

}
# 
# timeA <- as.matrix.cast_matrix(t(cast(timeA, formula = ~ Feature.ID, mean, value = 'Z.Score.Mean')))
# timeB <- as.matrix.cast_matrix(t(cast(timeB, formula = ~ Feature.ID, mean, value = 'Z.Score.Mean')))
# timeC <- as.matrix.cast_matrix(t(cast(timeC, formula = ~ Feature.ID, mean, value = 'Z.Score.Mean')))
timeA <- as.matrix.cast_matrix(t(cast(timeA, formula = ~ Feature.ID, mean, value = as.character(valueCol))))
timeB <- as.matrix.cast_matrix(t(cast(timeB, formula = ~ Feature.ID, mean, value = as.character(valueCol))))
if(numTimePoints > 2){
  timeC <- as.matrix.cast_matrix(t(cast(timeC, formula = ~ Feature.ID, mean, value = as.character(valueCol))))
}

timeA <- as.data.frame(cbind(rownames(timeA),timeA[,1]), stringsAsFactors=FALSE, row.names = FALSE)
colnames(timeA) <- c("Feature.ID", "Signal")

timeB <- as.data.frame(cbind(rownames(timeB),timeB[,1]), stringsAsFactors=FALSE, row.names = FALSE)
colnames(timeB) <- c("Feature.ID", "Signal")

if(numTimePoints > 2){
  timeC <- as.data.frame(cbind(rownames(timeC),timeC[,1]), stringsAsFactors=FALSE, row.names = FALSE)
  colnames(timeC) <- c("Feature.ID", "Signal")
}

if(numTimePoints > 2){
  timeABCFeatureIDs <- intersect(timeC[,idColumn],intersect(timeB[,idColumn],timeA[,idColumn]))
} else {
  timeABCFeatureIDs <- intersect(timeB[,idColumn],timeA[,idColumn])
}
#timeABFeatureIDs <- timeA[timeA[,idColumn] %in% timeB[,idColumn], idColumn]
#timeABFeatureIDs <- timeA[,idColumn] 

if(length(timeABCFeatureIDs) > maxFeatures) {
  
  timeABCFeatureIDs <- timeABCFeatureIDs[sample(length(timeABCFeatureIDs), maxFeatures, replace = FALSE)]
  
}

timeABCSignalDiff <- c()
timeABCSignalRatio <- c()
featureCount <- 1                              
                              
for(currFeatureID in timeABCFeatureIDs) {
  
  cat(paste0(featureCount, " of ", length(timeABCFeatureIDs),"\n"))
  featureCount <- featureCount + 1
  
  currASignal <- as.double(timeA[timeA$Feature.ID == currFeatureID, signalColumn])
  
  if(length(currASignal) > 1) {
    
    currASignal <- mean(currASignal)
    
  }
  
  if(currASignal < 0) {
    
    currASignal <- 0
    
  }
  
  currASignalRate <- currASignal / timeAMin
  
  currBSignal <- as.double(timeB[timeB$Feature.ID == currFeatureID, signalColumn])
  
  if(length(currBSignal) > 1) {
    
    currBSignal <- mean(currBSignal)
    
  }
  
  if(currBSignal < 0) {
    
    currBSignal <- 0
    
  }
  
  currBSignalRate <- currBSignal / timeBMin
  
  if(numTimePoints > 2){
    
    currCSignal <- as.double(timeC[timeC$Feature.ID == currFeatureID, signalColumn])
    
    if(length(currCSignal) > 1) {
      
      currCSignal <- mean(currCSignal)
      
    }
    
    if(currCSignal < 0) {
      
      currCSignal <- 0
      
    }
    
    currCSignalRate <- currCSignal / timeCMin
  } else { 
    
    currCSignal <- 0
    currCSignalRate <- 0
    
  }
  
  currZeroASignalDiff <- diff(c(time0, currASignal))
    
  currABSignalDiff <- diff(c(currASignal,currBSignal))
  currBCSignalDiff <- diff(c(currBSignal,currCSignal))
  
  currBASignalRatio <- currBSignal/currASignal
  currCBSignalRatio <- currCSignal/currBSignal
  
  if(logTransform) {
    
    currSignalRow <- cbind(currFeatureID,log10(currASignal),log10(currBSignal), log10(currCSignal), log10(currASignalRate), log10(currBSignalRate),log10(currCSignalRate), log10(currZeroASignalDiff), log10(currABSignalDiff), log10(currBCSignalDiff), log10(currBASignalRatio), log10(currCBSignalRatio))
     
  } else {
    
    currSignalRow <- cbind(currFeatureID,currASignal,currBSignal, currCSignal, currASignalRate, currBSignalRate,currCSignalRate, currZeroASignalDiff, currABSignalDiff, currBCSignalDiff, currBASignalRatio, currCBSignalRatio)
    #currSignalRow <- cbind(currFeatureID,currASignal,currBSignal, "NULL", currASignalRate, currBSignalRate,"NULL", currZeroASignalDiff, currABSignalDiff, "NULL")
    
  }
  
  if(ncol(currSignalRow) == 12) {
    
    if(length(timeABCSignalDiff) == 0) {
      
      timeABCSignalDiff <- as.data.frame(currSignalRow,stringsAsFactors = FALSE)
      
    } else {
      
      timeABCSignalDiff <- rbind(timeABCSignalDiff,as.data.frame(currSignalRow,stringsAsFactors = FALSE))
       
    }
  }
}

# timeABCSignalDiff <- cbind(timeABCSignalDiff, 
#                            as.numeric(abs(as.numeric(timeABCSignalDiff[[5]])))/max(abs(as.numeric(timeABCSignalDiff[[5]])),na.rm = TRUE),
#                            as.numeric(abs(as.numeric(timeABCSignalDiff[[6]])))/max(abs(as.numeric(timeABCSignalDiff[[6]])),na.rm = TRUE))

if(scaleToMax) {

  if(logTransform) {
    
    timeABCSignalDiff <- cbind(timeABCSignalDiff, 
                               log10(as.numeric(timeABCSignalDiff[[5]]))/max(log10(as.numeric(timeABCSignalDiff[[5]])),na.rm = TRUE),
                               log10(as.numeric(timeABCSignalDiff[[6]]))/max(log10(as.numeric(timeABCSignalDiff[[6]])),na.rm = TRUE),
                               log10(as.numeric(timeABCSignalDiff[[7]]))/max(log10(as.numeric(timeABCSignalDiff[[7]])),na.rm = TRUE))
    
  } else {
    
    timeABCSignalDiff <- cbind(timeABCSignalDiff, 
                               as.numeric(timeABCSignalDiff[[5]])/max(as.numeric(timeABCSignalDiff[[5]]),na.rm = TRUE),
                               as.numeric(timeABCSignalDiff[[6]])/max(as.numeric(timeABCSignalDiff[[6]]),na.rm = TRUE),
                               as.numeric(timeABCSignalDiff[[7]])/max(as.numeric(timeABCSignalDiff[[7]]),na.rm = TRUE))
  }
  maxXAxis <- 1
  
} else {
  
  if(logTransform) {
    
    timeABCSignalDiff <- cbind(timeABCSignalDiff, 
                               log10(as.numeric(timeABCSignalDiff[[5]])),
                               log10(as.numeric(timeABCSignalDiff[[6]])),
                               log10(as.numeric(timeABCSignalDiff[[7]])))
  } else {
    
    timeABCSignalDiff <- cbind(timeABCSignalDiff, 
                               as.numeric(timeABCSignalDiff[[5]]),
                               as.numeric(timeABCSignalDiff[[6]]),
                               as.numeric(timeABCSignalDiff[[7]]))
    
  }
}

#colnames(timeABCSignalDiff) <- c('Feature.ID', 'Time.A.Signal', 'Time.B.Signal', 'Time.C.Signal', 'Time.A.Signal.Rate', 'Time.B.Signal.Rate', 'Time.C.Signal.Rate','Time.0.Time.A.Diff', 'Time.A.Time.B.Diff', 'Time.B.Time.C.Diff', 'Time.A.Signal.Rate.Scaled', 'Time.B.Signal.Rate.Scaled', 'Time.C.Signal.Rate.Scaled')
colnames(timeABCSignalDiff) <- c('Feature.ID', 'Time.A.Signal', 'Time.B.Signal', 'Time.C.Signal', 'Time.A.Signal.Rate', 'Time.B.Signal.Rate', 'Time.C.Signal.Rate','Time.0.Time.A.Diff', 'Time.A.Time.B.Diff', 'Time.B.Time.C.Diff', 'Time.B.Time.A.Ratio', 'Time.C.Time.B.Ratio', 'Time.A.Signal.Rate.Scaled', 'Time.B.Signal.Rate.Scaled', 'Time.C.Signal.Rate.Scaled')

write.csv(timeABCSignalDiff, outputFile)

allRates <- as.data.frame(c(as.numeric(t(timeABCSignalDiff$Time.A.Signal.Rate.Scaled)), as.numeric(t(timeABCSignalDiff$Time.B.Signal.Rate.Scaled)), as.numeric(t(timeABCSignalDiff$Time.C.Signal.Rate.Scaled))))
colnames(allRates) <- c('Signal.Rates')

maxXAxis <- ceiling(max(allRates$Signal.Rates,na.rm = TRUE))

currPlot <- ggplot(data=allRates, aes(allRates$Signal.Rates, xmin=0, xmax=ceiling(max(allRates$Signal.Rates,na.rm = TRUE)))) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(BCplotAxis)

ggsave(paste0(dirname(outputFile),"/All_Rates_Combined_dist.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)

currPlot <- ggplot(data=timeABCSignalDiff, aes(as.numeric(timeABCSignalDiff$Time.A.Signal.Rate.Scaled), xmin=0, xmax=ceiling(max(as.numeric(timeABCSignalDiff$Time.A.Signal.Rate.Scaled),na.rm = TRUE)))) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(zeroAplotAxis)

ggsave(paste0(dirname(outputFile),"/",zeroAplotAxis,"_dist.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)

currPlot <- ggplot(data=timeABCSignalDiff, aes(as.numeric(timeABCSignalDiff$Time.B.Signal.Rate.Scaled), xmin=0, xmax=ceiling(max(as.numeric(timeABCSignalDiff$Time.B.Signal.Rate.Scaled),na.rm = TRUE)))) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(ABplotAxis)

ggsave(paste0(dirname(outputFile),"/",ABplotAxis,"_dist.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)

currPlot <- ggplot(data=timeABCSignalDiff, aes(as.numeric(timeABCSignalDiff$Time.C.Signal.Rate.Scaled), xmin=0, xmax=ceiling(max(as.numeric(timeABCSignalDiff$Time.C.Signal.Rate.Scaled),na.rm = TRUE)))) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(BCplotAxis)

ggsave(paste0(dirname(outputFile),"/",BCplotAxis,"_dist.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)

currPlot <- ggplot(data=timeABCSignalDiff, aes(as.numeric(timeABCSignalDiff$Time.B.Time.A.Ratio), xmin=floor(min(as.numeric(timeABCSignalDiff$Time.B.Time.A.Ratio))), xmax=ceiling(max(as.numeric(timeABCSignalDiff$Time.B.Time.A.Ratio),na.rm = TRUE)))) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(paste0(ABplotAxis, ' Ratio'))

ggsave(paste0(dirname(outputFile),"/",ABplotAxis,"_dist_ratio.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)

currPlot <- ggplot(data=timeABCSignalDiff, aes(as.numeric(timeABCSignalDiff$Time.C.Time.B.Ratio), xmin=floor(min(as.numeric(timeABCSignalDiff$Time.B.Time.A.Ratio))), xmax=ceiling(max(as.numeric(timeABCSignalDiff$Time.C.Time.B.Ratio),na.rm = TRUE)))) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(paste0(BCplotAxis, ' Ratio'))

ggsave(paste0(dirname(outputFile),"/",BCplotAxis,"_dist_ratio.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
