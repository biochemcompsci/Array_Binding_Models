# Read Sample Specific Lists for a Given Incubation and Wash Condition

library(ggplot2)

#workingDir <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Green_Channel_Chagas/60_min_incubation/10X_Wash/sample_level/"
#workingDir <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/12_min_incubation/10X_Wash/sample_level/"
workingDir <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/04DEC2015/GPR_Data/Green_Channel_Chagas/12_min_incubation/3X_Wash/sample_level/ver2/"
dataFilePrefix <- "Chagas_Twelve_Min_3X_Unscaled_Specific-"
#dataFilePrefix <- "Chagas_Sixty_Min_10X_Unscaled_Nonspecific-"
dataFileExtension <- ".csv"

sampleNames <- c("Chagas.G1.5969","Chagas.G1.6085","Chagas.Ctl.7276","Chagas.G2.6153","Chagas.ND53..Serum.","Chagas.Ctl.7300", "Chagas.G2.6066","Chagas.G1.6041","Chagas.G2.6051","Chagas.G2.6020","Chagas.G1.6133", "Chagas.Ctl.7278")

featureDF <- data.frame(stringsAsFactors = FALSE)
colLabels <- c()

filterSignals <- TRUE
filterColumn <- 2
filterFloorValue <- 250
filterCeilingValue <- 100000

#excludeIDs <- c('RHSVVGSG','RHSVVVGSG')
excludeIDs <- c()
idColumn <- 1

for(currSample in sampleNames) {
  
  currSpecificFeatures <- read.csv(paste0(workingDir,dataFilePrefix,currSample,dataFileExtension),header = TRUE, stringsAsFactors = FALSE)
  
  if(length(excludeIDs) > 0) {
    
    includeIdxs1 <- c(rep(1:nrow(currSpecificFeatures)))
    
    for(currExcludeID in excludeIDs) {
      
      includeIdxs1 <- includeIdxs1[which(!grepl(currExcludeID, currSpecificFeatures[includeIdxs1,idColumn]), arr.ind = TRUE)]
      
    }
    
    currSpecificFeatures <- currSpecificFeatures[includeIdxs1,]
    
  }
  
  if(filterSignals) {
    
    currSpecificFeatures <- currSpecificFeatures[which(currSpecificFeatures[,filterColumn] < filterCeilingValue),]
    currSpecificFeatures <- currSpecificFeatures[which(currSpecificFeatures[,filterColumn] > filterFloorValue),]
    
  }
  
  if(length(colLabels) == 0) {
    
    colLabels <- colnames(currSpecificFeatures)
    
  }
  
  currSpecificFeatures <- cbind(currSpecificFeatures, rep(currSample, times = nrow(currSpecificFeatures)))
  
  featureDF <- rbind(featureDF, currSpecificFeatures)
  
}

colnames(featureDF) <- c("Feature ID", colLabels[2:length(colLabels)], "Sample ID")

featureDF <- featureDF[which(!is.na(featureDF[,2])),]

write.csv(featureDF, paste0(workingDir,dataFilePrefix,"Combined_Specificity",dataFileExtension))

currPlot <- ggplot(data=featureDF, aes(featureDF[[2]]), xmin=0) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(dataFilePrefix)
ggsave(paste0(workingDir,dataFilePrefix,"Combined_Specificity.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)

currPlot <- ggplot(data=featureDF, aes(log10(featureDF[[2]])), xmin=0) +
  theme(axis.title = element_text(size = 5)) +
  geom_density(alpha=0.15, colour = 'black', fill = 'black') + xlab(dataFilePrefix)
ggsave(paste0(workingDir,dataFilePrefix,"Combined_Specificity_Log.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
