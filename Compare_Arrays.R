library(limma)
library(ggplot2)
library(stats)
library(data.table)
library(openxlsx)
library(dplyr)
library(eVenn)

verbose_progress <- as.logical(TRUE)

# GPR Analysis Definition File
gprAnalysisDefFile <- '/Volumes/Macintosh HD/Users/mattg/Dropbox/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/Green_Channel_Chagas/Chagas_Panel_Condition_Specificity-2.xlsx'

quantileCutoffs <- sort(as.numeric(c(1:10)))

gprDataDef <- read.xlsx(gprAnalysisDefFile, sheet = 1, colNames = TRUE)

analysisGroups <- unique(gprDataDef$Analysis_Group)
sampleIndices <- c()
sampleIDs <- c()

for(currentAnalysisGroup in analysisGroups) {
  
  outputFile <- unique(gprDataDef$CSV_Output[gprDataDef$Analysis_Group == currentAnalysisGroup])
  dataSetNames <- unique(gprDataDef$Data_Set_Name[gprDataDef$Analysis_Group == currentAnalysisGroup])
  
  arrayOutputFile <- paste0(outputFile, 'Analysis_Group_', currentAnalysisGroup, '_Array_Metrics.csv')
  seqOutputFile <- paste0(outputFile, 'Analysis_Group_', currentAnalysisGroup, '_Seq_Metrics.csv')
  
  if(!dir.exists(paste0(outputFile, '/plots')) ) {
    
    dir.create(paste0(outputFile, '/plots/'))
      
  }
  
  comparativeOutputDir <- paste0(outputFile, '/plots')
  
  random_sample <- as.logical(unique(gprDataDef$Random_Sampling[gprDataDef$Analysis_Group == currentAnalysisGroup]))
  
  gprFiles <- read.maimages(files = unique(gprDataDef$GPR_Input_File[gprDataDef$Analysis_Group == currentAnalysisGroup]), 
                            source = 'genepix', columns = list(E=unique(gprDataDef$Signal_Field[gprDataDef$Analysis_Group == currentAnalysisGroup]), 
                                                               Eb=unique(gprDataDef$Background_Field[gprDataDef$Analysis_Group == currentAnalysisGroup])), 
                            annotation = list(unique(gprDataDef$Name_Field[gprDataDef$Analysis_Group == currentAnalysisGroup]), 
                                              unique(gprDataDef$ID_Field[gprDataDef$Analysis_Group == currentAnalysisGroup])))
  
  gprData <- data.frame(stringsAsFactors = FALSE)
  gprBackground <- data.frame(stringsAsFactors = FALSE)
  
  # ID(s) of Feature(s) of Interest
    
  sampleSize <-  unique(gprDataDef$Sample_Size[gprDataDef$Analysis_Group == currentAnalysisGroup])
  
  if(length(sampleIndices) == 0) {
    
    if(sampleSize == -1) {
      
      sampleIndices = c(1:length(gprFiles@.Data[[4]]$ID))
      
    } else {
    
      if(random_sample) {
      
        sampleIndices <- sample(1:length(gprFiles@.Data[[4]]$ID), sampleSize, replace = FALSE)
        
      } else {
       
        sampleIndices <- c(1:sampleSize)
         
      }
      
    }
    
    sampleIDs <- gprFiles@.Data[[4]]$ID[sampleIndices]
  
  }
  
  if(verbose_progress) {
    
    cat("Head Indices: ", head(sampleIndices), "\nTail Indices: ", tail(sampleIndices))
    
  }
    
  dataSetIdx <- 1
  
  for(currDataSet in dataSetNames) {
    
    if(length(gprData) == 0) {
    
      gprData <- data.frame(gprFiles@.Data[[1]][sampleIndices,dataSetIdx],stringsAsFactors = FALSE)
      gprBackground <- data.frame(gprFiles@.Data[[2]][sampleIndices,dataSetIdx],stringsAsFactors = FALSE)
      
    } else {
      
      gprData <- cbind(gprData, data.frame(gprFiles@.Data[[1]][sampleIndices,dataSetIdx],stringsAsFactors = FALSE))
      gprBackground <- cbind(gprBackground, data.frame(gprFiles@.Data[[2]][sampleIndices,dataSetIdx],stringsAsFactors = FALSE))
      
    }
    
    dataSetIdx <- dataSetIdx + 1
    
  }
  
  colnames(gprData) <- dataSetNames
  colnames(gprBackground) <- dataSetNames
  rownames(gprData) <- paste0(sampleIDs, '_-_', sampleIndices)
  rownames(gprBackground) <- paste0(sampleIDs, '_-_', sampleIndices)
  
  gprDataZScoreScaled <- scale(gprData)
  gprBackgroundZScoreScaled <- scale(gprBackground)
  
  rownames(gprDataZScoreScaled) <- paste0(sampleIDs, '_-_', sampleIndices)
  rownames(gprBackgroundZScoreScaled) <- paste0(sampleIDs, '_-_', sampleIndices)

  gprDataZScoreScaledQuantiles <- as.data.frame(gprDataZScoreScaled, stringsAsFactors = FALSE)
  
  dataSetIdx <- 1
  
  for(currDataSet in dataSetNames) {
    
    gprDataZScoreScaledQuantiles[,dataSetIdx] <- ntile(gprDataZScoreScaledQuantiles[,dataSetIdx], length(quantileCutoffs))
    
    dataSetIdx <- dataSetIdx + 1
    
  }
  
  quantileBins <- quantileCutoffs
  vennDataPaths <- c()
  vennFeatureIDs <- c()
  vennFeatureIdxs <- c()
  
  for(currQuantile in quantileBins) {
    
    dataSetIdx <- 1
    
    currQuantileBinaryMatrix <- matrix(0, nrow = nrow(gprDataZScoreScaledQuantiles), ncol = length(dataSetNames))
    colnames(currQuantileBinaryMatrix) <- dataSetNames
    #rownames(currQuantileBinaryMatrix) <- sampleIndices
    rownames(currQuantileBinaryMatrix) <- paste0(sampleIDs, '_-_', sampleIndices)
    
#     currQuantileBinaryMatrixIDs <- matrix(0, nrow = nrow(gprDataZScoreScaledQuantiles), ncol = length(dataSetNames))
#     colnames(currQuantileBinaryMatrixIDs) <- dataSetNames
#     rownames(currQuantileBinaryMatrixIDs) <- sampleIndices
#     
    for(currDataSet in dataSetNames) {
  
      # Venn Diagrams and Comparative Stats Output
      quantileIdxs <- which(gprDataZScoreScaledQuantiles[,dataSetIdx] == currQuantile)
      currQuantileBinaryMatrix[quantileIdxs,dataSetIdx] <- 1
#       currQuantileBinaryMatrixIDs[quantileIdxs,dataSetIdx] <- sampleIDs[quantileIdxs]
#       
      dataSetIdx <- dataSetIdx + 1
      
    }
    
#     nonZeroRows <- which(rowSums(currQuantileBinaryMatrix) > 0)
#     
#     currQuantileBinaryMatrixNonZero <- currQuantileBinaryMatrix[nonZeroRows,]
#     currQuantileBinaryMatrixIDsNonZero <- currQuantileBinaryMatrixIDs[nonZeroRows,]
#     
#    evenn(matLists=currQuantileBinaryMatrixNonZero, prop=TRUE, annot = TRUE, display=TRUE, pathRes=comparativeOutputDir, CompName=paste0(currQuantile, '_Quantile_Bin'))
    
    evenn(matLists=currQuantileBinaryMatrix, prop=TRUE, annot = TRUE, display=TRUE, pathRes=comparativeOutputDir, CompName=paste0(currQuantile, '_Quantile_Bin'))
    
    #******************!!!!!!!!!!!!*******************
    #vennDataPaths <- paste0(comparativeOutputDir, '/Venn_', currQuantile, '_Quantile_Bin/VennMatrixBin.txt')
    currVennPath <- paste0(comparativeOutputDir, '/Venn_', currQuantile, '_Quantile_Bin/VennMatrixBin.txt')
    
    vennCSVTable <- read.csv(currVennPath,sep = "\t",header = TRUE)
    
    #vennCSVTableLogical <- cbind(sampleIDs[vennCSVTable[,1]], vennCSVTable)
    #vennCSVTableLogical <- cbind(gprFiles@.Data[[4]]$ID[vennCSVTable[,1]], vennCSVTable)
    vennCSVTableLogical <- vennCSVTable
    
    pairWiseColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
    pairWiseCombn <- combn(pairWiseColNames,2)
    pairWiseCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(pairWiseCombn))
    pairWiseCombnBinaryMatrixColNames <- c()
    
    for(currCombn in 1:ncol(pairWiseCombn)) {
      
      pairWiseCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,pairWiseCombn[,currCombn]],1,function(x) {all(x)}))
      pairWiseCombnBinaryMatrixColNames <- cbind(pairWiseCombnBinaryMatrixColNames, paste0(pairWiseCombn[,currCombn], '_',collapse = ''))
      
    }
    
    colnames(pairWiseCombnBinaryMatrix) <- pairWiseCombnBinaryMatrixColNames
    
    vennCSVTableLogical <- cbind(vennCSVTableLogical, pairWiseCombnBinaryMatrix)
    
    rm(pairWiseCombnBinaryMatrix)
    
    if(length(vennCSVTable)-3 >= 3) {
      
      triplicateWiseColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
      triplicateWiseCombn <- combn(triplicateWiseColNames,3)
      triplicateWiseCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(triplicateWiseCombn))
      triplicateWiseCombnBinaryMatrixColNames <- c()
      
      for(currCombn in 1:ncol(triplicateWiseCombn)) {
        
        triplicateWiseCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,triplicateWiseCombn[,currCombn]],1,function(x) {all(x)}))
        triplicateWiseCombnBinaryMatrixColNames <- cbind(triplicateWiseCombnBinaryMatrixColNames, paste0(triplicateWiseCombn[,currCombn], '_',collapse = ''))
        
      }
      
      colnames(triplicateWiseCombnBinaryMatrix) <- triplicateWiseCombnBinaryMatrixColNames
      
      vennCSVTableLogical <- cbind(vennCSVTableLogical, triplicateWiseCombnBinaryMatrix)
      
      rm(triplicateWiseCombnBinaryMatrix)
      
    }
    
    if(length(vennCSVTable)-3 > 3) {
      
      singleVsRemainderColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
      singleVsRemainderCombn <- combn(singleVsRemainderColNames,(length(vennCSVTable)-3))
      singleVsRemainderCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(singleVsRemainderCombn))
      singleVsRemainderCombnBinaryMatrixColNames <- c()
      
      for(currCombn in 1:ncol(singleVsRemainderCombn)) {
        
        singleVsRemainderCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,singleVsRemainderCombn[,currCombn]],1,function(x) {all(x)}))
        singleVsRemainderCombnBinaryMatrixColNames <- cbind(singleVsRemainderCombnBinaryMatrixColNames, paste0(singleVsRemainderCombn[,currCombn], '_',collapse = ''))
        
      }
      
      colnames(singleVsRemainderCombnBinaryMatrix) <- singleVsRemainderCombnBinaryMatrixColNames
      
      vennCSVTableLogical <- cbind(vennCSVTableLogical, singleVsRemainderCombnBinaryMatrix)
      
      rm(singleVsRemainderCombnBinaryMatrix)
      
    }
    
    write.csv(vennCSVTableLogical[which(vennCSVTableLogical['Total'] != 0),], file = paste0(currVennPath, '_Quantile_', currQuantile, '_logicals.csv'), sep = "\t", col.names = TRUE, row.names = FALSE)
    
  }

}