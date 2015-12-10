library(limma)
library(ggplot2)
library(stats)
library(data.table)
library(openxlsx)
library(dplyr)
library(eVenn)

verbose_progress <- as.logical(TRUE)

# GPR Analysis Definition File
gprAnalysisDefFile <- '/Volumes/Macintosh HD/Users/mattg/Dropbox/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/04DEC2015/GPR_Data/Chagas_Panel_AB1_Specificity-12_min_10X.xlsx'

quantileCutoffs <- sort(as.numeric(c(1:40)), decreasing = TRUE)
calculateSpecificSeqs <- as.logical(FALSE)
maxQuantilesToReport <- 10

quantilePercentLabels <- as.numeric(c(1:length(quantileCutoffs))) * (100 / max(quantileCutoffs))

gprDataDef <- read.xlsx(gprAnalysisDefFile, sheet = 1, colNames = TRUE)

analysisGroups <- unique(gprDataDef$Analysis_Group)
sampleIndices <- c()
sampleIDs <- c()

for(currentAnalysisGroup in analysisGroups) {
  
  outputFile <- unique(gprDataDef$CSV_Output[gprDataDef$Analysis_Group == currentAnalysisGroup])
  dataSetNames <- unique(gprDataDef$Data_Set_Name[gprDataDef$Analysis_Group == currentAnalysisGroup])
  
  arrayOutputFile <- paste0(outputFile, '/Analysis_Group_', currentAnalysisGroup, '_Array_Metrics.csv')
  seqOutputFile <- paste0(outputFile, '/Analysis_Group_', currentAnalysisGroup, '_Seq_Metrics.csv')
  quantileMeansOutputFile <- paste0(outputFile, '/Analysis_Group_', currentAnalysisGroup, '_Quantile_Means.csv')
  quantileMediansOutputFile <- paste0(outputFile, '/Analysis_Group_', currentAnalysisGroup, '_Quantile_Medians.csv')
  quantileCountsOutputFile <- paste0(outputFile, '/Analysis_Group_', currentAnalysisGroup, '_Quantile_Counts.csv')
  
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
    
  quantileMeans <- matrix(data = NA, nrow = length(dataSetNames), ncol = length(quantileCutoffs))
  colnames(quantileMeans) <- quantileCutoffs
  rownames(quantileMeans) <- dataSetNames
  
  quantileRawMeans <- matrix(data = NA, nrow = length(dataSetNames), ncol = length(quantileCutoffs))
  colnames(quantileRawMeans) <- quantileCutoffs
  rownames(quantileRawMeans) <- dataSetNames
  
  quantileMedians <- matrix(data = NA, nrow = length(dataSetNames), ncol = length(quantileCutoffs))
  colnames(quantileMedians) <- quantileCutoffs
  rownames(quantileMedians) <- dataSetNames
  
  quantileRawMedians <- matrix(data = NA, nrow = length(dataSetNames), ncol = length(quantileCutoffs))
  colnames(quantileRawMedians) <- quantileCutoffs
  rownames(quantileRawMedians) <- dataSetNames
  
  quantileCounts <- matrix(data = NA, nrow = length(dataSetNames), ncol = length(quantileCutoffs))
  colnames(quantileCounts) <- quantileCutoffs
  rownames(quantileCounts) <- dataSetNames
  
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
  
  write.csv(gprDataZScoreScaled, file = arrayOutputFile, sep = "\t", col.names = TRUE, row.names = TRUE)
  write.csv(gprData, file = paste0(arrayOutputFile, "_unscaled.csv"), sep = "\t", col.names = TRUE, row.names = TRUE)
  
  gprDataZScoreScaledQuantiles <- as.data.frame(gprDataZScoreScaled, stringsAsFactors = FALSE)
  
  dataSetIdx <- 1
  
  for(currDataSet in dataSetNames) {
    
    gprDataZScoreScaledQuantiles[,dataSetIdx] <- ntile(gprDataZScoreScaledQuantiles[,dataSetIdx], length(quantileCutoffs))
    
    dataSetIdx <- dataSetIdx + 1
    
  }
  
  quantileBins <- quantileCutoffs[1:maxQuantilesToReport]
  
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
      
      quantileMeans[dataSetIdx,as.character(currQuantile)] <- mean(gprDataZScoreScaled[quantileIdxs,dataSetIdx])
      quantileRawMeans[dataSetIdx,as.character(currQuantile)] <- mean(gprData[quantileIdxs,dataSetIdx])
      
      quantileMedians[dataSetIdx,as.character(currQuantile)] <- median(gprDataZScoreScaled[quantileIdxs,dataSetIdx])
      quantileRawMedians[dataSetIdx,as.character(currQuantile)] <- median(gprData[quantileIdxs,dataSetIdx])
      
      quantileCounts[dataSetIdx,as.character(currQuantile)] <- length(gprDataZScoreScaled[quantileIdxs,dataSetIdx])
      
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
    
    evenn(matLists=currQuantileBinaryMatrix[which(rowSums(currQuantileBinaryMatrix) != 0),], prop=TRUE, annot = TRUE, display=TRUE, pathRes=comparativeOutputDir, CompName=paste0(currQuantile, '_Quantile_Bin'))
    
    if(calculateSpecificSeqs) {
      
      #******************!!!!!!!!!!!!*******************
      #vennDataPaths <- paste0(comparativeOutputDir, '/Venn_', currQuantile, '_Quantile_Bin/VennMatrixBin.txt')
      currVennPath <- paste0(comparativeOutputDir, '/Venn_', currQuantile, '_Quantile_Bin/VennMatrixBin.txt')
      
      vennCSVTable <- read.csv(currVennPath,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
      
      #vennCSVTableLogical <- cbind(sampleIDs[vennCSVTable[,1]], vennCSVTable)
      #vennCSVTableLogical <- cbind(gprFiles@.Data[[4]]$ID[vennCSVTable[,1]], vennCSVTable)
      vennCSVTableLogical <- vennCSVTable
      
      pairWiseColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
      pairWiseCombn <- combn(pairWiseColNames,2)
      pairWiseCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(pairWiseCombn))
      pairWiseCombnBinaryMatrixColNames <- c()
      maxSpecificOverlapCombnCount <- rep_len(1, ncol(pairWiseCombn))
      maxSpecificOverlapCombnCountColNames <- c()
      
      for(currCombn in 1:ncol(pairWiseCombn)) {
        
        pairWiseCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,pairWiseCombn[,currCombn]],1,function(x) {all(x)}))
        pairWiseCombnBinaryMatrixColNames <- cbind(pairWiseCombnBinaryMatrixColNames, paste0(pairWiseCombn[,currCombn], '_',collapse = ''))
        
      }
      
      colnames(pairWiseCombnBinaryMatrix) <- pairWiseCombnBinaryMatrixColNames
      maxSpecificOverlapCombnCountColNames <- append(maxSpecificOverlapCombnCountColNames, pairWiseCombnBinaryMatrixColNames)
      
      vennCSVTableLogical <- cbind(vennCSVTableLogical, pairWiseCombnBinaryMatrix)
      
      rm(pairWiseCombnBinaryMatrix)
      
      if(length(vennCSVTable)-3 >= 3) {
        
        maxSpecificOverlap <- 2
        
        triplicateWiseColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
        triplicateWiseCombn <- combn(triplicateWiseColNames,3)
        triplicateWiseCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(triplicateWiseCombn))
        triplicateWiseCombnBinaryMatrixColNames <- c()
        maxSpecificOverlapCombnCount <- append(maxSpecificOverlapCombnCount, rep_len(2, ncol(triplicateWiseCombn)))
        
        for(currCombn in 1:ncol(triplicateWiseCombn)) {
          
          triplicateWiseCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,triplicateWiseCombn[,currCombn]],1,function(x) {all(x)}))
          triplicateWiseCombnBinaryMatrixColNames <- cbind(triplicateWiseCombnBinaryMatrixColNames, paste0(triplicateWiseCombn[,currCombn], '_',collapse = ''))
          
        }
        
        colnames(triplicateWiseCombnBinaryMatrix) <- triplicateWiseCombnBinaryMatrixColNames
        maxSpecificOverlapCombnCountColNames <- append(maxSpecificOverlapCombnCountColNames, triplicateWiseCombnBinaryMatrixColNames)
        
        vennCSVTableLogical <- cbind(vennCSVTableLogical, triplicateWiseCombnBinaryMatrix)
        
        rm(triplicateWiseCombnBinaryMatrix)
        
      }
      
      if(length(vennCSVTable)-3 > 3) {
        
        singleVsRemainderColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
        singleVsRemainderCombn <- combn(singleVsRemainderColNames,(length(vennCSVTable)-3))
        singleVsRemainderCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(singleVsRemainderCombn))
        singleVsRemainderCombnBinaryMatrixColNames <- c()
        maxSpecificOverlapCombnCount <- append(maxSpecificOverlapCombnCount, rep_len((length(dataSetNames) - 1), ncol(singleVsRemainderCombn)))
      
        for(currCombn in 1:ncol(singleVsRemainderCombn)) {
          
          singleVsRemainderCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,singleVsRemainderCombn[,currCombn]],1,function(x) {all(x)}))
          singleVsRemainderCombnBinaryMatrixColNames <- cbind(singleVsRemainderCombnBinaryMatrixColNames, paste0(singleVsRemainderCombn[,currCombn], '_',collapse = ''))
          
        }
        
        colnames(singleVsRemainderCombnBinaryMatrix) <- singleVsRemainderCombnBinaryMatrixColNames
        maxSpecificOverlapCombnCountColNames <- append(maxSpecificOverlapCombnCountColNames, singleVsRemainderCombnBinaryMatrixColNames)
        
        vennCSVTableLogical <- cbind(vennCSVTableLogical, singleVsRemainderCombnBinaryMatrix)
        
        rm(singleVsRemainderCombnBinaryMatrix)
        
      }
      
      maxSpecificOverlapCombnCount <- t(maxSpecificOverlapCombnCount)
      colnames(maxSpecificOverlapCombnCount) <- maxSpecificOverlapCombnCountColNames
      
      #write.csv(vennCSVTableLogical, file = paste0(currVennPath, '_Quantile_', currQuantile, '_logicals.csv'), sep = "\t", col.names = TRUE, row.names = FALSE)
      
      # Calculate Specificity Totals for Each Sample Combination
      
      # In this Quantile for all samples
      allSamplesSeqs <- matrix(data = NA, nrow = nrow(vennCSVTableLogical), ncol = 1)
      currCombnSeqList <- as.vector(vennCSVTableLogical$X[which(vennCSVTableLogical$Total == length(dataSetNames))])
      allSamplesSeqs <- append(currCombnSeqList, matrix(data = NA, nrow = (nrow(vennCSVTableLogical) - length(currCombnSeqList)), ncol = 1))
      
      # In this Quantile for only specific samples
      specificSampleCols <- colnames(vennCSVTableLogical)[2:(length(dataSetNames)+1)]
      specificSampleCols <- specificSampleCols[which(!is.na(specificSampleCols))]
      
      specificSampleLists <- matrix(data = NA, nrow = nrow(vennCSVTableLogical), ncol = length(dataSetNames))
      colnames(specificSampleLists) <- specificSampleCols
      
      maxSpecificSampleCount <- 0
      
      for(currSampleCol in specificSampleCols) {
        
        currSampleSpecificIdxs <- as.logical(vennCSVTableLogical[,currSampleCol] == 1) & as.logical(rowSums(vennCSVTableLogical[,specificSampleCols]) == 1)
        
        currSampleSpecificSeqList <- as.vector(vennCSVTableLogical$X[currSampleSpecificIdxs])
        
        if(length(currSampleSpecificSeqList) > 0) {
          
          specificSampleLists[1:length(currSampleSpecificSeqList), currSampleCol] <- currSampleSpecificSeqList
          
          if(length(currSampleSpecificSeqList) > maxSpecificSampleCount) {
            
            maxSpecificSampleCount <- length(currSampleSpecificSeqList)
            
          }
          
        }
        
      }
      
      # In this Quantile for only specific sample combinations
      combnCols <- colnames(vennCSVTableLogical)[length(dataSetNames)+3:ncol(vennCSVTableLogical)]
      combnCols <- combnCols[which(!is.na(combnCols))]
      
      specificCombnLists <- matrix(data = NA, nrow = nrow(vennCSVTableLogical), ncol = length(combnCols))
      colnames(specificCombnLists) <- combnCols
      
      maxCombnCount <- 0
      
      for(currCombnCol in combnCols) {
        
        currCombnIdxs <- as.logical(vennCSVTableLogical[,currCombnCol] == 1) & as.logical(rowSums(vennCSVTableLogical[,combnCols]) <= maxSpecificOverlapCombnCount[1,currCombnCol])
        
        currCombnSeqList <- as.vector(vennCSVTableLogical$X[currCombnIdxs])
        
        if(length(currCombnSeqList) > 0) {
          
          specificCombnLists[1:length(currCombnSeqList), currCombnCol] <- currCombnSeqList
          
          if(length(currCombnSeqList) > maxCombnCount) {
            
            maxCombnCount <- length(currCombnSeqList)
            
          }
          
        }
        
      }
      
      
      specificCombnLists <- specificCombnLists[1:max(maxCombnCount,maxSpecificSampleCount),]
      specificSampleLists <- specificSampleLists[1:max(maxCombnCount,maxSpecificSampleCount),]
      
      specificSeqsList <- cbind(specificSampleLists, specificCombnLists)
      
      write.csv(specificSeqsList, file = paste0(currVennPath, '_Quantile_', currQuantile, '_specific_seqs.csv'), sep = "\t", col.names = TRUE, row.names = FALSE)
      
    }
    
  }
  
  colnames(quantileMeans) <- quantilePercentLabels
  colnames(quantileRawMeans) <- quantilePercentLabels
  colnames(quantileMedians) <- quantilePercentLabels
  colnames(quantileRawMedians) <- quantilePercentLabels
  
  write.csv(quantileMeans, file = quantileMeansOutputFile, sep = "\t", col.names = TRUE, row.names = TRUE)
  write.csv(quantileRawMeans, file = paste0(quantileMeansOutputFile, "_raw.csv"), sep = "\t", col.names = TRUE, row.names = TRUE)
  
  write.csv(quantileMedians, file = quantileMediansOutputFile, sep = "\t", col.names = TRUE, row.names = TRUE)
  write.csv(quantileRawMedians, file = paste0(quantileMediansOutputFile, "_raw.csv"), sep = "\t", col.names = TRUE, row.names = TRUE)
  
  write.csv(quantileCountsOutputFile, file = quantileMediansOutputFile, sep = "\t", col.names = TRUE, row.names = TRUE)
  
}