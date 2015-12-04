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

# quantileCutoff_A <- 0.05
# quantileCutoff_B <- 0.25
# quantileCutoff_C <- 0.75
# quantileCutoff_D <- 0.95
# 
# quantileCutoffs <- sort(as.numeric(c(quantileCutoff_A, quantileCutoff_B, quantileCutoff_C, quantileCutoff_D)))
quantileCutoffs <- sort(as.numeric(c(1:10)))

gprDataDef <- read.xlsx(gprAnalysisDefFile, sheet = 1, colNames = TRUE)

#dataSetName_A <- '815_S4'
#dataSetName_B <- '815_S5'

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
  
  #sampleIndices <- c()
  
  gprDataZScoreScaled <- scale(gprData)
  gprBackgroundZScoreScaled <- scale(gprBackground)
  
  rownames(gprDataZScoreScaled) <- sampleIndices
  rownames(gprBackgroundZScoreScaled) <- sampleIndices
  
#  gprGSGIdx <- data.frame('GSG_Index' = which(gprFiles@.Data[[4]]$ID == 'GSG'), stringsAsFactors = FALSE)
  
#   array_A_Quantiles <- quantile(gprFiles@.Data[[1]][,1], quantileCutoffs)
#   array_B_Quantiles <- quantile(gprFiles@.Data[[1]][,2], quantileCutoffs)
#   
  
  gprDataZScoreScaledQuantiles <- as.data.frame(gprDataZScoreScaled, stringsAsFactors = FALSE)
  
  dataSetIdx <- 1
  
  for(currDataSet in dataSetNames) {
    
    gprDataZScoreScaledQuantiles[,dataSetIdx] <- ntile(gprDataZScoreScaledQuantiles[,dataSetIdx], length(quantileCutoffs))
    
    dataSetIdx <- dataSetIdx + 1
    
  }
  
  #quantileBins <- sort(unique(gprDataZScoreScaledQuantiles[,dataSetIdx]))
  quantileBins <- quantileCutoffs
  vennDataPaths <- c()
  vennFeatureIDs <- c()
  vennFeatureIdxs <- c()
  
  for(currQuantile in quantileBins) {
    
    dataSetIdx <- 1
    
    currQuantileBinaryMatrix <- matrix(0, nrow = nrow(gprDataZScoreScaledQuantiles), ncol = length(dataSetNames))
    colnames(currQuantileBinaryMatrix) <- dataSetNames
    
    currQuantileBinaryMatrixIDs <- matrix(0, nrow = nrow(gprDataZScoreScaledQuantiles), ncol = length(dataSetNames))
    colnames(currQuantileBinaryMatrixIDs) <- dataSetNames
    
    for(currDataSet in dataSetNames) {
  
      # Venn Diagrams and Comparative Stats Output
      quantileIdxs <- which(gprDataZScoreScaledQuantiles[,dataSetIdx] == currQuantile)
      currQuantileBinaryMatrix[quantileIdxs,dataSetIdx] <- 1
      currQuantileBinaryMatrixIDs[quantileIdxs,dataSetIdx] <- sampleIDs[quantileIdxs]
      
      dataSetIdx <- dataSetIdx + 1
      
#       q1ZScoreB <- which(gprDataZScoreScaled[,2] >= array_B_ZScore_Quantiles[4])
#       
#       q1ZScoreValsA <- cbind(q1ZScoreA, gprDataZScoreScaled[q1ZScoreA,1])
#       q1ZScoreValsB <- cbind(q1ZScoreB, gprDataZScoreScaled[q1ZScoreB,2])
#       
#       q1ABZScoreIdxs <- sort(unique(as.numeric(rbind(q1ZScoreA, q1ZScoreB))))
#       listLength <- length(q1ABZScoreIdxs)
#       
#       q1ABZScoreBinaryMat <- as.matrix(cbind(as.vector(rep(0, listLength)), as.vector(rep(0, listLength))))
#       colnames(q1ABZScoreBinaryMat) <- c(dataSetName_A, dataSetName_B)
#       rownames(q1ABZScoreBinaryMat) <- q1ABZScoreIdxs
#       q1ABZScoreBinaryMat[as.numeric(rownames(q1ABZScoreBinaryMat)) %in% q1ZScoreA, 1] <- 1
#       q1ABZScoreBinaryMat[as.numeric(rownames(q1ABZScoreBinaryMat)) %in% q1ZScoreB, 2] <- 1
      
    }
    
    nonZeroRows <- which(rowSums(currQuantileBinaryMatrix) > 0)
    
    currQuantileBinaryMatrixNonZero <- currQuantileBinaryMatrix[nonZeroRows,]
    currQuantileBinaryMatrixIDsNonZero <- currQuantileBinaryMatrixIDs[nonZeroRows,]
    
    rownames(currQuantileBinaryMatrixNonZero) <- sampleIndices[nonZeroRows]
    colnames(currQuantileBinaryMatrixNonZero) <- colnames(currQuantileBinaryMatrix)
    #names(currQuantileBinaryMatrix) <- currQuantileBinaryMatrixRows
    
    rownames(currQuantileBinaryMatrixIDsNonZero) <- sampleIndices[nonZeroRows]
    colnames(currQuantileBinaryMatrixIDsNonZero) <- colnames(currQuantileBinaryMatrix)
    #names(currQuantileBinaryMatrixIDs) <- currQuantileBinaryMatrixRows
    
    evenn(matLists=currQuantileBinaryMatrixNonZero, prop=TRUE, annot = TRUE, display=TRUE, pathRes=comparativeOutputDir, CompName=paste0(currQuantile, '_Quantile_Bin'))
    
    #******************!!!!!!!!!!!!*******************
    vennDataPaths <- paste0(comparativeOutputDir, '/Venn_', currQuantile, '_Quantile_Bin/VennMatrixBin.txt')
    currVennPath <- vennDataPaths
    vennCSVTable <- read.csv(currVennPath,sep = "\t",header = TRUE)
    
    #vennCSVTableLogical <- cbind(sampleIDs[vennCSVTable[,1]], vennCSVTable)
    vennCSVTableLogical <- cbind(gprFiles@.Data[[4]]$ID[vennCSVTable[,1]], vennCSVTable)
      
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
    
    write.csv(vennCSVTableLogical, file = paste0(currVennPath,'_logicals.csv'), sep = "\t", col.names = TRUE, row.names = FALSE)
    
    #******************!!!!!!!!!!!!*******************
    #vennDataPaths <- cbind(vennDataPaths, paste0(comparativeOutputDir, '/Venn_', currQuantile, '_Quantile_Bin/VennMatrixBin.txt'))
    
    
    #vennFeatureIDs <- cbind(vennFeatureIDs, sampleIDs[currQuantileBinaryMatrixRows])
    #vennFeatureIdxs <- cbind(vennFeatureIdxs, currQuantileBinaryMatrixRows)
    
  }
  
#   currVennCount <- 1
#   
#   for(currVennPath in vennDataPaths) {
#     
#     vennCSVTable <- read.csv(currVennPath,sep = "\t",header = TRUE)
#     
#     vennCSVTableLogical <- cbind(vennFeatureIDs[,currVennCount], vennCSVTable)
#     currVennCount <- currVennCount + 1
#     
#     pairWiseColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
#     pairWiseCombn <- combn(pairWiseColNames,2)
#     pairWiseCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(pairWiseCombn))
#     pairWiseCombnBinaryMatrixColNames <- c()
#     
#     for(currCombn in 1:ncol(pairWiseCombn)) {
#       
#       pairWiseCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,pairWiseCombn[,currCombn]],1,function(x) {all(x)}))
#       pairWiseCombnBinaryMatrixColNames <- cbind(pairWiseCombnBinaryMatrixColNames, paste0(pairWiseCombn[,currCombn], '_',collapse = ''))
#       
#     }
#     
#     colnames(pairWiseCombnBinaryMatrix) <- pairWiseCombnBinaryMatrixColNames
#     
#     vennCSVTableLogical <- cbind(vennCSVTableLogical, pairWiseCombnBinaryMatrix)
#     
#     rm(pairWiseCombnBinaryMatrix)
#     
#     if(length(vennCSVTable)-3 >= 3) {
#       
#       triplicateWiseColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
#       triplicateWiseCombn <- combn(triplicateWiseColNames,3)
#       triplicateWiseCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(triplicateWiseCombn))
#       triplicateWiseCombnBinaryMatrixColNames <- c()
#       
#       for(currCombn in 1:ncol(triplicateWiseCombn)) {
#         
#         triplicateWiseCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,triplicateWiseCombn[,currCombn]],1,function(x) {all(x)}))
#         triplicateWiseCombnBinaryMatrixColNames <- cbind(triplicateWiseCombnBinaryMatrixColNames, paste0(triplicateWiseCombn[,currCombn], '_',collapse = ''))
#         
#       }
#       
#       colnames(triplicateWiseCombnBinaryMatrix) <- triplicateWiseCombnBinaryMatrixColNames
#       
#       vennCSVTableLogical <- cbind(vennCSVTableLogical, triplicateWiseCombnBinaryMatrix)
#       
#       rm(triplicateWiseCombnBinaryMatrix)
#       
#     }
#     
#     if(length(vennCSVTable)-3 > 3) {
#       
#       singleVsRemainderColNames <- colnames(vennCSVTable)[2:(length(vennCSVTable)-1)]
#       singleVsRemainderCombn <- combn(singleVsRemainderColNames,(length(vennCSVTable)-3))
#       singleVsRemainderCombnBinaryMatrix <- matrix(0, nrow = nrow(vennCSVTable), ncol = ncol(singleVsRemainderCombn))
#       singleVsRemainderCombnBinaryMatrixColNames <- c()
#       
#       for(currCombn in 1:ncol(singleVsRemainderCombn)) {
#         
#         singleVsRemainderCombnBinaryMatrix[,currCombn] <- as.numeric(apply(vennCSVTable[,singleVsRemainderCombn[,currCombn]],1,function(x) {all(x)}))
#         singleVsRemainderCombnBinaryMatrixColNames <- cbind(singleVsRemainderCombnBinaryMatrixColNames, paste0(singleVsRemainderCombn[,currCombn], '_',collapse = ''))
#         
#       }
#       
#       colnames(singleVsRemainderCombnBinaryMatrix) <- singleVsRemainderCombnBinaryMatrixColNames
#     
#       vennCSVTableLogical <- cbind(vennCSVTableLogical, singleVsRemainderCombnBinaryMatrix)
#       
#       rm(singleVsRemainderCombnBinaryMatrix)
#       
#     }
#     
#     write.csv(vennCSVTableLogical, file = paste0(currVennPath,'_logicals.csv'), sep = "\t", col.names = TRUE, row.names = FALSE)
#     
#     rm(vennCSVTable)
#     
#   }
  
  #evenn(matLists=q1ABZScoreBinaryMat, prop=TRUE, annot = TRUE, display=TRUE, pathRes=comparativeOutputDir, CompName=paste0(dataSetName_A, '_vs_', dataSetName_B, '_', '95th_Percentile'))
  
#   q1A <- which(gprData$Array_A <= array_A_Quantiles[1])
#   q1B <- which(gprData$Array_B <= array_B_Quantiles[1])
#   q1A_GSG <- gprData$Array_A[intersect(gprGSGIdx$GSG_Index, which(gprData$Array_A <= array_A_Quantiles[1]))]
#   q1B_GSG <- gprData$Array_B[intersect(gprGSGIdx$GSG_Index, which(gprData$Array_B <= array_B_Quantiles[1]))]
#   q1Intersect <- intersect(q1A,q1B)
#   q1DiffSet <- setdiff(q1A,q1B)
#   q1IntersectDiffRatio <- length(q1Intersect) / (length(q1Intersect) + length(q1DiffSet))  
#   
#   q2A <- which(between(gprData$Array_A, array_A_Quantiles[1], array_A_Quantiles[2], incbounds=TRUE))
#   q2B <- which(between(gprData$Array_B, array_B_Quantiles[1], array_B_Quantiles[2], incbounds=TRUE))
#   q2A_GSG <- gprData$Array_A[intersect(gprGSGIdx$GSG_Index, which(between(gprData$Array_A, array_A_Quantiles[1], array_A_Quantiles[2], incbounds=TRUE)))]
#   q2B_GSG <- gprData$Array_B[intersect(gprGSGIdx$GSG_Index, which(between(gprData$Array_B, array_B_Quantiles[1], array_B_Quantiles[2], incbounds=TRUE)))]
#   q2Intersect <- intersect(q2A,q2B)
#   q2DiffSet <- setdiff(q2A,q2B)
#   q2IntersectDiffRatio <- length(q2Intersect) / (length(q2Intersect) + length(q2DiffSet))  
#   
#   q3A <- which(between(gprData$Array_A, array_A_Quantiles[2], array_A_Quantiles[3], incbounds=TRUE))
#   q3B <- which(between(gprData$Array_B, array_B_Quantiles[2], array_B_Quantiles[3], incbounds=TRUE))
#   q3A_GSG <- gprData$Array_A[intersect(gprGSGIdx$GSG_Index, which(between(gprData$Array_A, array_A_Quantiles[2], array_A_Quantiles[3], incbounds=TRUE)))]
#   q3B_GSG <- gprData$Array_B[intersect(gprGSGIdx$GSG_Index, which(between(gprData$Array_B, array_B_Quantiles[2], array_B_Quantiles[3], incbounds=TRUE)))]
#   q3Intersect <- intersect(q3A,q3B)
#   q3DiffSet <- setdiff(q3A,q3B)
#   q3IntersectDiffRatio <- length(q3Intersect) / (length(q3Intersect) + length(q3DiffSet)) 
#   
#   q4A <- which(between(gprData$Array_A, array_A_Quantiles[3], array_A_Quantiles[4], incbounds=TRUE))
#   q4B <- which(between(gprData$Array_B, array_B_Quantiles[3], array_B_Quantiles[4], incbounds=TRUE))
#   q4A_GSG <- gprData$Array_A[intersect(gprGSGIdx$GSG_Index, which(between(gprData$Array_A, array_A_Quantiles[3], array_A_Quantiles[4], incbounds=TRUE)))]
#   q4B_GSG <- gprData$Array_B[intersect(gprGSGIdx$GSG_Index, which(between(gprData$Array_B, array_B_Quantiles[3], array_B_Quantiles[4], incbounds=TRUE)))]
#   q4Intersect <- intersect(q4A,q4B)
#   q4DiffSet <- setdiff(q4A,q4B)
#   q4IntersectDiffRatio <- length(q4Intersect) / (length(q4Intersect) + length(q4DiffSet)) 
#   
#   q5A <- which(gprData$Array_A >= array_A_Quantiles[4])
#   q5B <- which(gprData$Array_B >= array_B_Quantiles[4])
#   q5A_GSG <- gprData$Array_A[intersect(gprGSGIdx$GSG_Index, which(gprData$Array_A >= array_A_Quantiles[4]))]
#   q5B_GSG <- gprData$Array_B[intersect(gprGSGIdx$GSG_Index, which(gprData$Array_B >= array_B_Quantiles[4]))]
#   q5Intersect <- intersect(q5A,q5B)
#   q5DiffSet <- setdiff(q5A,q5B)
#   q5IntersectDiffRatio <- length(q5Intersect) / (length(q5Intersect) + length(q5DiffSet))
#   
#   array_Comparison_Summary_A <- data.frame(median(gprData$Array_A),
#                                            mean(gprData$Array_A),
#                                            min(gprData$Array_A),
#                                            max(gprData$Array_A),
#                                            median(gprBackground$Array_A),
#                                            mean(gprBackground$Array_A),
#                                            median(gprData$Array_A[gprGSGIdx$GSG_Index]),
#                                            mean(gprData$Array_A[gprGSGIdx$GSG_Index]),
#                                            min(gprData$Array_A[gprGSGIdx$GSG_Index]),
#                                            max(gprData$Array_A[gprGSGIdx$GSG_Index]),
#                                            length(which(gprData$Array_A > 65500)),
#                                            array_A_Quantiles[1],
#                                            array_A_Quantiles[2],
#                                            array_A_Quantiles[3],
#                                            array_A_Quantiles[4],
#                                            length(q1A),
#                                            length(q2A),
#                                            length(q3A),
#                                            length(q4A),
#                                            length(q5A),
#                                            length(q1A_GSG),
#                                            length(q2A_GSG),
#                                            length(q3A_GSG),
#                                            length(q4A_GSG),
#                                            length(q5A_GSG),
#                                            stringsAsFactors = FALSE, row.names = c(dataSetName_A))
#   
#   array_Comparison_Summary_B <- data.frame(median(gprData$Array_B),
#                                            mean(gprData$Array_B),
#                                            min(gprData$Array_B),
#                                            max(gprData$Array_B),
#                                            median(gprBackground$Array_B),
#                                            mean(gprBackground$Array_B),
#                                            median(gprData$Array_B[gprGSGIdx$GSG_Index]),
#                                            mean(gprData$Array_B[gprGSGIdx$GSG_Index]),
#                                            min(gprData$Array_B[gprGSGIdx$GSG_Index]),
#                                            max(gprData$Array_B[gprGSGIdx$GSG_Index]),
#                                            length(which(gprData$Array_B > 65500)),
#                                            array_B_Quantiles[1],
#                                            array_B_Quantiles[2],
#                                            array_B_Quantiles[3],
#                                            array_B_Quantiles[4], 
#                                            length(q1B),
#                                            length(q2B),
#                                            length(q3B),
#                                            length(q4B),
#                                            length(q5B),
#                                            length(q1B_GSG),
#                                            length(q2B_GSG),
#                                            length(q3B_GSG),
#                                            length(q4B_GSG),
#                                            length(q5B_GSG),
#                                            stringsAsFactors = FALSE, row.names = c(dataSetName_B))
#   
#   colHeaders <- as.character(c('(Intensity) Median', 
#                                '(Intensity) Mean',
#                                '(Intensity) Min',
#                                '(Intensity) Max',
#                                '(Background) Median', 
#                                '(Background) Mean',
#                                '(GSG) Median',
#                                '(GSG) Mean',
#                                '(GSG) Min',
#                                '(GSG) Max',
#                                '(Count) Saturated',
#                                paste0('(Intensity) ',quantileCutoff_A, ' Quantile'), 
#                                paste0('(Intensity) ',quantileCutoff_B, ' Quantile'), 
#                                paste0('(Intensity) ',quantileCutoff_C, ' Quantile'), 
#                                paste0('(Intensity) ',quantileCutoff_D, ' Quantile'),
#                                paste0('(Count) <= ',quantileCutoff_A, ' Quantile'), 
#                                paste0('(Count) <> ',quantileCutoff_A, ' - ',quantileCutoff_B, ' Quantile'), 
#                                paste0('(Count) <> ',quantileCutoff_B, ' - ',quantileCutoff_C, ' Quantile'), 
#                                paste0('(Count) <> ',quantileCutoff_C, ' - ',quantileCutoff_D, ' Quantile'),
#                                paste0('(Count) >= ',quantileCutoff_D, ' Quantile'),
#                                paste0('(GSG Count) <= ',quantileCutoff_A, ' Quantile'), 
#                                paste0('(GSG Count) <> ',quantileCutoff_A, ' - ',quantileCutoff_B, ' Quantile'), 
#                                paste0('(GSG Count) <> ',quantileCutoff_B, ' - ',quantileCutoff_C, ' Quantile'), 
#                                paste0('(GSG Count) <> ',quantileCutoff_C, ' - ',quantileCutoff_D, ' Quantile'),
#                                paste0('(GSG Count) >= ',quantileCutoff_D, ' Quantile')))
#   
#   names(array_Comparison_Summary_A) <- colHeaders
#   names(array_Comparison_Summary_B) <- colHeaders
#   
#   array_Comparison_Summary <- rbind(array_Comparison_Summary_A,array_Comparison_Summary_B)
#   
#   array_Diff_Summary <- data.frame(length(q1Intersect),
#                                    length(q2Intersect),
#                                    length(q3Intersect),
#                                    length(q4Intersect),
#                                    length(q5Intersect),
#                                    length(q1DiffSet),
#                                    length(q2DiffSet),
#                                    length(q3DiffSet),
#                                    length(q4DiffSet),
#                                    length(q5DiffSet), stringsAsFactors = FALSE, row.names = NULL)
#   
#   colHeaders <- as.character(c(paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Intersection Count) <= ',quantileCutoff_A,' Quantile'), 
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Intersection Count) <> ',quantileCutoff_A, ' - ',quantileCutoff_B, ' Quantile'), 
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Intersection Count) <> ',quantileCutoff_B, ' - ',quantileCutoff_C, ' Quantile'), 
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Intersection Count) <> ',quantileCutoff_C, ' - ',quantileCutoff_D, ' Quantile'),
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Intersection Count) >= ',quantileCutoff_D,' Quantile'),
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Difference Count) <= ',quantileCutoff_A,' Quantile'), 
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Difference Count) <> ',quantileCutoff_A, ' - ',quantileCutoff_B, ' Quantile'), 
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Difference Count) <> ',quantileCutoff_B, ' - ',quantileCutoff_C, ' Quantile'), 
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Difference Count) <> ',quantileCutoff_C, ' - ',quantileCutoff_D, ' Quantile'),
#                                            paste0(dataSetName_A, ' vs ', dataSetName_B, ' (Difference Count) >= ',quantileCutoff_D,' Quantile')))
#   
#   names(array_Diff_Summary) <- colHeaders
#   
#   array_Comparison_Summary <- cbind(array_Comparison_Summary, array_Diff_Summary)
#   
#   write.csv(array_Comparison_Summary, file = arrayOutputFile)
#   
#   co.var.df <- function(x) ( apply(x,1,sd)/rowMeans(x) )
#   
#   sequenceMetrics <- data.frame('ID' = gprFiles@.Data[[4]]$ID, 
#                                 'Quantile Array 1' = rep(0,nrow(gprFiles@.Data[[4]])), 
#                                 'Quantile Array 2' = rep(0,nrow(gprFiles@.Data[[4]])), 
#                                 'FG Median Array 1' = gprData$Array_A,
#                                 'FG Median Array 2' = gprData$Array_B,
#                                 'BG Median Array 1' = gprBackground$Array_A,
#                                 'BG Median Array 2' = gprBackground$Array_B, 
#                                 'Array 1 vs 2 FG Median CV' =  co.var.df(cbind(gprData$Array_A,gprData$Array_B)))
#   
#   sequenceMetrics$Quantile.Array.1[q1A] <- paste0('<= ',quantileCutoff_A)
#   sequenceMetrics$Quantile.Array.1[q2A] <- paste0('<> ',quantileCutoff_A, ' - ',quantileCutoff_B)
#   sequenceMetrics$Quantile.Array.1[q3A] <- paste0('<> ',quantileCutoff_B, ' - ',quantileCutoff_C)  
#   sequenceMetrics$Quantile.Array.1[q4A] <- paste0('<> ',quantileCutoff_C, ' - ',quantileCutoff_D) 
#   sequenceMetrics$Quantile.Array.1[q5A] <- paste0('>= ',quantileCutoff_D)
#   sequenceMetrics$Quantile.Array.2[q1B] <- paste0('<= ',quantileCutoff_A)
#   sequenceMetrics$Quantile.Array.2[q2B] <- paste0('<> ',quantileCutoff_A, ' - ',quantileCutoff_B)
#   sequenceMetrics$Quantile.Array.2[q3B] <- paste0('<> ',quantileCutoff_B, ' - ',quantileCutoff_C)  
#   sequenceMetrics$Quantile.Array.2[q4B] <- paste0('<> ',quantileCutoff_C, ' - ',quantileCutoff_D) 
#   sequenceMetrics$Quantile.Array.2[q5B] <- paste0('>= ',quantileCutoff_D)
#   
#   write.csv(sequenceMetrics, file = seqOutputFile)        
  
#   xyScatterPlotFull <- ggplot(gprData, aes(x = Array_A, y = Array_B)) + geom_point() + ggtitle("All Feature Scatter 691-13 vs 692-13, ND53 Sample")
#   show(xyScatterPlotFull)
  
}