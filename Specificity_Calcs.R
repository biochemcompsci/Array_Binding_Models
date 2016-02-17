library(Biobase)

outputDir <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/04DEC2015/GPR_Data/Green_Channel_Chagas/12_min_incubation/10X_Wash/sample_level/ver2/"
#outputDir <- "~/Dropbox (Nextval)/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/120_min_incubation/10X_Wash/sample_level/"
inputFile <- "Analysis_Group_1_Array_Metrics.csv_unscaled.csv"
#inputFile <- "Analysis_Group_1_Array_Metrics.csv_unscaled.csv"
#sampleCondition <- "Chagas_OneTwenty_Min_10X_Unscaled"
sampleCondition <- "Chagas_Twelve_Min_10X_Unscaled_Specific"

#outputColNames <- c('Feature ID (Mean)', 'Raw Mean', 'Specificity Mean', 'Feature ID (Median)', 'Raw Median', 'Specificity Median')
outputColNames <- c('Feature ID (Mean)', 'Z-Score Mean', 'Specificity Mean', 'Feature ID (Median)', 'Z-Score Median', 'Specificity Median')

zScoreCols <- c(rep(2:13))

Z_Cutoff <- 200
#Z_Cutoff <- 2
Specificity_Cutoff <- 2
Calculate_Nonspecific <- FALSE

currData <- read.csv(paste0(outputDir, inputFile), header = TRUE, stringsAsFactors=FALSE)

#analysisCols <- c('Chagas.G2.6153')
#analysisCols <- c('Chagas.Ctl.7278')
#analysisCols <- c('Chagas.ND53..Serum.')
analysisCols <- colnames(currData[,zScoreCols])

combine.df <- function(x, y, fillVal) {
  rows.x <- nrow(x)
  rows.y <- nrow(y)
  if (rows.x > rows.y) {
    diff <- rows.x - rows.y
    df.na <- matrix(fillVal, diff, ncol(y))
    colnames(df.na) <- colnames(y)
    cbind(x, rbind(y, df.na))
  } else {
    diff <- rows.y - rows.x
    df.na <- matrix(fillVal, diff, ncol(x))
    colnames(df.na) <- colnames(x)
    cbind(rbind(x, df.na), y)
  }
}

rowIds <- currData[,1]
featureMeans <- rowMeans(currData[,zScoreCols])
featureMedians <- rowMedians(as.matrix(currData[,zScoreCols]))

for(analysisCol in analysisCols) {

  analyteData <- currData[,analysisCol]
  numRows <- nrow(currData)
  sampleCols <- analysisCol
  
  specificityRatioCurrData <- as.matrix(analyteData) %o% 1/abs(featureMeans[1:length(featureMeans)])
  specificityRatioCurrData <- matrix(specificityRatioCurrData, ncol = length(sampleCols), nrow = numRows)
  colnames(specificityRatioCurrData) <- sampleCols
  rownames(specificityRatioCurrData) <- rowIds
  
  specificityMedianRatioCurrData <- as.matrix(analyteData) %o% 1/abs(featureMedians[1:length(featureMedians)])
  specificityMedianRatioCurrData <- matrix(specificityMedianRatioCurrData, ncol = length(sampleCols), nrow = numRows)
  colnames(specificityMedianRatioCurrData) <- sampleCols
  rownames(specificityMedianRatioCurrData) <- rowIds
  
  ZScore_Mean <- rowMeans(as.matrix(analyteData))
  ZScore_Median <- rowMedians(as.matrix(analyteData))
  Specificity_Mean <- rowMeans(cbind(specificityRatioCurrData))
  Specificity_Median <- rowMedians(cbind(specificityMedianRatioCurrData))
  
  attr(Specificity_Median,"names") <- rowIds
  
  if(Calculate_Nonspecific) {
    
    Z_Specificity_Idxs <- intersect(which(ZScore_Mean >= Z_Cutoff, arr.ind = TRUE), which(Specificity_Mean < Specificity_Cutoff, arr.ind = TRUE))
    Z_Specificity_Median_Idxs <- intersect(which(ZScore_Median >= Z_Cutoff, arr.ind = TRUE), which(Specificity_Median < Specificity_Cutoff, arr.ind = TRUE))
    
  } else {
    
    Z_Specificity_Idxs <- intersect(which(ZScore_Mean >= Z_Cutoff, arr.ind = TRUE), which(Specificity_Mean >= Specificity_Cutoff, arr.ind = TRUE))
    Z_Specificity_Median_Idxs <- intersect(which(ZScore_Median >= Z_Cutoff, arr.ind = TRUE), which(Specificity_Median >= Specificity_Cutoff, arr.ind = TRUE))
    
  }
  
  Z_Specificity_Hits <- cbind(ZScore_Mean[Z_Specificity_Idxs],Specificity_Mean[Z_Specificity_Idxs])
  Z_Specificity_Median_Hits <- cbind(ZScore_Median[Z_Specificity_Median_Idxs],Specificity_Median[Z_Specificity_Median_Idxs])
  
  #Z_Specificity <- combine.df(Z_Specificity_Hits, Z_Specificity_Median_Hits, "NA")
  Z_Specificity <- combine.df(cbind(rownames(Z_Specificity_Hits),Z_Specificity_Hits), cbind(rownames(Z_Specificity_Median_Hits),Z_Specificity_Median_Hits), "NA")
  colnames(Z_Specificity) <- outputColNames
  
  write.csv(Z_Specificity,paste0(outputDir,sampleCondition,"-",analysisCol,".csv"), row.names = FALSE)
  
}
