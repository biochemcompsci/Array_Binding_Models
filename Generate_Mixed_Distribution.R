generateBimodalDistribution <- function(distType, sampleSize, dist1Mean, dist2Mean, dist1SD, dist2SD, dist1Weight, minVal, maxVal) {
  
  dist2Weight <- 1 - dist1Weight
  
  dist1SampleSize <- ceiling(sampleSize*dist1Weight)
  dist2SampleSize <- sampleSize - dist1SampleSize
  
  mixedDist1Vals <- rep(NA, dist1SampleSize)
  mixedDist2Vals <- rep(NA, dist2SampleSize)
  
  if(distType == "LogNorm") {
    
    minVal <- log(minVal)
    maxVal <- log(maxVal)
    
  }
  
  if(dist1Weight > 0) {
    
    for(i in 1:dist1SampleSize) {
      
      if(distType == "LogNorm") {
        
        currDistVal <- log(rlnorm(1, meanlog = dist1Mean, sdlog = dist1SD))
        
      } else {
        
        currDistVal <- rnorm(1, dist1Mean, dist1SD)
        
      }
      
      while(is.na(currDistVal) || currDistVal > maxVal || currDistVal < minVal) {
        
        if(distType == "LogNorm") {
          
          currDistVal <- log(rlnorm(1, meanlog = dist1Mean, sdlog = dist1SD))
          
        } else {
          
          currDistVal <- rnorm(1, dist1Mean, dist1SD)
          
        }
        
        #cat("minVal: ", minVal, "maxVal: ", maxVal, "currDist1Val: ", currDistVal, "\n")
      }
      
      mixedDist1Vals[i] <- currDistVal
      
    }
    
  }
  
  if(dist2Weight > 0) {
    
    for(j in 1:dist2SampleSize) {
      
      if(distType == "LogNorm") {
        
        currDistVal <- log(rlnorm(1, meanlog = dist2Mean, sdlog = dist2SD))
        
      } else {
        
        currDistVal <- rnorm(1, dist2Mean, dist2SD)
        
      }
      
      while(is.na(currDistVal) || currDistVal > maxVal || currDistVal < minVal) {
        
        if(distType == "LogNorm") {
          
          currDistVal <- log(rlnorm(1, meanlog = dist2Mean, sdlog = dist2SD))
          
        } else {
          
          currDistVal <- rnorm(1, dist2Mean, dist2SD)
          
        }
        
        #cat("minVal: ", minVal, "maxVal: ", maxVal, "currDist2Val: ", currDistVal, "\n")
        
      }
      
      mixedDist2Vals[j] <- currDistVal
      
    }
    
  }
  
  mixedDistVals <- append(mixedDist1Vals, mixedDist2Vals)
  
  if(distType == "LogNorm") {
    
    exp(mixedDistVals)
    
  } else {
  
    mixedDistVals
    
  }
  
} 