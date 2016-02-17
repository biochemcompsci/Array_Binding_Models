distributionsDistanceMetric <- function(distribution1Vector, distribution2Vector) {
  
  dist(rbind(distribution1Vector, distribution2Vector))
  
}