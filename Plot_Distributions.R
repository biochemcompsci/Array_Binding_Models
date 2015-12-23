twelveMin_Intersect <- intersect(Twelve_Min_3X$X, Twelve_Min_10X$X)
Twelve_Min_3X_Intersect_Idx <- match(twelveMin_Intersect, Twelve_Min_3X$X)
Twelve_Min_10X_Intersect_Idx <- match(twelveMin_Intersect, Twelve_Min_10X$X)
twelveMin_5969_10_3X <- cbind(Twelve_Min_3X$X[Twelve_Min_3X_Intersect_Idx], Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X_Intersect_Idx], Twelve_Min_10X$X[Twelve_Min_10X_Intersect_Idx], Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X_Intersect_Idx])

twelveMin_5969_10_3X <- as.data.frame(cbind(Twelve_Min_3X$X[Twelve_Min_3X_Intersect_Idx], as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X_Intersect_Idx]), Twelve_Min_10X$X[Twelve_Min_10X_Intersect_Idx], as.numeric(Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X_Intersect_Idx])), stringsAsFactors = FALSE)
colnames(twelveMin_5969_10_3X) <- c("3X ID", "12 min 3X","10X ID", "12 min 10X")

randIndices <- sample(1:length(twelveMin_5969_10_3X$V4),ceiling(length(twelveMin_5969_10_3X$V4)/2),replace=FALSE)

twelveMin_Scatter <- ggplot(data=twelveMin_5969_10_3X[randIndices,], aes(x=twelveMin_5969_10_3X$V2[randIndices],y=twelveMin_5969_10_3X$V4[randIndices])) + geom_point(na.rm = TRUE, alpha = 1/10) + theme(axis.ticks=element_blank(), axis.text = element_blank())

y_axis_ticks <- seq(floor(min(as.numeric(twelveMin_5969_10_3X$`12 min 10X`[randIndices]),na.rm = TRUE)),ceiling(max(as.numeric(twelveMin_5969_10_3X$`12 min 10X`[randIndices]),na.rm = TRUE)), 1)
x_axis_ticks <- seq(floor(min(as.numeric(twelveMin_5969_10_3X$`12 min 3X`[randIndices]),na.rm = TRUE)),ceiling(max(as.numeric(twelveMin_5969_10_3X$`12 min 3X`[randIndices]),na.rm = TRUE)), 1)
  
twelveMin_Scatter <- ggplot(data=twelveMin_5969_10_3X[randIndices,], aes(x=twelveMin_5969_10_3X$`12 min 3X`[randIndices],y=twelveMin_5969_10_3X$`12 min 10X`[randIndices])) + 
  geom_point(na.rm = TRUE, alpha = 1/10) + 
  coord_cartesian(ylim=c(min(y_axis_ticks), max(y_axis_ticks)),
                  xlim=c(min(x_axis_ticks), max(x_axis_ticks))) +
  scale_y_continuous(breaks = y_axis_ticks) + 
  scale_x_continuous(breaks = x_axis_ticks)