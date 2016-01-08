library(ggplot2)
library(reshape2)
library(dplyr)

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

currData1 <- OneTwenty_Min_3X
currData2 <- OneTwenty_Min_10X

currXLabel <- " AB1 120 min Incubation, 3X & 10X Cycle Wash"
plotDir <- "/Volumes/Macintosh HD/Users/mattg/Dropbox/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/24NOV2015/GPR_Data/Red_Channel_AB1/120_min_incubation/dist_plots/"

sampleCols <- colnames(currData1)[2:ncol(currData1)]

# Calculate Matrix of Signal vs. Mean for each Sample and Feature Mean (average across samples)
rowMeansCurrData1 <- abs(rowMeans(currData1[2:ncol(currData1)]))
rowMeansCurrData2 <- abs(rowMeans(currData2[2:ncol(currData2)]))

specificityRatioCurrDatatmp <- as.matrix(currData1[1:nrow(currData1),2:ncol(currData1)])
specificityRatioCurrData1 <- specificityRatioCurrDatatmp %o% 1/rowMeansCurrData1[1:length(rowMeansCurrData1)]
specificityRatioCurrData1 <- matrix(specificityRatioCurrData1, ncol = ncol(specificityRatioCurrDatatmp), nrow = nrow(specificityRatioCurrDatatmp))
colnames(specificityRatioCurrData1) <- sampleCols
rownames(specificityRatioCurrData1) <- currData1[1:nrow(currData1),1]

specificityRatioCurrDatatmp <- as.matrix(currData2[1:nrow(currData2),2:ncol(currData2)])
specificityRatioCurrData2 <- specificityRatioCurrDatatmp %o% 1/rowMeansCurrData2[1:length(rowMeansCurrData2)]
specificityRatioCurrData2 <- matrix(specificityRatioCurrData2, ncol = ncol(specificityRatioCurrDatatmp), nrow = nrow(specificityRatioCurrDatatmp))
colnames(specificityRatioCurrData2) <- sampleCols
rownames(specificityRatioCurrData2) <- currData2[1:nrow(currData2),1]

rm(specificityRatioCurrDatatmp)

for(currCol in sampleCols) {
  
  GT1_0 <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[currData1[,as.character(currCol)] >= 0])
  GT1_1 <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[currData1[,as.character(currCol)] >= 1])
  GT1_2 <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[currData1[,as.character(currCol)] >= 2])
  GT1_4 <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[currData1[,as.character(currCol)] >= 4])
  GT1_8 <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[currData1[,as.character(currCol)] >= 8])
  GT2_0 <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[currData2[,as.character(currCol)] >= 0])
  GT2_1 <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[currData2[,as.character(currCol)] >= 1])
  GT2_2 <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[currData2[,as.character(currCol)] >= 2])
  GT2_4 <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[currData2[,as.character(currCol)] >= 4])
  GT2_8 <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[currData2[,as.character(currCol)] >= 8])
  
  GT_0 <- combine.df(GT1_0,GT2_0, -1)
  GT_1 <- combine.df(GT1_1,GT2_1, -1)
  GT_2 <- combine.df(GT1_2,GT2_2, -1)
  GT_4 <- combine.df(GT1_4,GT2_4, -1)
  GT_8 <- combine.df(GT1_8,GT2_8, -1)
  
  rm(GT1_0,GT2_0,GT1_1,GT2_1,GT1_2,GT2_2,GT1_4,GT2_4,GT1_8,GT2_8)
  
  GT1_0_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 0, arr.ind = TRUE))])
  GT1_1_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 1, arr.ind = TRUE))])
  GT1_2_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 2, arr.ind = TRUE))])
  GT1_4_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 4, arr.ind = TRUE))])
  GT1_8_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 8, arr.ind = TRUE))])
  GT2_0_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 0, arr.ind = TRUE))])
  GT2_1_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 1, arr.ind = TRUE))])
  GT2_2_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 2, arr.ind = TRUE))])
  GT2_4_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 4, arr.ind = TRUE))])
  GT2_8_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 2, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 8, arr.ind = TRUE))])
  
  GT_0_Specific <- combine.df(GT1_0_Specific,GT2_0_Specific, -1)
  GT_1_Specific <- combine.df(GT1_1_Specific,GT2_1_Specific, -1)
  GT_2_Specific <- combine.df(GT1_2_Specific,GT2_2_Specific, -1)
  GT_4_Specific <- combine.df(GT1_4_Specific,GT2_4_Specific, -1)
  GT_8_Specific <- combine.df(GT1_8_Specific,GT2_8_Specific, -1)
  
  rm(GT1_0_Specific,GT2_0_Specific,GT1_1_Specific,GT2_1_Specific,GT1_2_Specific,GT2_2_Specific,GT1_4_Specific,GT2_4_Specific,GT1_8_Specific,GT2_8_Specific)
  
  GT_0 <- combine.df(GT_0,GT_0_Specific, -1)
  GT_1 <- combine.df(GT_1,GT_1_Specific, -1)
  GT_2 <- combine.df(GT_2,GT_2_Specific, -1)
  GT_4 <- combine.df(GT_4,GT_4_Specific, -1)
  GT_8 <- combine.df(GT_8,GT_8_Specific, -1)
  
  rm(GT_0_Specific,GT_1_Specific,GT_2_Specific,GT_4_Specific,GT_8_Specific)
  
  colnames(GT_0) <- c(paste0(currCol,"_3X"),paste0(currCol,"_10X"),paste0(currCol,"_3X_Specific"),paste0(currCol,"_10X_Specific"))
  colnames(GT_1) <- c(paste0(currCol,"_3X"),paste0(currCol,"_10X"),paste0(currCol,"_3X_Specific"),paste0(currCol,"_10X_Specific"))
  colnames(GT_2) <- c(paste0(currCol,"_3X"),paste0(currCol,"_10X"),paste0(currCol,"_3X_Specific"),paste0(currCol,"_10X_Specific"))
  colnames(GT_4) <- c(paste0(currCol,"_3X"),paste0(currCol,"_10X"),paste0(currCol,"_3X_Specific"),paste0(currCol,"_10X_Specific"))
  colnames(GT_8) <- c(paste0(currCol,"_3X"),paste0(currCol,"_10X"),paste0(currCol,"_3X_Specific"),paste0(currCol,"_10X_Specific"))
  
  GT_0_Long <- melt(GT_0, na.rm = TRUE)
  GT_1_Long <- melt(GT_1, na.rm = TRUE)
  GT_2_Long <- melt(GT_2, na.rm = TRUE)
  GT_4_Long <- melt(GT_4, na.rm = TRUE)
  GT_8_Long <- melt(GT_8, na.rm = TRUE)
  
  rm(GT_0,GT_1,GT_2,GT_4,GT_8)
  
  currPlot <- ggplot(data=GT_0_Long, aes(GT_0_Long[[2]], xmin=0, xmax=75, fill=variable)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 0"))
  ggsave(paste0(plotDir,currCol, "_GT0_3X_10X_Specificity_Small.png"), plot = currPlot, height = 3.25, width = 6.5, units = "in", dpi=200)
  
  currPlot <- ggplot(data=GT_1_Long, aes(GT_1_Long[[2]], xmin=0, xmax=75, fill=variable)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 1"))
  ggsave(paste0(plotDir,currCol, "_GT1_3X_10X_Specificity_Small.png"), plot = currPlot, height = 3.25, width = 6.5, units = "in", dpi=200)
  
  currPlot <- ggplot(data=GT_2_Long, aes(GT_2_Long[[2]], xmin=0, xmax=75, fill=variable)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 2"))
  ggsave(paste0(plotDir,currCol, "_GT2_3X_10X_Specificity_Small.png"), plot = currPlot, height = 3.25, width = 6.5, units = "in", dpi=200)
  
  currPlot <- ggplot(data=GT_4_Long, aes(GT_4_Long[[2]], xmin=0, xmax=75, fill=variable)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 4"))
  ggsave(paste0(plotDir,currCol, "_GT4_3X_10X_Specificity_Small.png"), plot = currPlot, height = 3.25, width = 6.5, units = "in", dpi=200)
  
  currPlot <- ggplot(data=GT_8_Long, aes(GT_8_Long[[2]], xmin=0, xmax=75, fill=variable)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 8"))
  ggsave(paste0(plotDir,currCol, "_GT8_3X_10X_Specificity_Small.png"), plot = currPlot, height = 3.25, width = 6.5, units = "in", dpi=200)
  
}


twelveMin_3X_GT_0 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 0]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_0) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_0, aes(twelveMin_3X_GT_0$Chagas.G1.5969, xmin=0, fill="12 Min. Inc., 3X Wash"))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_density(alpha = 0.3) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

twelveMin_3X_GT_1 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 1]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_1) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_1, aes(twelveMin_3X_GT_1$Chagas.G1.5969, xmin=1, fill="12 Min. Inc., 3X Wash"))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_density(alpha = 0.3) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

twelveMin_3X_GT_2 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 2]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_2) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_2, aes(twelveMin_3X_GT_2$Chagas.G1.5969, xmin=2, fill="12 Min. Inc., 3X Wash"))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_density(alpha = 0.3) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

twelveMin_10X_GT_0 <- as.data.frame(as.numeric(Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X$Chagas.G1.5969 >= 0]), stringsAsFactors = FALSE)
colnames(twelveMin_10X_GT_0) <- c("Chagas.G1.5969")
twelveMin_10X_G1_5969 <- ggplot(data=twelveMin_10X_GT_0, aes(twelveMin_10X_GT_0$Chagas.G1.5969, xmin=0, fill="12 Min. Inc., 10X Wash"))
twelveMin_10X_G1_5969 <- twelveMin_10X_G1_5969 + geom_density(alpha = 0.3) + xlab("Chagas G1 #5969: 12 min Incubation, 10X Cycle Wash")
twelveMin_10X_G1_5969

twelveMin_10X_GT_1 <- as.data.frame(as.numeric(Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X$Chagas.G1.5969 >= 1]), stringsAsFactors = FALSE)
colnames(twelveMin_10X_GT_1) <- c("Chagas.G1.5969")
twelveMin_10X_G1_5969 <- ggplot(data=twelveMin_10X_GT_1, aes(twelveMin_10X_GT_1$Chagas.G1.5969, xmin=1, fill="12 Min. Inc., 10X Wash"))
twelveMin_10X_G1_5969 <- twelveMin_10X_G1_5969 + geom_density(alpha = 0.3) + xlab("Chagas G1 #5969: 12 min Incubation, 10X Cycle Wash")
twelveMin_10X_G1_5969

twelveMin_10X_GT_2 <- as.data.frame(as.numeric(Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X$Chagas.G1.5969 >= 2]), stringsAsFactors = FALSE)
colnames(twelveMin_10X_GT_2) <- c("Chagas.G1.5969")
twelveMin_10X_G1_5969 <- ggplot(data=twelveMin_10X_GT_2, aes(twelveMin_10X_GT_2$Chagas.G1.5969, xmin=2, fill="12 Min. Inc., 10X Wash"))
twelveMin_10X_G1_5969 <- twelveMin_10X_G1_5969 + geom_density(alpha = 0.3) + xlab("Chagas G1 #5969: 12 min Incubation, 10X Cycle Wash")
twelveMin_10X_G1_5969


twelveMin_3X_GT_0 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 0]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_0) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_0, aes(twelveMin_3X_GT_0$Chagas.G1.5969, ..density.., xmin=0))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_histogram(binwidth = 0.25) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

twelveMin_3X_GT_1 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 1]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_1) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_1, aes(twelveMin_3X_GT_1$Chagas.G1.5969, ..density.., xmin=1))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_histogram(binwidth = 0.25) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

twelveMin_3X_GT_2 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 2]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_2) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_2, aes(twelveMin_3X_GT_2$Chagas.G1.5969, ..density.., xmin=2))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_histogram(binwidth = 0.25) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

twelveMin_3X_GT_4 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 4]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_4) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_4, aes(twelveMin_3X_GT_4$Chagas.G1.5969, ..density.., xmin=4))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_histogram(binwidth = 0.25) + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969



twelveMin_3X_GT_4 <- as.data.frame(as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X$Chagas.G1.5969 >= 4]), stringsAsFactors = FALSE)
colnames(twelveMin_3X_GT_4) <- c("Chagas.G1.5969")
twelveMin_3X_G1_5969 <- ggplot(data=twelveMin_3X_GT_4, aes(twelveMin_3X_GT_4$Chagas.G1.5969, xmin=4))
twelveMin_3X_G1_5969 <- twelveMin_3X_G1_5969 + geom_density() + xlab("Chagas G1 #5969: 12 min Incubation, 3X Cycle Wash")
twelveMin_3X_G1_5969

# twelveMin_Intersect <- intersect(Twelve_Min_3X$X, Twelve_Min_10X$X)
# Twelve_Min_3X_Intersect_Idx <- match(twelveMin_Intersect, Twelve_Min_3X$X)
# Twelve_Min_10X_Intersect_Idx <- match(twelveMin_Intersect, Twelve_Min_10X$X)
# twelveMin_5969_10_3X <- cbind(Twelve_Min_3X$X[Twelve_Min_3X_Intersect_Idx], Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X_Intersect_Idx], Twelve_Min_10X$X[Twelve_Min_10X_Intersect_Idx], Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X_Intersect_Idx])
# 
# twelveMin_5969_10_3X <- as.data.frame(cbind(Twelve_Min_3X$X[Twelve_Min_3X_Intersect_Idx], as.numeric(Twelve_Min_3X$Chagas.G1.5969[Twelve_Min_3X_Intersect_Idx]), Twelve_Min_10X$X[Twelve_Min_10X_Intersect_Idx], as.numeric(Twelve_Min_10X$Chagas.G1.5969[Twelve_Min_10X_Intersect_Idx])), stringsAsFactors = FALSE)
# colnames(twelveMin_5969_10_3X) <- c("3X ID", "12 min 3X","10X ID", "12 min 10X")
# 
# randIndices <- sample(1:length(twelveMin_5969_10_3X$V4),ceiling(length(twelveMin_5969_10_3X$V4)/2),replace=FALSE)
# 
# twelveMin_Scatter <- ggplot(data=twelveMin_5969_10_3X[randIndices,], aes(x=twelveMin_5969_10_3X$V2[randIndices],y=twelveMin_5969_10_3X$V4[randIndices])) + geom_point(na.rm = TRUE, alpha = 1/10) + theme(axis.ticks=element_blank(), axis.text = element_blank())
# 
# y_axis_ticks <- seq(floor(min(as.numeric(twelveMin_5969_10_3X$`12 min 10X`[randIndices]),na.rm = TRUE)),ceiling(max(as.numeric(twelveMin_5969_10_3X$`12 min 10X`[randIndices]),na.rm = TRUE)), 0.5)
# x_axis_ticks <- seq(floor(min(as.numeric(twelveMin_5969_10_3X$`12 min 3X`[randIndices]),na.rm = TRUE)),ceiling(max(as.numeric(twelveMin_5969_10_3X$`12 min 3X`[randIndices]),na.rm = TRUE)), 0.5)
#   
# twelveMin_Scatter <- ggplot(data=twelveMin_5969_10_3X[randIndices,], aes(x=twelveMin_5969_10_3X$`12 min 3X`[randIndices],y=twelveMin_5969_10_3X$`12 min 10X`[randIndices])) + 
#   geom_point(na.rm = TRUE, alpha = 1/10) + 
#   ylab(colnames(twelveMin_5969_10_3X)[4]) +
#   xlab(colnames(twelveMin_5969_10_3X)[2]) +
#   theme(axis.ticks=element_blank(), axis.text = element_blank()) +
#   scale_y_discrete(breaks = y_axis_ticks) + 
#   scale_x_discrete(breaks = x_axis_ticks)
# 
# chart.Correlation(as.matrix(cbind(as.numeric(twelveMin_5969_10_3X$`12 min 3X`[randIndices]),as.numeric(twelveMin_5969_10_3X$`12 min 10X`[randIndices]))),histogram = TRUE, method="spearman")
# 
# chart.Histogram(as.vector(as.numeric(twelveMin_5969_10_3X$`12 min 3X`)),methods = c("add.density", "add.rug"))