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

currData1 <- Sixty_Min_3X
currData2 <- Sixty_Min_10X

excludeSeqs <- c('RHSVV')
#excludeSeqs <- c()

currXLabel <- " AB1 60 min Incubation, 3X & 10X Cycle Wash w/o Epitope >100kDa"
plotDir <- "/Volumes/Macintosh HD/Users/mattg/Dropbox/HealthTell_MPG_PS/Research/Modeling/Experimental/Chagas/04DEC2015/GPR_Data/Red_Channel_AB1/60_min_incubation/dist_plots/no_epitope/"

includeIdxs1 <- c(rep(1:nrow(currData1)))
includeIdxs2 <- c(rep(1:nrow(currData2)))

for(currExcludeSeq in excludeSeqs) {
  
  includeIdxs1 <- includeIdxs1[which(!grepl(currExcludeSeq, currData1[includeIdxs1,1]), arr.ind = TRUE)]
  includeIdxs2 <- includeIdxs2[which(!grepl(currExcludeSeq, currData2[includeIdxs2,1]), arr.ind = TRUE)]
  
}

currData1 <- currData1[includeIdxs1,]
currData2 <- currData2[includeIdxs2,]

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
  
  GT1_0_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 0, arr.ind = TRUE))])
  GT1_1_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 1, arr.ind = TRUE))])
  GT1_2_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 2, arr.ind = TRUE))])
  GT1_4_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 4, arr.ind = TRUE))])
  GT1_8_Specific <- as.data.frame(as.numeric(currData1[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData1[,as.character(currCol)] >= 8, arr.ind = TRUE))])
  GT2_0_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 0, arr.ind = TRUE))])
  GT2_1_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 1, arr.ind = TRUE))])
  GT2_2_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 2, arr.ind = TRUE))])
  GT2_4_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 4, arr.ind = TRUE))])
  GT2_8_Specific <- as.data.frame(as.numeric(currData2[,as.character(currCol)])[intersect(which(specificityRatioCurrData1[,as.character(currCol)] >= 1.5, arr.ind = TRUE),which(currData2[,as.character(currCol)] >= 8, arr.ind = TRUE))])
  
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
  
  sampleGroupNames <- colnames(GT_0)
  
  GT_0_Long <- melt(GT_0, na.rm = TRUE)
  GT_1_Long <- melt(GT_1, na.rm = TRUE)
  GT_2_Long <- melt(GT_2, na.rm = TRUE)
  GT_4_Long <- melt(GT_4, na.rm = TRUE)
  GT_8_Long <- melt(GT_8, na.rm = TRUE)
  
  GT_0_Long <- GT_0_Long[GT_0_Long$value > 0,]
  GT_1_Long <- GT_1_Long[GT_1_Long$value > 0,]
  GT_2_Long <- GT_2_Long[GT_2_Long$value > 0,]
  GT_4_Long <- GT_4_Long[GT_4_Long$value > 0,]
  GT_8_Long <- GT_8_Long[GT_8_Long$value > 0,]
  
  rm(GT_0,GT_1,GT_2,GT_4,GT_8)
  
  currPlot <- ggplot(data=GT_0_Long, aes(GT_0_Long[[2]], xmin=0, xmax=75, fill=variable, linetype = variable)) +
    scale_linetype_manual(breaks = sampleGroupNames, values = c(7,5,4,3,2)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.15) + xlab(paste0(currCol, currXLabel, " Z > 0"))
  ggsave(paste0(plotDir,currCol, "_GT0_3X_10X_Specificity_Small.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
  
  currPlot <- ggplot(data=GT_1_Long, aes(GT_1_Long[[2]], xmin=0, xmax=75, fill=variable, linetype = variable)) +
    scale_linetype_manual(breaks = sampleGroupNames, values = c(7,5,4,3,2)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 1"))
  ggsave(paste0(plotDir,currCol, "_GT1_3X_10X_Specificity_Small.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
  
  currPlot <- ggplot(data=GT_2_Long, aes(GT_2_Long[[2]], xmin=0, xmax=75, fill=variable, linetype = variable)) +
    scale_linetype_manual(breaks = sampleGroupNames, values = c(7,5,4,3,2)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 2"))
  ggsave(paste0(plotDir,currCol, "_GT2_3X_10X_Specificity_Small.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
  
  currPlot <- ggplot(data=GT_4_Long, aes(GT_4_Long[[2]], xmin=0, xmax=75, fill=variable, linetype = variable)) +
    scale_linetype_manual(breaks = sampleGroupNames, values = c(7,5,4,3,2)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 4"))
  ggsave(paste0(plotDir,currCol, "_GT4_3X_10X_Specificity_Small.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
  
  currPlot <- ggplot(data=GT_8_Long, aes(GT_8_Long[[2]], xmin=0, xmax=75, fill=variable, linetype = variable)) +
    scale_linetype_manual(breaks = sampleGroupNames, values = c(7,5,4,3,2)) +
    theme(axis.title = element_text(size = 5)) +
    geom_density(alpha=0.2) + xlab(paste0(currCol, currXLabel, " Z > 8"))
  ggsave(paste0(plotDir,currCol, "_GT8_3X_10X_Specificity_Small.png"), plot = currPlot, height = 6, width = 9, units = "in", dpi=300)
  
}