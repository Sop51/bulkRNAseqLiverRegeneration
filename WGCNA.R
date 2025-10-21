library(WGCNA)
library(tidyverse)
library(DESeq2)

# organize the metadata ----
# read in the metadata
metaData <- read.csv('/Users/sophiemarcotte/Desktop/bulkRNAdata/metaData.csv', header = TRUE, sep = ",")
# edit the name of the samples
metaData$sample <- gsub("20240430_(\\d+_\\d+_[CT])_.*", "\\1", metaData$sample)
# remove the samples that are outliers
metaData <- metaData[!(metaData$sample %in% c("2_16_C", "2_16_T")), ]
# convert treatment and timepoint to factors
metaData$treatment <- factor(metaData$treatment)
metaData$time <- factor(metaData$time)


# organize the counts data ----
# read in the counts data
cpm <- read.csv('/Users/sophiemarcotte/Desktop/bulkRNAdata/normalized_counts.csv', row.names=1, check.names = FALSE)
dataExpr <- t(cpm)

# choose a softhold threshold power
powers = c(c(1:20), seq(from = 22, to = 30, by = 2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose=5)

par(mfrow = c(1,2));
cex1 = 0.9;

# plot to visualize scale independence and mean connectivity
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# choose power from there ^
softPower = 3
temp_cor <- cor       
cor <- WGCNA::cor 

# run
netwk <- blockwiseModules(dataExpr,
                        power = softPower, 
                        networkType = "unsigned",
                        deepSplit = 2,
                        pamRespectsDendro = F,
                        minModuleSize = 30,
                        maxBlockSize = 4000,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        verbose = 3)

cor <- temp_cor

MEs <- netwk$MEs
moduleTraitCor <- cor(MEs, sampleTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(dataExpr))

# Heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sampleTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module-trait relationships")
