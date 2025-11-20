library(edgeR)
library(ggplot2)

############## INITIAL EDGER RUN ###############
# organize the metadata for edge R ----
# read in the metadata
metaData <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/metaData.csv', header = TRUE, sep = ",")
# remove the samples that are outliers
metaData <- metaData[!(metaData$sample %in% c("2_16_C", "2_16_T")), ]
# convert treatment and timepoint to factors
metaData$treatment <- factor(metaData$treatment)
metaData$time <- factor(metaData$time)
metaData$treat_time <- factor(metaData$treat_time)

# organize the counts data ----
# read in the counts data
countData <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/countMatrix.csv', header = TRUE, sep = ",")
colnames(countData) <- gsub("^X", "", colnames(countData))
# Set the first column as row names
rownames(countData) <- countData[[1]] 
# Remove the first column
countData <- countData[, -1]  

# running edge R as a time series ----
# create the DGElist object
y <- DGEList(counts=countData,group=metaData$treat_time)
# filter by expression
keep <- filterByExpr(y)
# filter by library size
y <- y[keep,,keep.lib.sizes=FALSE]
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
# create the model design
design <- model.matrix(~0+metaData$treat_time)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(metaData$treat_time)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

# compute the glm QLF test for different contrasts
qlf_treatvscontrol <- glmQLFTest(fit, contrast = c(-1, -1, -1, 1, 1, 1))
qlf_timepoint2vs1 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 1, 0))
qlf_timepoint3vs1 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 0, 1))
qlf_timepoint3vs2 <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, -1, 1))
qlf_timepoint1 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 1, 0, 0))
qlf_timepoint2 <- glmQLFTest(fit, contrast = c(0, -1, 0, 0, 1, 0))
qlf_timepoint3 <- glmQLFTest(fit, contrast = c(0, 0, -1, 0, 0, 1))
qlf_timepointControl2vsControl1 <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0))
qlf_timepointControl3vsControl2 <- glmQLFTest(fit, contrast = c(0, -1, 1, 0, 0, 0))

# export results from edgeR to CSV
write.csv(topTags(qlf_timepoint2vs1, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint1.5vs0.csv")
write.csv(topTags(qlf_timepoint3vs1, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint4.5vs0.csv")
write.csv(topTags(qlf_timepoint3vs2, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint4.5vs1.5.csv")
write.csv(topTags(qlf_timepoint1, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint0vsAll.csv")
write.csv(topTags(qlf_timepoint2, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint1.5vsAll.csv")
write.csv(topTags(qlf_timepoint3, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint4.5vsAll.csv")
write.csv(topTags(qlf_treatvscontrol, n=Inf)$table, file="/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_AllTvsAllC.csv")
write.csv(topTags(qlf_timepointControl2vsControl1, n=Inf)$table, file="/Users/sm2949/Desktop/bulkRNAdata/qlf_timepointControl1.5vsControl0.csv")
write.csv(topTags(qlf_timepointControl3vsControl2, n=Inf)$table, file="/Users/sm2949/Desktop/bulkRNAdata/qlf_timepointControl4.5vsControl1.5.csv")


# look at the results of all treated vs control
results <- as.data.frame(qlf_treatvscontrol)
resultsT3 <- as.data.frame(qlf_timepoint3)
resultsT2 <- as.data.frame(qlf_timepoint2)
resultsT1 <- as.data.frame(qlf_timepoint1)

############## NORMALIZE RAW DATA TO WORK WITH IN CPM ##############
cpms <- edgeR::cpm(y, log = TRUE)
cpms <- as.data.frame(cpms)
write.csv(cpms, file="/Users/sm2949/Desktop/bulkRNAdata/normalized_counts.csv")

############### PCA PLOT ##################
#run pca
pca <- prcomp(t(cpms), scale. = TRUE)
# make data frame for plotting
pca_df <- data.frame(pca$x[,1:2], group = y$samples$group)
# plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))
