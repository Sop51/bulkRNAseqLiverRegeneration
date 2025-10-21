library(WGCNA)
library(tidyverse)
library(DESeq2)
library(CorLevelPlot)
library(biomaRt)

# organize the metadata ----
# read in the metadata
metaData <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/metaData.csv', header = TRUE, sep = ",")
# edit the name of the samples
metaData$sample <- gsub("20240430_(\\d+_\\d+_[CT])_.*", "\\1", metaData$sample)
# remove the samples that are outliers
metaData <- metaData[!(metaData$sample %in% c("2_16_C", "2_16_T")), ]
# convert treatment and timepoint to factors
metaData$treatment <- factor(metaData$treatment)
metaData$time <- factor(metaData$time)
rownames(metaData) <- metaData$sample


# organize the counts data ----
# read in the counts data
countData <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/countMatrix.csv', header = TRUE, sep = ",")
colnames(countData) <- gsub("^X", "", colnames(countData))
# Set the first column as row names
rownames(countData) <- countData[[1]] 
# Remove the first column
countData <- countData[, -1]

# check for outliers
gsg <- goodSamplesGenes(t(countData))
summary(gsg)
gsg$allOK

# create a deseq2 dataset
all(rownames(metaData) %in% colnames(countData))
all(rownames(metaData) == colnames(countData))

# create dds
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ 1)

# filter
dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]

# vst
dds_norm <- vst(dds75)

# reverse format
norm.counts <- assay(dds_norm) %>%
  t()

# network construction
power <- c(c(1:10), seq(from = 12, to = 50, by =2))

sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = 'signed',
                  verbose = 5)

# create ajacency matrix
norm.counts[] <- sapply(norm.counts, as.numeric)
soft_power <- 16
temp_cor <- cor

cor <- WGCNA::cor

bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 5000,
                 TOMType = 'signed',
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

cor <- temp_cor

# module eigengenes
module_eigengenes <- bwnet$MEs

# genes per module
table(bwnet$colors)

# relate modules to traits
metaData <- metaData %>% 
            mutate(treat_bin = ifelse(grepl('T', treatment), 1, 0)) %>% 
            mutate(time0_bin =ifelse(grepl(1, time), 1, 0)) %>% 
            mutate(time1.5_bin =ifelse(grepl(2, time), 1, 0)) %>%
            mutate(time4.5_bin =ifelse(grepl(3, time), 1, 0)) %>%
            select(5:8)

# define the numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, metaData, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module trait association as a heatmap
heatmap.data <- merge(module_eigengenes, metaData, by = 'row.names')
heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:21],
             y = names(heatmap.data)[1:17])


# genes in a module
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>%
  filter(bwnet$colors == 'turquoise') %>%
  rownames()

# identifying driver genes
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pval <- corPvalueStudent(module.membership.measure, nSamples)

# top 10 sig driver genes
module_name = 'MEbrown' 
module_pvals <- module.membership.measure.pval[module_name, ]
sorted_genes <- sort(module_pvals)
top_50_genes <- names(sorted_genes)[1:100]

ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
# query biomart for gene symbols
genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = top_50_genes,
  mart = ensembl
)
print(genes$external_gene_name)


# top sig genes associated with a metadata col
gene.signf.corr <- cor(norm.counts, metaData$treat_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>%
  arrange(V1) %>%
  head(25)
