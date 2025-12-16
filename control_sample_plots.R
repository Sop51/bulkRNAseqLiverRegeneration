library(ggplot2)
library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggrepel)
library(RColorBrewer)

# read in the comparison files between the control groups
twovsone <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/qlf_timepointControl1.5vsControl0.csv")
threevstwo <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/qlf_timepointControl4.5vsControl1.5.csv")

# filter only for significant genes
twovsonesig <- twovsone[abs(twovsone$logFC) > 1 & twovsone$FDR < 0.05,]
threevstwosig <- threevstwo[abs(threevstwo$logFC) > 1 & threevstwo$FDR < 0.05,]

############################## DE GO and KEGG ###############################
# prepare input for twovsone
entrez_ids_twovsone <- bitr(
  twovsonesig$X,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Dr.eg.db
)$ENTREZID

# go two vs one
twovsone_go <- enrichGO(
  gene          = entrez_ids_twovsone,
  OrgDb         = org.Dr.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",     # biological process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
# plot
p_twovsone_go <- dotplot(twovsone_go)
ggsave(filename = "/Users/sm2949/Desktop/bulkRNAdata/plots/GO_upregulated_and_downregulated_Control1.5dpavsControl0dpa.png", plot = p_twovsone_go, width = 8, height = 6)

# kegg two vs one
twovsone_kegg <- enrichKEGG(
  gene          = entrez_ids_twovsone,
  organism      = "dre",
  pvalueCutoff  = 0.05
)
# plot
p_twovsone_kegg <- dotplot(twovsone_kegg)
ggsave(filename = "/Users/sm2949/Desktop/bulkRNAdata/plots/KEGG_upregulated_and_downregulated_Control1.5dpavsControl0dpa.png", plot = p_twovsone_kegg, width = 8, height = 6)

# prepare input for threevstwo
entrez_ids_threevstwo <- bitr(
  threevstwosig$X,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Dr.eg.db
)$ENTREZID

# go threevstwo
threevstwo_go <- enrichGO(
  gene          = entrez_ids_threevstwo,
  OrgDb         = org.Dr.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",     # biological process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
# plot
p_threevstwo_go <- dotplot(threevstwo_go)
ggsave(filename = "/Users/sm2949/Desktop/bulkRNAdata/plots/GO_upregulated_and_downregulated_Control4.5dpavsControl1.5dpa.png", plot = p_threevstwo_go, width = 8, height = 6)


# kegg threevstwo
threevstwo_kegg <- enrichKEGG(
  gene          = entrez_ids_threevstwo,
  organism      = "dre",
  pvalueCutoff  = 0.05
)
# plot
p_threevstwo_kegg <- dotplot(threevstwo_kegg)
ggsave(filename = "/Users/sm2949/Desktop/bulkRNAdata/plots/KEGG_upregulated_and_downregulated_Control4.5dpavsControl1.5dpa.png", plot = p_threevstwo_kegg, width = 8, height = 6)

############################## pca plot! ###############################
# read in the normalized counts matrix
cpms <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/normalized_counts.csv", row.names = 1)
colnames(cpms) <- gsub("^X", "", colnames(cpms)) # fix format
cpms_control <- cpms[, endsWith(colnames(cpms), "C")] # only control

# read in the metadata
metaData <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/metaData.csv', header = TRUE, sep = ",")
# remove the samples that are outliers
metaData <- metaData[!(metaData$sample %in% c("2_16_C", "2_16_T")), ]
metadata_control <-metaData[metaData$treatment == "C",] #only control

# run prcomp
pca <- prcomp(t(cpms_control), scale. = TRUE)

# df for plotting
df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = rownames(pca$x),
  timepoint = metadata_control$time
)

# calculate % variance explained
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# create a properly ordered factor with more informative labels
df$timepoint <- factor(metadata_control$time,
                       levels = c(1, 2, 3),            
                       labels = c("0dpa", "1.5dpa", "4.5dpa")) 

# plot
ggplot(df, aes(x = PC1, y = PC2, color = timepoint)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(type = "norm", geom = "polygon", alpha = 0.15, fill = "gray50", color = NA) +  # confidence ellipses
  stat_ellipse(type = "norm", linewidth = 1.2) + 
  scale_color_brewer(palette = "Dark2") +         
  labs(
    title = "PCA of Control Samples",
    subtitle = paste0("Timepoint explains the majority of variation (PC1: ", percentVar[1], "%, PC2: ", percentVar[2], "%)"),
    color = "Timepoint",
    x = paste0("PC1 (", percentVar[1], "% of variance)"),
    y = paste0("PC2 (", percentVar[2], "% of variance)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

