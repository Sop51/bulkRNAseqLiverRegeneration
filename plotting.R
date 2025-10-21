library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(biomaRt)
library(org.Dr.eg.db)
library(clusterProfiler)

# read in the normalized counts
cpm <- read.csv('/Users/sophiemarcotte/Desktop/bulkRNAdata/normalized_counts.csv', row.names=1, check.names = FALSE)

# read in all edgeR result tables
timepoint1.5vs0 <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint1.5vs0.csv", row.names = 1)
timepoint4.5vs0 <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint4.5vs0.csv", row.names = 1)
timepoint4.5vs1.5 <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint4.5vs1.5.csv", row.names = 1)

timepoint0vsAll <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint0vsAll.csv", row.names = 1)
timepoint1.5vsAll <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint1.5vsAll.csv", row.names = 1)
timepoint4.5vsAll <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_timepoint4.5vsAll.csv", row.names = 1)

AllTvsAllC <- read.csv("/Users/sophiemarcotte/Desktop/bulkRNAdata/qlf_AllTvsAllC.csv", row.names = 1)

########## PLOTTING A SPECIFIC GENE ACROSS TIME ############
# plot a specific gene
gene <- "ENSDARG00000019949"
gene_symbol <- "serpinh1b"
gene_counts <- cpm[gene, ]

# convert the data to long format for plotting
data_long <- pivot_longer(gene_counts, 
                          cols = everything(), 
                          names_to = "sample", 
                          values_to = "gene_count")

# create a new column for timepoint based on sample names
data_long$timepoint <- gsub("^(\\d)_.*", "\\1", data_long$sample)
data_long$group <- ifelse(grepl("_C$", data_long$sample), "Control", "Treatment")

# plot the genes
ggplot(data_long, aes(x = timepoint, y = gene_count, fill = group)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = group), size = 3, shape = 21, position = position_dodge(width = 0.8)) +
  labs(
    title = paste(gene_symbol),
    x = "Timepoint",
    y = "CPM Normalized Expression"
  ) +
  scale_fill_manual(values = c("Control" = "#69b3a2", "Treatment" = "#e6b800")) +  # custom colors for groups
  scale_color_manual(values = c("Control" = "#69b3a2", "Treatment" = "#e6b800")) +  # match point color to fill color
  scale_x_discrete(labels = c("1" = "0dpa", "2" = "1.5dpa", "3" = "4.5dpa")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),        # increase text size
    axis.title = element_text(size = 16),  # larger axis titles
    axis.text = element_text(size = 12),   # adjust axis text size
    plot.title = element_text(size = 18, hjust = 0.5),  # center title
    legend.title = element_blank()         
  ) +
  stat_pwc(
    method = 't_test',
    label = "p.adj.signif",
    group.by = 'timepoint'
  )

################# VOLCANO PLOTS ################
results <- timepoint4.5vs1.5
title_p <- "DE Genes: 4.5 dpa vs 1.5 dpa"

# connect to zebrafish database
ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# get gene symbols for Ensembl IDs
gene_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values = rownames(results),
  mart = ensembl
)

# add gene symbols to results
results$ensembl_id <- rownames(results)
results <- merge(results, gene_annotations, 
                 by.x = "ensembl_id", by.y = "ensembl_gene_id", 
                 all.x = TRUE)

# use gene symbol if available, otherwise use Ensembl ID
results$gene_label <- ifelse(is.na(results$external_gene_name) | results$external_gene_name == "", 
                             results$ensembl_id, 
                             results$external_gene_name)

rownames(results) <- results$ensembl_id

# add a column for -log10(FDR)
results$neg_log10_FDR <- -log10(results$FDR)

# qdd a column to classify genes as UP, DOWN, or NS (not significant)
results$DE_status <- "NS"
results$DE_status[results$logFC > 2.5 & results$FDR < 0.001] <- "UP"
results$DE_status[results$logFC < -2.5 & results$FDR < 0.001] <- "DOWN"

# get top 10 upregulated and top 10 downregulated genes
top_up <- head(results[results$DE_status == "UP", ][order(results[results$DE_status == "UP", ]$FDR), ], 10)
top_down <- head(results[results$DE_status == "DOWN", ][order(results[results$DE_status == "DOWN", ]$FDR), ], 10)
top_genes <- rbind(top_up, top_down)

# create volcano plot with gene symbol labels
volcano_plot_labeled <- ggplot(results, aes(x = logFC, y = neg_log10_FDR, color = DE_status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("UP" = "#F9C784", "DOWN" = "#70B0EC", "NS" = "grey"),
                     labels = c("UP" = "Upregulated", "DOWN" = "Downregulated", "NS" = "Not Significant")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = gene_label),
                  size = 3, max.overlaps = 20, 
                  box.padding = 0.5, point.padding = 0.3,
                  color = "black", show.legend = FALSE) +
  labs(title = title_p,
       x = "Log2 Fold Change",
       y = "-Log10(FDR)",
       color = "Differential Expression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

ggsave("/Users/sophiemarcotte/Desktop/bulkRNAdata/plots/volcano_plot_timepoint4.5vs1.5.png", plot = volcano_plot_labeled, width = 9, height = 5, dpi = 300, bg = "white")

############### KEGG and GO PLOTS #################
results <- timepoint4.5vs1.5
title_p <- "4.5 vs 1.5 dpa"

# add a column for -log10(FDR)
results$neg_log10_FDR <- -log10(results$FDR)

# add a column to classify genes as UP, DOWN, or NS (not significant)
results$DE_status <- "NS"
results$DE_status[results$logFC > 1 & results$FDR < 0.05] <- "UP"
results$DE_status[results$logFC < -1 & results$FDR < 0.05] <- "DOWN"

# get Ensembl IDs for upregulated and downregulated genes
up_ensembl <- rownames(results)[results$DE_status == "UP"]
down_ensembl <- rownames(results)[results$DE_status == "DOWN"]

# convert Ensembl IDs to Entrez IDs for zebrafish
up_entrez <- bitr(up_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", 
                  OrgDb = org.Dr.eg.db)$ENTREZID
down_entrez <- bitr(down_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", 
                    OrgDb = org.Dr.eg.db)$ENTREZID

# GO enrichment for upregulated genes
go_up <- enrichGO(gene = up_entrez,
                  OrgDb = org.Dr.eg.db,
                  ont = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# GO enrichment for downregulated genes
go_down <- enrichGO(gene = down_entrez,
                    OrgDb = org.Dr.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# KEGG enrichment for upregulated genes
kegg_up <- enrichKEGG(gene = up_entrez,
                      organism = "dre",  
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# KEGG enrichment for downregulated genes
kegg_down <- enrichKEGG(gene = down_entrez,
                        organism = "dre",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

# GO Dotplot - Upregulated
if (!is.null(go_up) && nrow(go_up@result) > 0) {
  go_up_plot <- dotplot(go_up, showCategory = 10, title = paste("GO Enrichment - Upregulated Genes", title_p)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(go_up_plot)
  ggsave("/Users/sophiemarcotte/Desktop/bulkRNAdata/plots/GO_upregulated_timepoint4.5vs1.5.png", plot = go_up_plot, width = 10, height = 7, dpi = 300, bg = "white")
}

# GO Dotplot - Downregulated
if (!is.null(go_down) && nrow(go_down@result) > 0) {
  go_down_plot <- dotplot(go_down, showCategory = 10, title = paste("GO Enrichment - Downregulated Genes", title_p)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(go_down_plot)
  ggsave("/Users/sophiemarcotte/Desktop/bulkRNAdata/plots/GO_downregulated_timepoint4.5vs1.5.png", plot = go_down_plot, width = 10, height = 7, dpi = 300, bg = "white")
}

# KEGG Dotplot - Upregulated
if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
  kegg_up_plot <- dotplot(kegg_up, showCategory = 10, title = paste("KEGG Pathways - Upregulated Genes", title_p)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(kegg_up_plot)
  ggsave("/Users/sophiemarcotte/Desktop/bulkRNAdata/plots/KEGG_upregulated_timepoint4.5vs1.5.png", plot = kegg_up_plot, width = 10, height = 7, dpi = 300, bg = "white")
}

# KEGG Dotplot - Downregulated
if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
  kegg_down_plot <- dotplot(kegg_down, showCategory = 10, title = paste("KEGG Pathways - Downregulated Genes", title_p)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(kegg_down_plot)
  ggsave("/Users/sophiemarcotte/Desktop/bulkRNAdata/plots/KEGG_downregulated_timepoint4.5vs1.5.png", plot = kegg_down_plot, width = 10, height = 7, dpi = 300, bg = "white")
}
