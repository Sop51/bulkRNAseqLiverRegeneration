library("biomaRt")
library(ggplot2)
library(ggrepel)
library(dplyr)

# read in inflammation genes
inflammation_genes <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/GO_inflammatory_response.txt', col.names = 'genes')

# add ensembl gene names 
ensembl = useMart("ensembl",dataset="drerio_gene_ensembl")
map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = inflammation_genes$genes,
  mart = ensembl
)

# read in DE from timepoints
timepoint <- read.csv('/Users/sm2949/Desktop/bulkRNAdata/qlf_timepoint4.5vsAll.csv')

# merge with symbol
timepoint <- inner_join(
  timepoint,
  map,
  by = c("X" = "ensembl_gene_id")
)

df <- timepoint %>%
  mutate(
    neglogP = -log10(PValue),
    sig = case_when(
      FDR < 0.05 & logFC > 1  ~ "Up",
      FDR < 0.05 & logFC < -1 ~ "Down",
      TRUE                    ~ "NS"
    )
  )

# top 10 up + top 10 down
top_labels <- df %>%
  filter(sig != "NS") %>%
  arrange(FDR) %>%
  group_by(sig) %>%
  slice_head(n = 10) %>%
  ungroup()

# color palette 
volcano_colors <- c(
  "Up" = "#D55E00",
  "Down" = "#0072B2",
  "NS" = "grey70"
)

#plot
ggplot(df, aes(x = logFC, y = neglogP)) +
  geom_point(aes(color = sig), alpha = 0.8, size = 1.8) +
  scale_color_manual(values = volcano_colors) +
  
  # threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  
  # gene labels
  geom_text_repel(
    data = top_labels,
    aes(label = external_gene_name),
    size = 3.5,
    box.padding = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  
  # theme styling
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ) +
  labs(
    title = "4.5dpa Regen vs. Control GO Inflammatory Response Genes",
    x = "log2 Fold Change",
    y = "-log10(P-Value)"
  )

