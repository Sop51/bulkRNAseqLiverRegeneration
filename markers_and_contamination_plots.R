library(tidyverse)
library(ggplot2)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(pheatmap)

# read in the normalized counts matrix
cpms <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/normalized_counts.csv", row.names = 1)
colnames(cpms) <- gsub("^X", "", colnames(cpms)) # fix format

# define the lists of markers 
liver_markers   <- c("fabp10a", "cp", "apoa2")
pancreas_markers <- c("pdx1", "ins", "nkx6.1")
heart_markers <- c("myh6", "myl7", "vmhc")
spleen_kidney_markers <- c("gata1a", "pax5", "trbc1")

# combine
all_marker_genes <- c(liver_markers, pancreas_markers, 
                      heart_markers, spleen_kidney_markers)

# map to ensembl/symbol format
gene_map <- AnnotationDbi::select(org.Dr.eg.db,
                                  keys = rownames(cpms),
                                  columns = c("ENSEMBL", "SYMBOL"),
                                  keytype = "ENSEMBL")

# subset CPMs to genes that have a symbol in marker list
plot_mat <- cpms %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  inner_join(gene_map %>% dplyr::select(ENSEMBL, SYMBOL), by = "ENSEMBL") %>%
  filter(SYMBOL %in% all_marker_genes) %>%
  column_to_rownames("SYMBOL") %>%
  dplyr::select(-ENSEMBL) %>%
  as.matrix()

# add missing markers as rows of zeros (no expression if not mapped)
missing <- setdiff(all_marker_genes, rownames(plot_mat))
if(length(missing) > 0){
  zero_rows <- matrix(0, nrow = length(missing), ncol = ncol(plot_mat),
                      dimnames = list(missing, colnames(plot_mat)))
  plot_mat <- rbind(plot_mat, zero_rows)
}

# order rows for plotting
gene_order <- c(liver_markers, pancreas_markers, heart_markers, spleen_kidney_markers)
plot_mat <- plot_mat[gene_order, ]

# create the annotation df
anno_row <- data.frame(
  Marker_For = factor(c(rep("Liver", length(liver_markers)),
                        rep("Pancreas", length(pancreas_markers)),
                        rep("Heart", length(heart_markers)),
                        rep("Spleen/Kidney", length(spleen_kidney_markers))),
                      levels = c("Liver", "Pancreas", "Heart", "Spleen/Kidney")),
  row.names = gene_order
)

# colors for annotations
anno_colors <- list(
  Marker_For = c(Liver = "#A53247", Pancreas = "#F4F1DE", 
                 Heart = "#81B29A", `Spleen/Kidney` = "#3D405B")
)

# plot heatmap
pheatmap(plot_mat,
         color = colorRampPalette(c("#013B7D", "white", "#FFA43C"))(100),
         cluster_rows = FALSE,         
         cluster_cols = FALSE,       
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "column",
         annotation_row = anno_row,     
         annotation_colors = anno_colors,
         annotation_legend = TRUE,
         annotation_names_row = FALSE,
         fontsize_row = 10,
         fontsize_col = 9,
         border_color = NA,
         main = "Cross-Tissue Contamination Marker Heatmap",
         legend = TRUE,
         width = 10, height = 9)
