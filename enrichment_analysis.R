library(clusterProfiler)
library(org.Dr.eg.db)
library(patchwork)
library(tidyverse)
library(enrichplot)

# set organism
organism = "org.Dr.eg.db"

# read in the DE results
timepoint0vsAll <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/qlf_timepoint0vsAll.csv", row.names = 1)
timepoint1.5vsAll <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/qlf_timepoint1.5vsAll.csv", row.names = 1)
timepoint4.5vsAll <- read.csv("/Users/sm2949/Desktop/bulkRNAdata/qlf_timepoint4.5vsAll.csv", row.names = 1)

# define timepoint to use
df <- timepoint0vsAll

# create ranked list
gene_df <- df %>%
  mutate(ENSEMBL = rownames(df)) %>%
  select(ENSEMBL, logFC) %>%
  filter(!is.na(logFC))

# convert ensembl to entrez
ids <- bitr(
  gene_df$ENSEMBL,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = organism
)

# join ensembl and entrez
gene_df2 <- gene_df %>%
  inner_join(ids, by = "ENSEMBL")

# collapse duplicates
gene_df2 <- gene_df2 %>%
  group_by(ENTREZID) %>%
  summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop")

# prepare gene list, rank
kegg_gene_list <- gene_df2$logFC
names(kegg_gene_list) <- gene_df2$ENTREZID
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

# run enrichment
kk2 <- gseKEGG(
  geneList      = kegg_gene_list,
  organism      = "dre",
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "none",
  keyType       = "kegg",
  scoreType = 'pos'
)

# network plot
cnetplot(kk2, categorySize="pvalue", foldChange=kegg_gene_list)

# pathway plot
#dme <- pathview(gene.data=kegg_gene_list, pathway.id="dre04350", species = "dre")

kk2 <- pairwise_termsim(kk2)

emapplot(kk2)


# run on go terms ----

# define timepoint to use
df <- timepoint0vsAll

# subset to only positive and significant
df <- df[(df$logFC > 1 & df$FDR < 0.05), ]

gene <- bitr(
  rownames(df),
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = organism
)

ego <- enrichGO(gene = gene$ENTREZID,
                OrgDb = organism,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

ego <- pairwise_termsim(ego)

emapplot(ego)


