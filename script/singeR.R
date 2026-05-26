
# umi-count from dev cell 2022 paper
# 1. Load libraries
library(SummarizedExperiment)
library(scuttle)            # for logNormCounts (optional but recommended)

# 2. Read your expression matrix
#    - assume first column is gene IDs, headers are sample names
expr_mat <- read.delim(
  "data/dev_cell_2022/hspc_raw_matrix.txt",
  header    = TRUE,
  row.names = 1,     # gene IDs in col 1
  check.names = FALSE
)

# 3. (Optional) if your matrix came in as samples × genes, transpose:
# expr_mat <- t(expr_mat)

# 4. Read in your sample‐to‐label table
#    assume a two‐column file with “sample” and “label”
sample_df <- read.delim(
  "data/dev_cell_2022/Cell_metadata.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)
# make sure rownames match colnames(expr_mat)
rownames(sample_df) <- sample_df$cell
# sample_df <- sample_df[colnames(expr_mat), , drop=FALSE]
stopifnot(all.equal(colnames(expr_mat), rownames(sample_df)))

# 5. Build the SummarizedExperiment
ref_se <- SummarizedExperiment(
  assays = list(counts = as.matrix(expr_mat)),
  colData = sample_df
)

# 6. Log‐normalize (so you have a “logcounts” assay)
ref_se <- logNormCounts(ref_se)
# now assays(ref_se) has both “counts” and “logcounts”

# 7. Double‐check
ref_se
#   > class: SummarizedExperiment 
#   > dim: <nrow(expr_mat)> <ncol(expr_mat)> 
#   > assays(0): counts logcounts
#   > rownames(0): ...
#   > colnames(0): your sample names
#   > colData(0):  data.frame with your “label” column

saveRDS(ref_se, 'data/dev_cell_2022/ref_se.rds')

# 8. Pass into SingleR
library(SingleR)
pred <- SingleR(
  test            = filtered2[['RNA']]$counts,
  ref             = ref_se,
  labels          = colData(ref_se)$RNA_clusters,
  de.method       = "wilcox",
  # de.n            = 50, 
  assay.type.ref  = "logcounts",
  assay.type.test = "logcounts"
)

sankey_df <- data.frame(
  cluster = filtered2$RNA_snn_res.0.8,
  # cluster2 = filtered2$RNA_snn_res.0.8,
  devcell = pred$labels,
  FACS = filtered2$celltype_merge,
  rowname = colnames(filtered2), 
  source = filtered2$Cell_Type, 
  treatment = filtered2$Treatment,
  group = filtered2$group) 

## compare dev cell vs transcriptome based clusters
sankey_df %>%
  filter(source != 'MD') %>%
  group_by(source, treatment) %>%
  group_modify(~make_long(.x, cluster:devcell)) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node),
             label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) 

sankey_df %>%
  mutate(devcell = factor(devcell, 
                          levels = c("HSC","MPP1","MPP2","LMPP","MLP","MEP","Ery","Neu1","Mk","MD","EBM","B_NK1","Neu2"))) %>%
  ggplot(aes(x = cluster)) +
  geom_bar(aes(fill = devcell), position = "fill") +
  # scale_fill_manual(
  #   values   = col_facs,
  #   na.value = na_grey,
  #   name     = "Cell type"
  # ) +
  scale_fill_manual(values = pal13 <- viridis(13, option = "C")  ) +  
  labs(x = "Cluster", y = "Proportion") +
  force_panelsizes(rows = unit(4, 'inch'), cols = unit(3, 'inch')) 

sankey_df %>%
  filter(source != 'MDSL') %>%
  group_by(source, treatment) %>%
  group_modify(~make_long(.x, group, cluster, devcell)) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node),
             label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) 


ggplot(data = sankey_df %>%
         filter(source != 'MDSL') ) +
  geom_mosaic(aes(x = product(group,
                              devcell,
                              cluster
                              ), 
                  fill=group), divider=mosaic("h")) 
