library(matrixStats)
library(ggfortify)
library(Seurat)
library(edgeR)
library(tidyverse)
library(DT)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(EnhancedVolcano)
library(ggh4x)
library(viridis)
library(patchwork)
library(readxl)
library(SingleR)
library(RColorBrewer)
library(UpSetR)
library(data.table)
library(ggsankey)
library(ggalluvial)
library(ggrepel)
library(scales)
library(ggsci)

seuratobj <- readRDS('output/seurat_zumi_rpkm_unfilt2.rds')

## load FACS data
All_cell_assignments <- read_excel("data/All_cell_assignments_allCD34pos.xlsx") %>%
  separate(Plate_Well, c('Plate', NA,NA,'Well_ID'), remove = F) %>%
  mutate(Plate = str_remove(Plate, '^P') %>% as.numeric())

## load LINE1 data
pcr <- read.csv('output/LINE1_qPCR_corrected.csv', row.names = 1) 

## load amplicond data, new version
amplicon <- read.csv('output/genotyped_data_2025-12-02.csv') %>%
  dplyr::select(c(1:2,11:19)) %>%
  mutate(Plate2 = Plate, Well_ID = Well) %>%
  mutate(Plate = str_remove(Plate2, '^P') %>% as.numeric()) %>%
  mutate(across(matches("wt_|mut_"), ~ log1p(.x), 
                .names = "{.col}_log1p"))

master <- left_join(pcr, All_cell_assignments,
                    by = c('Plate','Well_ID')) %>%
  left_join(amplicon, by = c('Plate','Well_ID')) %>%
  mutate( # use healthy cell type
    celltype_merge = factor(celltype_Healthy,
                            levels = c("HSC","MPP","LMPP","CMP","CLP","MEP","GMP","other",'DEAD')),
    group = factor(group, 
                   levels = c("Primary_CB34_DMSO","Primary_CB34_AZA","Primary_Patient_DMSO","Primary_Patient_AZA" )))

rownames(master) <- master$Transcriptome.ID

seuratobj <- AddMetaData(seuratobj, 
                         master)

seuratobj@meta.data %>% 
  filter(Source == 'MD') %>%
  ggplot(aes(x = col, y = row, fill = nFeature_RNA, col = empty_well)) + 
  geom_tile(linewidth=1) +
  facet_wrap(.~Plate, ncol = 7, nrow = 4) +
  scale_colour_manual(values=c('TRUE'='red', 'FALSE'='black')) +
  scale_fill_viridis() +
  ggtitle('Number of genes faceted by plate') +
  # --- New additions below ---
  scale_y_discrete(limits = rev) +     # Flips the Y-axis so 'A' is at the top
  scale_x_continuous(breaks = 1:12)    # Forces X-axis ticks to show integers 1 through 12
ggsave('output/mdsl_wellqc_nfeature_facet.pdf', width = 9, height = 3)

seuratobj@meta.data %>% 
  filter(Source == 'MD') %>% 
  # Group by the grid coordinates (and empty_well if the layout is constant)
  group_by(row, col, empty_well) %>% 
  # Calculate the average for each well position
  summarise(avg_nFeature = mean(nFeature_RNA, na.rm = TRUE)) %>%
  ggplot(aes(x = col, y = row, fill = avg_nFeature, col = empty_well)) + 
  geom_tile(linewidth=1) +
  scale_colour_manual(values=c('TRUE'='red', 'FALSE'='black')) +
  scale_fill_viridis() +
  ggtitle('Average number of genes across plates') +
  # --- New additions below ---
  scale_y_discrete(limits = rev) +     # Flips the Y-axis so 'A' is at the top
  scale_x_continuous(breaks = 1:12)    # Forces X-axis ticks to show integers 1 through 12
ggsave('output/mdsl_wellqc_nfeature_avg.pdf', width = 5, height = 3)

filtered2 <- subset(seuratobj, subset = nFeature_RNA > 1000 & Source == 'MD')

filtered2@meta.data %>% dplyr::count(Source, Treatment)
# Source Treatment  n
# 1     MD       AZA 92
# 2     MD      DMSO 69

col_phase <- c(G1='#663399', G2M='#E66100', S='#2A934E')

col_treatment <- setNames(brewer.pal(2, "Set2")[1:2], c('DMSO','AZA'))

col_genotype <- c(
  "unidentifiable" = '#d9d9d9',
  "wt/wt" = '#253494',
  "mut/wt" = '#dd1c77',
  # "mut/mut" = '#980043',
  "mut/mut" = '#dd1c77'
)

filtered2 <- FindVariableFeatures(filtered2, 
                                  selection.method = "vst", 
                                  nfeatures = 2000)

filtered2 <- ScaleData(filtered2, features = rownames(filtered2))

filtered2 <- RunPCA(filtered2, features = VariableFeatures(object = filtered2))
VizDimLoadings(filtered2, dims = 1:2, reduction = "pca")

DimPlot(filtered2, reduction = "pca") 

filtered2 <- CellCycleScoring(object = filtered2, 
                              g2m.features = cc.genes.updated.2019$g2m.genes,
                              s.features = cc.genes.updated.2019$s.genes)

filtered2 <- FindNeighbors(filtered2, dims = 1:10)
filtered2 <- FindClusters(filtered2, resolution = 0.8)

filtered2 <- RunUMAP(filtered2, dims = 1:10, 
                     n.neighbors = 30)

cluster.markers <- FindAllMarkers(filtered2, only.pos = TRUE)

flat_list <- lapply(0:3, function(x) {
  cluster.markers %>% 
    dplyr::filter(cluster == x,
                  p_val_adj < 0.05) %>%
    pull(gene) %>%
    bitr(., fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db) %>%
    pull(ENTREZID)
})
names(flat_list) <- paste('cluster',0:3, sep = '_')

kegg <- compareCluster(flat_list, fun = 'enrichKEGG') %>%
  setReadable(., OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
dotplot(kegg) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
write.csv(kegg@compareClusterResult, 'output/mdsl_kegg_cc.csv')
ggsave('output/mdsl_kegg_cc.pdf', width = 6, height = 5)

gobp <- compareCluster(flat_list, fun = 'enrichGO', ont = 'BP', 
                       OrgDb = org.Hs.eg.db,
                       readable = T) %>%
  simplify()
dotplot(gobp, label_format = 50) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
write.csv(gobp@compareClusterResult, 'output/mdsl_gobp_cc.csv')
ggsave('output/mdsl_gobp_cc.pdf', width = 7, height = 5)

m_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = 'H') %>%
  dplyr::select(gs_name, ncbi_gene)
hm <- compareCluster(flat_list, fun = 'enricher', TERM2GENE=m_df) 
dotplot(hm, label_format = 50) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
write.csv(hm@compareClusterResult, 'output/mdsl_hm_cc.csv')
ggsave('output/mdsl_hm_cc.pdf', width = 7, height = 5)


selected_marker <- cluster.markers %>%
  group_by(cluster) %>% 
  arrange(p_val) %>% 
  slice_head(n = 5) %>% 
  ungroup() %>% pull(gene)

pdf('output/MDSL_umap.pdf', width = 7, height = 5)
DimPlot(filtered2, reduction = "umap", label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )

DotPlot(filtered2, features = selected_marker %>% unique())  + 
  RotatedAxis() +
  force_panelsizes(rows = unit(2, 'inch'), cols = unit(length(selected_marker)/6, 'inch')) +
  scale_color_viridis(direction = 1) +
  ggtitle('Top 5 markers in each cluster')

DimPlot(filtered2, reduction = "umap",  group.by = 'Phase',
        cols = col_phase) + 
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )

DimPlot(filtered2, reduction = "umap",  group.by = 'Treatment',
        cols = col_treatment) + 
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )

DimPlot(filtered2, reduction = "umap",
        group.by = 'genotype_srsf2',
        cols = col_genotype) +
  plot_layout(ncol = 1, nrow = 1,
              widths = unit(3, "in"),
              heights = unit(3, "in") )

DimPlot(filtered2, reduction = "umap",
        group.by = 'genotype_tet2a',
        cols = col_genotype) +
  plot_layout(ncol = 1, nrow = 1,
              widths = unit(3, "in"),
              heights = unit(3, "in") )

DimPlot(filtered2, reduction = "umap",
        group.by = 'genotype_tet2b',
        cols = col_genotype) +
  plot_layout(ncol = 1, nrow = 1,
              widths = unit(3, "in"),
              heights = unit(3, "in") )

patchwork::wrap_plots( FeaturePlot( filtered2,
                                    features=c("dCt_adj"),
                                    combine=FALSE)) +
  plot_layout(ncol = 1, nrow = 1,
              widths = unit(c(3), "in"),
              heights = unit(3, "in") ) &
  theme(legend.position = "right") &
  scale_color_viridis_c(direction = -1,
                        limits = c(10.5, 16),
                        oob = squish, # squish values outside limits
                        na.value = "grey90") &
  ggtitle('UMAP of LINE1 and RNAseq passed cells')

dev.off()

filtered2@meta.data %>% write.csv('output/mdsl_filtered_metadata.csv')
