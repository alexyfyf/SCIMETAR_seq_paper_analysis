
# SCIMETAR-seq data analysis
## Single-cell RNA-seq Alignment and Quantification
Code used to perform analysis can be found in https://github.com/alexyfyf/SCIMETAR. In brief, raw FASTQ files from single-cell RNA-seq were aligned to the GRCh38 reference genome using STAR (v2.7.11b) in 2-pass mode (--twopassMode Basic) [[1]](#ref-1). A splice junction overhang of 150 bp (--sjdbOverhang 150) was used and we also set --outSAMmultNmax 1 and --outFilterMultimapNmax 50 as used by zUMIs pipeline [[2]](#ref-2). Gene expression quantification was performed on the aligned BAM files using featureCounts (v2.0.8) with Gencode v44 annotation, parameters -p (for paired-end data) and -s 0 (unstranded) [[3]](#ref-3). The output was a gene-by-cell count matrix.
## Targeted Amplicon Sequencing Processing
### Demultiplexing:
Raw amplicon data was demultiplexed using Matchbox (v0.1.0) [[4]](#ref-4). Reads were assigned to samples based on unique combinations of forward and reverse barcodes identified at the 5' ends of Read 1 and the 3' ends of Read 2 (reverse complement), allowing for a 1 bp mismatch. 
### Genotyping:
To genotype each sequencing fragments, we perform mutation calling separately in Read1 (R1) and Read2 reverse complement (R2rc) in 2 rounds. In the first round, we look for an exact 5bp match centered around the mutation site (table below) in both R1 and the R2rc for TET2a, TET2b and SRSF2. Any reads that didn't match in round 1 are searched again specifically for SRSF2, using a 10bp flanking sequences and allowing for a 1bp mismatch.
For each read, we then genotype them as
(1) TET2a WT or MUT only if the sequence was found in R1 alone, or in both R1 and R2rc. If it was only found in R2rc, it is discarded.
(2) TET2b WT or MUT only if the sequence was found in R2rc alone, or in both R1 and R2rc.
(3) SRSF2 WT or MUT if the sequence was found in R1 alone, both R1 and R2rc, or for in vivo data only if it was recovered during the less stringent Round 2 search. 

| Gene | Name | Mutation (hg19) | hg38 liftover |
| --- | --- | --- | --- |
| SRSF2 | SRSF2 | SRSF2_chr17:74732959_G>T | chr17:76736877 |
| TET2 | TET2a | TET2_chr4:106158367_AAGAC>A | chr4:105237210 |
| TET2 | TET2b | TET2_chr4:106164802_GC>G | chr4:105243645 |

To determine the genotype of each cell for a particular variant, cells with less than 20 or 10 genotyped reads, for in vitro and in vivo respectively, are labelled "unidentifiable". We then calculated the ratio of mutant reads / (mutant + wildtype reads), and if the ratio is > 0.99, the cell is labelled "mut/mut", if the ratio is < 0.05, the cell is labelled "wt/wt", and rest labelled "mut/wt".
## Single-cell RNA-seq Data Analysis
The analysis was performed using R version 4.4.2. The gene expression count matrix was loaded into the Seurat R package (v5.1.0) and RPKM normalized [[5]](#ref-5). Cells were filtered to retain high-quality cells. Cells were included if they had more than 1000, 2000 or 1000 detected features (for MDSL cell line, in vitro and in vivo respectively) and were not marked as empty wells or dead cells based on matched FACS index sorting data. Raw counts were normalized and variance stabilized using the vst selection method in Seurat to identify the top 2,000 variable features. Data was scaled, and Principal Component Analysis (PCA) was performed on the variable features. The top 20 (10 for MDSL) principal components were used for downstream clustering and dimensionality reduction. A Shared Nearest Neighbor (SNN) graph was constructed with k=50 (30 for MDSL) neighbors, and clusters were identified using the Louvain algorithm with a resolution of 1.0 (0.8 for MDSL). Uniform Manifold Approximation and Projection (UMAP) was calculated on the same PCA dimensions to visualize the data in two dimensions. Cell cycle scoring was performed to assess the G2M and S phase scores for each cell.
## Cell Type Annotation
Automated annotation was performed using a reference dataset from reference paper (Dev Cell 2022) using Azimuth (v0.5.0) to align with established hematopoietic hierarchies [[6]](#ref-6), [[7]](#ref-7).
Cluster-specific marker genes were identified using the FindAllMarkers function in Seurat with a Wilcoxon Rank Sum test. The resulting markers were compared against canonical markers for validation (e.g., AVP, CRHBP, DLK1 for HSC/MPP; MPO, CDK6 for Neu1; EBM markers). Pathway analysis of marker gene was performed using clusterProfiler package and Gene Ontology database [[8]](#ref-8), [[9]](#ref-9).
## Genotype-Phenotype Integration
Targeted amplicon sequencing data for SRSF2, TET2A, and TET2B mutations was integrated with the single-cell transcriptomic data. Genotyping results were matched to single cells using plate and well information. Cells were categorized based on their mutational status (e.g., WT/WT, Mut/WT, Mut/Mut).
We correlated genotypes with cell phenotypes by constructing Sankey/Alluvial plots to visualize the flow between FACS immunophenotypes, Seurat clusters, and reference-based annotations across different experimental conditions (e.g., DMSO vs. AZA treatment).
## LINE-1 qPCR Analysis
LINE-1 retrotransposon expression was quantified using qPCR. Cycle threshold (Ct) values were obtained for the target (methylated 5’ of LINE-1, FAM channel) and a reference LINE-1 ORF2 sequence (HEX channel).
Quality control was performed by analyzing the distribution of reference gene Ct values (HEX_Ref) across plates and experimental groups. Outlier wells were identified and removed using a Median Absolute Deviation (MAD) approach, retaining wells where HEX_Ref was within the median ± 2 * MAD (group-wise for in vitro and plate-wise for in vivo).
For the past wells,  dCT was calculated as Ct (FAM_Tar) – Ct (Hex_Ref).
For the in vivo dataset, we use MspI wells as internal control to normalize across plates and remove batch effect, limma::removeBatchEffect was used, and we used plate as batch and MspI wells as design matrix [[10]](#ref-10).

## Publication
Bhuyan, G. S., Yan, A., Nguyen, M. N. T., Zou, X., Gullapalli, V., Vaughan, L., Stonehouse, O., Hampton, H. R., Shen, S., Truong, P., Dissanayake, R., Ghodousi, E. S., Joshi, S., Koch, F. C., Chung, H. M., Zanini, F., Vafaee, F., Huang, Y., Thoms, J. A. I., Faridani, O., Jolly, C. J. & Pimanda, J. E. **SCIMETAR-seq tracks immunophenotype, demethylation, mutations, and transcriptomes in single cells undergoing HMA therapy.** *bioRxiv* 2026.04.26.720516 (2026). https://doi.org/10.64898/2026.04.26.720516

## References
1. <a id="ref-1"></a>Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29, 15–21 (2013).
2. <a id="ref-2"></a>Parekh, S., Ziegenhain, C., Vieth, B., Enard, W. & Hellmann, I. zUMIs - A fast and flexible pipeline to process RNA sequencing data with UMIs. *Gigascience* 7, giy059 (2018).
3. <a id="ref-3"></a>Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30, 923–930 (2014).
4. <a id="ref-4"></a>Schuster, J. et al. Powerful read processing with matchbox. 2025.11.09.685711 Preprint at https://doi.org/10.1101/2025.11.09.685711 (2026).
5. <a id="ref-5"></a>Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol* 42, 293–304 (2024).
6. <a id="ref-6"></a>Hao, Y. et al. Integrated analysis of multimodal single-cell data. *Cell* 184, 3573-3587.e29 (2021).
7. <a id="ref-7"></a>Zhang, Y. et al. Temporal molecular program of human hematopoietic stem and progenitor cells after birth. *Developmental Cell* 57, 2745-2760.e6 (2022).
8. <a id="ref-8"></a>Yu, G., Wang, L. G., Han, Y. & He, Q. Y. ClusterProfiler: An R package for comparing biological themes among gene clusters. *OMICS A Journal of Integrative Biology* 16, 284–287 (2012).
9. <a id="ref-9"></a>The Gene Ontology Consortium. Expansion of the Gene Ontology knowledgebase and resources. *Nucleic Acids Research* 45, D331–D338 (2017).
10. <a id="ref-10"></a>Ritchie, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic acids research* 43, e47 (2015).
