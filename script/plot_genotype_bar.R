plot_genotype_bar <- function(x, gene = "srsf2") {
  # tidy-eval setup
  gcol_name <- paste0("genotype_", tolower(gene))
  gsym <- rlang::sym(gcol_name)
  
  df_srsf2 <- x %>%
    # split `group` safely; keep original in case
    tidyr::separate(
      group, into = c(NA, "Cell_Type", "Treatment"),
      sep = "_", remove = FALSE, fill = "right", extra = "merge"
    ) %>% 
    filter(.data[[gcol_name]] != "unidentifiable")
  
  title_facets <- paste(toupper(gene), "mutation in FACS subtype")
  title_clusters <- paste(toupper(gene), "mutation in Seurat clusters")
  title_devcell <- paste(toupper(gene), "mutation in Dev Cell annotation")
  
  # ---------- FACS ----------
  props <- df_srsf2 %>%
    dplyr::filter(
      !celltype_merge %in% c("DEAD", "other"),
      !is.na(Cell_Type), !is.na(Treatment), !is.na(celltype_merge)
    ) %>%
    # count AFTER collapsing genotype, so each row is a single genotype bin
    dplyr::count(Cell_Type, celltype_merge, Treatment, !!gsym, name = "n") %>%
    dplyr::group_by(Cell_Type, celltype_merge, Treatment) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(total > 0) %>%
    dplyr::mutate(
      prop = n / total,
      Treatment = factor(Treatment, levels = treat_levels),
      !!gsym := factor(.data[[gcol_name]], levels = geno_levels)
    ) %>%
    dplyr::group_by(Cell_Type, celltype_merge, Treatment) %>%
    tidyr::complete(!!gsym, fill = list(n = 0, prop = 0)) %>%
    dplyr::ungroup()
  
  p1 <- ggplot(props, aes(x = Treatment, y = prop, fill = .data[[gcol_name]])) +
    geom_col(width = 0.8) +
    labs(x = "Cell type", y = "Proportion", fill = "Genotype", title = title_facets) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_nested(. ~ Cell_Type + celltype_merge, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = col_genotype, name = "Genotype")
  
  p2 <- ggplot(props, aes(x = Treatment, y = n, fill = .data[[gcol_name]])) +
    geom_col(width = 0.8) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5), size = 3.5, col = "white") +
    labs(x = "Cell type", y = "Number of cells", fill = "Genotype", title = title_facets) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_nested(. ~ Cell_Type + celltype_merge, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = col_genotype, name = "Genotype")
  
  # ---------- Seurat clusters ----------
  props <- df_srsf2 %>%
    dplyr::count(Cell_Type, seurat_clusters, Treatment, !!gsym, name = "n") %>%
    dplyr::group_by(Cell_Type, seurat_clusters, Treatment) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(total > 0) %>%
    dplyr::mutate(
      prop = n / total,
      Treatment = factor(Treatment, levels = treat_levels),
      !!gsym := factor(.data[[gcol_name]], levels = geno_levels)
    ) %>%
    dplyr::group_by(Cell_Type, seurat_clusters, Treatment) %>%
    tidyr::complete(!!gsym, fill = list(n = 0, prop = 0)) %>%
    dplyr::ungroup() %>%
    mutate(seurat_clusters = factor(seurat_clusters, levels = ord))
  
  p3 <- ggplot(props, aes(x = Treatment, y = prop, fill = .data[[gcol_name]])) +
    geom_col(width = 0.8) +
    labs(x = "Cell type", y = "Proportion", fill = "Genotype", title = title_clusters) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_nested(. ~ Cell_Type + seurat_clusters, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = col_genotype, name = "Genotype")
  
  p4 <- ggplot(props, aes(x = Treatment, y = n, fill = .data[[gcol_name]])) +
    geom_col(width = 0.8) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5), size = 3.5, col = "white") +
    labs(x = "Cell type", y = "Number of cells", fill = "Genotype", title = title_clusters) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_nested(. ~ Cell_Type + seurat_clusters, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = col_genotype, name = "Genotype")
  
  # ---------- DevCell ----------
  props <- df_srsf2 %>%
    dplyr::count(Cell_Type, SingleR.DevCell, Treatment, !!gsym, name = "n") %>%
    dplyr::group_by(Cell_Type, SingleR.DevCell, Treatment) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(total > 0) %>%
    dplyr::mutate(
      prop = n / total,
      Treatment = factor(Treatment, levels = treat_levels),
      !!gsym := factor(.data[[gcol_name]], levels = geno_levels)
    ) %>%
    dplyr::group_by(Cell_Type, SingleR.DevCell, Treatment) %>%
    tidyr::complete(!!gsym, fill = list(n = 0, prop = 0)) %>%
    dplyr::ungroup() %>%
    mutate(SingleR.DevCell = factor(SingleR.DevCell, levels = names(col_ref)))
  
  p5 <- ggplot(props, aes(x = Treatment, y = prop, fill = .data[[gcol_name]])) +
    geom_col(width = 0.8) +
    labs(x = "Cell type", y = "Proportion", fill = "Genotype", title = title_devcell) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_nested(. ~ Cell_Type + SingleR.DevCell, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = col_genotype, name = "Genotype")
  
  p6 <- ggplot(props, aes(x = Treatment, y = n, fill = .data[[gcol_name]])) +
    geom_col(width = 0.8) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5), size = 3.5, col = "white") +
    labs(x = "Cell type", y = "Number of cells", fill = "Genotype", title = title_devcell) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_nested(. ~ Cell_Type + SingleR.DevCell, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = col_genotype, name = "Genotype")
  
  return(p = list(p1,p2,p3,p4,p5,p6))
}
