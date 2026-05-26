calculate_amplicon_prop <- function(x) {
  
  props2 <- filtered2@meta.data %>%
    # collapse to a single TET2 genotype label
    mutate(
      genotype_tet2 = case_when(
        genotype_tet2a == "wt/wt"  & genotype_tet2b == "wt/wt"  ~ "wt/wt",
        genotype_tet2a == "mut/wt" & genotype_tet2b == "mut/wt" ~ "mut/wt",
        TRUE ~ "unidentifiable"
      )
    ) %>%
    # split `group` safely; keep original in case
    tidyr::separate(
      group, into = c(NA, "Cell_Type", "Treatment"),
      sep = "_", remove = FALSE, fill = "right", extra = "merge"
    ) %>%
    filter(
      !celltype_merge %in% c("DEAD", "other"),
      !is.na(Cell_Type), !is.na(Treatment), !is.na(celltype_merge)
    ) %>%
    # count AFTER collapsing genotype, so each row is a single genotype bin
    count(Cell_Type, celltype_merge, Treatment, genotype_tet2, name = "n") %>%
    group_by(Cell_Type, celltype_merge, Treatment) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    filter(total > 0) %>%                # drop panels with no cells (avoids 0/0)
    mutate(
      prop = n / total,
      Treatment      = factor(Treatment, levels = treat_levels),
      genotype_tet2  = factor(genotype_tet2, levels = geno_levels)
    ) %>%
    # ensure missing genotype bins show as 0-height bars, but keep only non-empty panels
    group_by(Cell_Type, celltype_merge, Treatment) %>%
    complete(genotype_tet2, fill = list(n = 0, prop = 0)) %>%
    ungroup()
  
  return(props2)
}