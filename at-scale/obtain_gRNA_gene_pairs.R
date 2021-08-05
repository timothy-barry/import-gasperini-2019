# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/")
raw_dir <- paste0(gasp_offsite, "raw/")
processed_dir <- paste0(gasp_offsite, "processed/")

# read gRNA-gene pairs and gRNA info
all_results <- readr::read_tsv(file = paste0(raw_dir, "GSE120861_all_deg_results.at_scale.txt"), col_names = TRUE, col_types = "c")
gRNA_groups <- readr::read_tsv(file = paste0(raw_dir, "GSE120861_grna_groups.at_scale.txt"), col_names = FALSE, col_types = "c") %>% dplyr::rename("gRNA_group" = "X1", "barcode" = "X2")

# obtain the pairs to analyze from the results
pairs_to_analyze <- all_results %>% dplyr::select(gene_id = ENSG, gRNA_group, gRNA_quality = quality_rank_grna, site_type)

# We see that the gene-gRNA pairs of "pairs_to_analyze" are unique. However, the gRNA groups are duplicated in "gRNA_groups."
pairs_to_analyze %>% dplyr::summarize(paste0(gene_id, "-", gRNA_group)) %>% dplyr::pull() %>% duplicated() %>% any()
gRNA_groups$gRNA_group %>% duplicated() %>% any()

# merge the gRNA_groups with all_results
pairs_to_analyze_expanded <- dplyr::left_join(x = pairs_to_analyze, y = gRNA_groups, by = "gRNA_group")

# finally, convert to factors
pairs_to_analyze_expanded <- dplyr::mutate_all(pairs_to_analyze_expanded, factor)

# save pairs to analyze
saveRDS(object = pairs_to_analyze_expanded, file = paste0(processed_dir, "gRNA_gene_pairs.rds"))
