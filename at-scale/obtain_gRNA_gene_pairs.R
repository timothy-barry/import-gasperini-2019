#######################
# 0. Load data; set fps
#######################
require(magrittr)

gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/")
raw_dir <- paste0(gasp_offsite, "raw/")
processed_dir <- paste0(gasp_offsite, "processed/")
processed_dir_grouped <- paste0(processed_dir, "gRNA_grouped/")
processed_dir_ungrouped <- paste0(processed_dir, "gRNA_ungrouped/")

# read gRNA-gene pairs and gRNA info
all_results <- readr::read_tsv(file = paste0(raw_dir, "GSE120861_all_deg_results.at_scale.txt"),
                               col_names = TRUE, col_types = "c")
gRNA_groups <- readr::read_tsv(file = paste0(raw_dir, "GSE120861_grna_groups.at_scale.txt"),
                               col_names = FALSE, col_types = "c") %>% dplyr::rename("gRNA_group" = "X1", "barcode" = "X2")
#######################
# 1. Save grouped pairs
#######################
grouped_pairs <- all_results %>% 
  dplyr::select(gene_id = ENSG, gRNA_id = gRNA_group, gRNA_quality = quality_rank_grna, site_type) %>%
  dplyr::mutate_all(factor)
saveRDS(grouped_pairs, paste0(processed_dir_grouped, "pairs_grouped.rds"))

#########################
# 2. Save ungrouped pairs
#########################
pairs_to_analyze <- all_results %>%
  dplyr::select(gene_id = ENSG, gRNA_group, gRNA_quality = quality_rank_grna, site_type)
pairs_to_analyze_expanded <- dplyr::left_join(x = pairs_to_analyze, y = gRNA_groups, by = "gRNA_group") %>% 
  dplyr::mutate_all(factor)
saveRDS(pairs_to_analyze_expanded, paste0(processed_dir_ungrouped, "pairs_ungrouped.rds"))
