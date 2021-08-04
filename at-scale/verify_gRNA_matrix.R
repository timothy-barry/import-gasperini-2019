# test correctness of gRNA matrix
require(dplyr)

# set directories
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")
raw_data_dir <- paste0(gasp_offsite, "raw/")

# load gRNA barcodes, count matrix, and indicator matrix
gRNA_barcodes_in_use <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.at_scale.txt"),
                                        col_names = c("group_name", "gRNA_barcode"), col_types = "cc")
gRNA_count_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))
gRNA_indicator_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_indicators.rds"))

# check equality between thresholded count matrix and indicator matrix
set.seed(10)
grp_names <- sample(gRNA_barcodes_in_use$group_name, 50) 

# verify that our nonzero entries are a superset of Gasperini's
for (grp_name in grp_names) {
  barcodes <- gRNA_barcodes_in_use %>% dplyr::filter(group_name %in% grp_name) %>% dplyr::pull(gRNA_barcode)
  our_threshold <- (Matrix::colSums(gRNA_count_matrix[barcodes,] >= 5) >= 1)
  gasp_threshold <- gRNA_indicator_matrix[,grp_name]
  our_threshold_nonzero <- which(our_threshold)
  gasp_threshold_nonzero <- which(gasp_threshold)
  print(all(gasp_threshold_nonzero %in% our_threshold_nonzero))
}
