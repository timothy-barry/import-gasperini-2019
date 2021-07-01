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

#### check inconsistency of specific gRNA and cell
if (FALSE) {
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
raw_data_dir <- paste0(gasp_offsite, "raw/")
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")
library(magrittr)
gRNA_barcodes <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.at_scale.txt"),
                                        col_names = c("group_name", "gRNA_barcode"), col_types = "cc")
gRNA_counts <- readr::read_tsv(paste0(raw_data_dir, "all_libraries.gRNAcaptured.aggregated.txt"),
                               col_names = TRUE,
                               col_types = "cccc") %>% dplyr::rename(cell_barcode = cell, gRNA_barcode = barcode)
gRNA_indicator_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_indicators.rds"))

gRNA_grp <- "chr10.2575_top_two"
gRNA_grp_barcodes <- dplyr::filter(gRNA_barcodes, group_name == gRNA_grp) %>% dplyr::pull(gRNA_barcode)
cell_barcode <- "GTCAAGTTCAGCGACC-1_1B_2_SI-GA-F3"
gRNA_counts_my_cell <- gRNA_counts[gRNA_counts$cell_barcode == gsub('.{9}$', '', cell_barcode),]

gRNA_counts_my_cell[gRNA_counts$gRNA_barcode %in% gRNA_grp_barcodes,] 
gRNA_indicator_matrix[cell_barcode, gRNA_grp] # no gRNA in indicator matrix; inconsistency
}