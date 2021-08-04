# source config file to get gasperini offsite location
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")

# create the intermediate data directory; set raw directory
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")
if (!dir.exists(intermediate_data_dir)) dir.create(path = intermediate_data_dir, recursive = TRUE)
raw_data_dir <- paste0(gasp_offsite, "raw/")

# Obtain binary gRNA matrix and cell metadata from monocole object
library(monocle)
library(magrittr)
monocle_obj <- readRDS(paste0(raw_data_dir, "/GSE120861_at_scale_screen.cds.rds"))
cell_metadata <- pData(monocle_obj)
rm(monocle_obj); gc()
covariates_cols <- 1:18
gRNA_cols <- 19:ncol(cell_metadata)

gRNA_indicators <- cell_metadata[,gRNA_cols]
cell_covariates <- cell_metadata[,covariates_cols]

# save the gRNA indicators and cell covariates
saveRDS(gRNA_indicators, paste0(intermediate_data_dir, "gRNA_indicators.rds"))
saveRDS(cell_covariates, paste0(intermediate_data_dir, "cell_covariates.rds"))

# load the gRNA "groups" at scale and gRNA count matrix
cell_barcodes_in_use_long <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_at_scale_screen.cells.txt"),
                                        col_names = FALSE, col_types = "c") %>% dplyr::pull()
# all cell barcodes have 32 characters; strip the last 9
cell_barcodes_in_use <- gsub('.{9}$', '', cell_barcodes_in_use_long)
gRNA_barcodes_in_use <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.at_scale.txt"),
                                        col_names = c("group_name", "gRNA_barcode"), col_types = "cc")

gRNA_counts <- readr::read_tsv(paste0(raw_data_dir, "all_libraries.gRNAcaptured.aggregated.txt"),
                               col_names = TRUE,
                               col_types = "cccc") %>% dplyr::rename(cell_barcode = cell, gRNA_barcode = barcode)
# keep rows with cells in use, gRNA barcode in use, and nonzero UMI counts
gRNA_counts_sub <- gRNA_counts %>% dplyr::filter(umi_count > 0,
                                                 cell_barcode %in% cell_barcodes_in_use,
                                                 gRNA_barcode %in% gRNA_barcodes_in_use$gRNA_barcode)
# assign integer labels to the cell_barcodes and gRNA_barcodes
cell_idxs <- match(x = gRNA_counts_sub$cell_barcode, table = cell_barcodes_in_use)
gRNA_idxs <- match(x = gRNA_counts_sub$gRNA_barcode, table = gRNA_barcodes_in_use$gRNA_barcode)

m <- Matrix::sparseMatrix(i = gRNA_idxs,
                          j = cell_idxs,
                          x = as.integer(gRNA_counts_sub$umi_count),
                          dims = c(length(unique(gRNA_barcodes_in_use$gRNA_barcode)), length(unique(cell_barcodes_in_use))),
                          repr = "T")
row.names(m) <- gRNA_barcodes_in_use$gRNA_barcode
colnames(m) <- cell_barcodes_in_use_long

# save the count matrix to the intermediate file directory
saveRDS(object = m, file = paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))

