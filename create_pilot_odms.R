# step 3 convert everything to ondisk matrix ###
# source config file to get gasperini offsite location
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
# gasp_offsite <- c("/Users/sijia_work/Documents/gasperini-2019/")

# create the processed data directory; set the raw and intermediate directories
processed_data_dir <- paste0(gasp_offsite, "processed/")
if (!dir.exists(processed_data_dir)) dir.create(path = processed_data_dir, recursive = TRUE)
raw_data_dir <- paste0(gasp_offsite, "raw/")
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")

# load ondisc
library(ondisc)

# gene count matrix
gene_odm_to_save_fp <- paste0(processed_data_dir, "gene_expression_odm.rds")
if (!file.exists(gene_odm_to_save_fp)) {
  mtx_fp <- paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.exprs.mtx")
  barcodes_fp <- paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.cells.txt")
  gene_ids_fp <- paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.genes.txt")

  gene_expression_metadata_odm <- create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp,
                                                     barcodes_fp = barcodes_fp,
                                                     features_fp = gene_ids_fp,
                                                     file_name = "odm_gene_expression",
                                                     progress = TRUE,
                                                     on_disk_dir = processed_data_dir,
                                                     return_metadata_ondisc_matrix = FALSE)
  gasp_cell_covariates <- readRDS(paste0(intermediate_data_dir, "cell_covariates.rds"))
  saveRDS(object = gene_expression_metadata_odm,
          file = paste0(processed_data_dir, "/gene_expression_metadata_odm.rds"))
}