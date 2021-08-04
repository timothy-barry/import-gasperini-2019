# source config file to get gasperini offsite location
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")

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
  mtx_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.exprs.mtx")
  barcodes_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.cells.txt")
  gene_ids_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.genes.txt")

  gene_expression_metadata_odm <- create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp,
                                                     barcodes_fp = barcodes_fp,
                                                     features_fp = gene_ids_fp,
                                                     file_name = "odm_gene_expression",
                                                     progress = TRUE,
                                                     on_disk_dir = processed_data_dir,
                                                     return_metadata_ondisc_matrix = FALSE)
  # Add p_mito and batch from cell_covariates data frame
  gasp_cell_covariates <- readRDS(paste0(intermediate_data_dir, "cell_covariates.rds"))
  gene_expression_metadata_odm@cell_covariates$p_mito <- gasp_cell_covariates$percent.mito
  gene_expression_metadata_odm@cell_covariates$batch <- factor(gasp_cell_covariates$prep_batch)
  saveRDS(object = gene_expression_metadata_odm,
          file = paste0(processed_data_dir, "/gene_expression_metadata_odm.rds"))
}

# gRNA count matrix
gRNA_odm_h5 <- paste0(processed_data_dir, "gRNA_odm.h5")
if (!file.exists(gRNA_odm_h5)) {
  gRNA_count_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))
  cell_barcodes <- colnames(gRNA_count_matrix)
  feature_barcodes <- data.frame(row.names(gRNA_count_matrix))
  gRNA_odm_list <- create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_count_matrix,
                                                      barcodes = cell_barcodes,
                                                      features_df = feature_barcodes,
                                                      on_disk_dir = processed_data_dir,
                                                      file_name = "gRNA_odm.h5",
                                                      return_metadata_ondisc_matrix = FALSE)
}


# the joint covariate matrix is to be computed later. Include group name and group type as columns in the feature df.
