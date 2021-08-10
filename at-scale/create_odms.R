# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/")

# create the processed data directory; set the raw and intermediate directories
processed_data_dir <- paste0(gasp_offsite, "processed/")
if (!dir.exists(processed_data_dir)) dir.create(path = processed_data_dir, recursive = TRUE)
raw_data_dir <- paste0(gasp_offsite, "raw/")
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")

# load ondisc
library(ondisc)

# gene count matrix
mtx_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.exprs.mtx")
barcodes_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.cells.txt")
gene_ids_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.genes.txt")

odm_fp <- paste0(processed_data_dir, "gasp_scale_gene_expressions.odm")
metadata_fp <- paste0(processed_data_dir, "gasp_scale_gene_metadata.rds")

# create the odm 
gene_odm <- create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp, barcodes_fp = barcodes_fp, 
                                          features_fp = gene_ids_fp, odm_fp = odm_fp,
                                          metadata_fp = metadata_fp, progress = TRUE)

# Add p_mito and batch from cell_covariates data frame
gasp_cell_covariates <- readRDS(paste0(intermediate_data_dir, "cell_covariates.rds"))
gene_odm_plus_pmito_batch <- mutate_cell_covariates(gene_odm, p_mito = gasp_cell_covariates$percent.mito,
                                                    batch = factor(gasp_cell_covariates$prep_batch))

# save the metadata (overwriting the original metadata file)
save_odm(odm = gene_odm_plus_pmito_batch,
         metadata_fp = metadata_fp)


# next, load the gRNA count matrix. Write the backing .odm file and metadata file to disk.
odm_fp <- paste0(processed_data_dir, "gasp_scale_gRNA_counts.odm")
metadata_fp <- paste0(processed_data_dir, "gasp_scale_gRNA_metadata.rds")
gRNA_count_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))
cell_barcodes <- colnames(gRNA_count_matrix)
feature_barcodes <- data.frame(row.names(gRNA_count_matrix))
gRNA_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_count_matrix,
                                               barcodes = cell_barcodes,
                                               features_df = feature_barcodes,
                                               odm_fp = odm_fp,
                                               metadata_fp = metadata_fp)
