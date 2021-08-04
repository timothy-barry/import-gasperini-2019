### step 2 generate intermediate data ###
# source config file to get gasperini offsite location
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
# gasp_offsite <- c("/Users/sijia_work/Documents/gasperini-2019/")

library(readr)
# create the intermediate data directory; set raw directory
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")
if (!dir.exists(intermediate_data_dir)) dir.create(path = intermediate_data_dir, recursive = TRUE)
raw_data_dir <- paste0(gasp_offsite, "raw/")
processed_data_dir <- paste0(gasp_offsite, "processed/")
if (!dir.exists(processed_data_dir)) dir.create(path = processed_data_dir, recursive = TRUE)


library(ondisc)
# Obtain binary gRNA matrix and cell metadata from monocole object
library(monocle)
library(magrittr)
monocle_obj <- readRDS(paste0(raw_data_dir, "/GSE120861_pilot_highmoi_screen.cds.rds"))
cell_metadata <- pData(monocle_obj)
rm(monocle_obj); gc()
covariates_cols <- 1:14
gRNA_cols <- 15:ncol(cell_metadata)

gRNA_indicators <- cell_metadata[,gRNA_cols]
cell_covariates <- cell_metadata[,covariates_cols]
gRNA_names <- colnames(cell_metadata[,gRNA_cols])

### added from vignette ###
gRNA_features <- data.frame(gRNA_id = gRNA_names)
write_tsv(x = gRNA_features, file = paste0(raw_data_dir, "/gRNAs.tsv"), col_names = FALSE)
gRNA_indics <- as.matrix(gRNA_indicators)
colnames(gRNA_indics) <- row.names(gRNA_indics) <- NULL
sparse_gRNA_indics <- t(Matrix(gRNA_indics, sparse = TRUE))
writeMM(obj = sparse_gRNA_indics, file = paste0(raw_data_dir, "/perturbations.mtx"))

features_fp <- paste0(raw_data_dir, "/gRNAs.tsv")
cell_fp <- paste0(raw_data_dir, "/GSE120861_pilot_highmoi_screen.cells.txt")
mtx_fp <- paste0(raw_data_dir, "/perturbations.mtx")

cell_barcodes <- read_tsv(paste0(raw_data_dir, "/GSE120861_pilot_highmoi_screen.cells.txt"),col_names = FALSE)

library(dplyr)
perturbations <- create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_indics,
                                               barcodes = as.character(dplyr::pull(cell_barcodes)),
                                               features_df = as.data.frame(as.character(dplyr::pull(gRNA_features)),stringsAsFactors=FALSE),
                                               on_disk_dir = processed_data_dir,
                                               return_metadata_ondisc_matrix = TRUE)
saveRDS(object = perturbations, file = paste0(raw_data_dir, "/perturbations.rds"))

saveRDS(cell_covariates, paste0(intermediate_data_dir, "cell_covariates.rds"))
saveRDS(gRNA_indicators, paste0(intermediate_data_dir, "gRNA_indicators.rds"))

