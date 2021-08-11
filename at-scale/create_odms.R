###########################
# 0. Load packages; set fps
###########################
# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/")

# create the processed data directories
processed_data_dir <- paste0(gasp_offsite, "processed/")
processed_data_dir_gene <- paste0(processed_data_dir, "gene/")
processed_data_dir_grouped <- paste0(processed_data_dir, "gRNA_grouped/")
processed_data_dir_ungrouped <- paste0(processed_data_dir, "gRNA_ungrouped/")
dirs_to_create <- c(processed_data_dir_gene,
                    processed_data_dir_grouped,
                    processed_data_dir_ungrouped)
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) dir.create(path = dir, recursive = TRUE)
} 

# set raw directories
raw_data_dir <- paste0(gasp_offsite, "raw/")
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")

# load packages
library(magrittr)
library(ondisc)

###########################
# 1. gene expression matrix
###########################
# gene count matrix
mtx_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.exprs.mtx")
barcodes_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.cells.txt")
gene_ids_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.genes.txt")

odm_fp <- paste0(processed_data_dir_gene, "gasp_scale_gene_expressions.odm")
metadata_fp <- paste0(processed_data_dir_gene, "gasp_scale_gene_metadata.rds")

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

##################################
# 2. gRNA count matrix (ungrouped)
##################################
# next, load the gRNA count matrix. Write the backing .odm file and metadata file to disk.
odm_fp <- paste0(processed_data_dir_ungrouped, "gasp_scale_gRNA_counts_ungrouped.odm")
metadata_fp <- paste0(processed_data_dir_ungrouped, "gasp_scale_gRNA_metadata_ungrouped.rds")
gRNA_count_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))
gRNA_groups_df <- readr::read_tsv(file = paste0(raw_data_dir,
                                             "GSE120861_grna_groups.at_scale.txt"),
                               col_types = "cc", col_names = c("gRNA_group", "barcodes"))
# confirm that ordering of gRNA groups matches that of barcodes
identical(row.names(gRNA_count_matrix), gRNA_groups_df$barcodes)
cell_barcodes <- colnames(gRNA_count_matrix)
# create the gRNA odm
gRNA_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_count_matrix,
                                               barcodes = cell_barcodes,
                                               features_df = gRNA_groups_df[,c(2,1)],
                                               odm_fp = odm_fp,
                                               metadata_fp = metadata_fp)

################################
# 3. gRNA count matrix (grouped)
################################
tab <- table(gRNA_groups_df$gRNA_group)
single_gRNA_groups <- names(tab)[which(tab == 1)]
double_gRNA_groups <- names(tab)[which(tab == 2)]

# extract the first and second sets of barcodes for the doubles
doubles_grouped <- dplyr::filter(gRNA_groups_df, gRNA_group %in% double_gRNA_groups) %>%
  dplyr::group_by(gRNA_group) %>% dplyr::arrange(gRNA_group) 
doubles_barcode_1 <- doubles_grouped %>% dplyr::summarize(barcode = barcodes[1]) %>%
  dplyr::pull(barcode)
doubles_barcode_2 <- doubles_grouped %>% dplyr::summarize(barcode = barcodes[2]) %>%
  dplyr::arrange(gRNA_group) %>% dplyr::pull(barcode)
ordered_groups <- unique(doubles_grouped$gRNA_group)
# verify on a few of the groups that the extraction was correct
for (i in sample(seq(1, length(ordered_groups)), 15, FALSE)) {
  grp <- ordered_groups[i]
  bc_check <- c(doubles_barcode_1[i], doubles_barcode_2[i])
  bc_test <- dplyr::filter(doubles_grouped, gRNA_group == grp) %>% dplyr::pull(barcodes)
  testthat::expect_true((bc_check %in% bc_test) && (bc_test %in% bc_check))
}
# sum the count vectors
double_m1 <- gRNA_count_matrix[doubles_barcode_1,]
double_m2 <- gRNA_count_matrix[doubles_barcode_2,]
double_sum <- double_m1 + double_m2
row.names(double_sum) <- ordered_groups
# append "double_sum" to the matrix of counts of singlet groups
singlet_map <- dplyr::filter(gRNA_groups_df, gRNA_group %in% single_gRNA_groups)
singlet_m <- gRNA_count_matrix[singlet_map$barcodes,]
row.names(singlet_m) <- singlet_map$gRNA_group
# create the final, "combined" matrix
combined_matrix <- rbind(double_sum, singlet_m)
combined_matrix <- as(combined_matrix, "dgTMatrix")
# set the file paths
odm_fp <- paste0(processed_data_dir_grouped, "gasp_scale_gRNA_counts_grouped.odm")
metadata_fp <- paste0(processed_data_dir_grouped, "gasp_scale_gRNA_metadata_grouped.rds")
gRNA_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = combined_matrix,
                                               barcodes = colnames(combined_matrix),
                                               features_df = data.frame(gRNA_id = row.names(combined_matrix)),
                                               odm_fp = odm_fp,
                                               metadata_fp = metadata_fp)
