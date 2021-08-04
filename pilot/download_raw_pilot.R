# source config file to get gasperini offsite location
gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
# load R.utils; increase timeout to 5 hours
library(R.utils)
options(timeout = 5 * 60 * 60)

# create raw directory
raw_data_dir_gasp <- paste0(gasp_offsite, "raw")
if (!dir.exists(raw_data_dir_gasp)) dir.create(path = raw_data_dir_gasp, recursive = TRUE)

# URL of data
remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file="

# names of genes -- ordered gene IDs for use in conjunction with expression mtx file
genes_filename <- "GSE120861_pilot_highmoi_screen.genes.txt"

# names of cells -- cell barcodes for use in conjunction with expression mtx file
cells_filename <- "GSE120861_pilot_highmoi_screen.cells.txt"

# all (gRNA, gene) pairs 
gRNAgroup_pair_table_filename <- "GSE120861_gene_gRNAgroup_pair_table.pilot.txt"

# list of gRNA groups used
# gRNA_groups_filename <- "GSE120861_grna_groups.pilot.txt"

# monocle Cell Data Set object with binary gRNA data
cds_filename <- "GSE120861_pilot_highmoi_screen.cds.rds"

# gene expression matrix in mtx format
expression_filename <- "GSE120861_pilot_highmoi_screen.exprs.mtx"

# cell phenotype data 
cell_phenodata_filename <- "GSE120861_pilot_highmoi_screen.phenoData.txt"

# list of files to download
filenames <- c(genes_filename,
               cells_filename,
               cds_filename,
               expression_filename,
               gRNAgroup_pair_table_filename,
               cell_phenodata_filename)

# download files if not already present from GEO
for (filename in filenames) {
  if (!file.exists(paste0(raw_data_dir_gasp, "/", filename))) {
    print(paste0("Downloading ", filename))
    source <- paste0(remote, filename, ".gz")
    dest <- paste0(raw_data_dir_gasp, "/", filename, ".gz")
    download.file(source, dest)
    gunzip(paste0(dest))
  }
}


