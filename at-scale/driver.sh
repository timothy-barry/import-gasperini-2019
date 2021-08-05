# download raw data
Rscript download_raw.R
Rscript  process_data.R
Rscript create_odms.R
Rscript obtain_gRNA_gene_pairs.R

# note that verify_gRNA_matrix.R is not actually part of the pipeline; this script ensures that we have handled the gRNA data correctly.