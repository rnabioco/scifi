# qc_check_df.R is similar to qc_check.R, but stores qc metrics in a data.frame
# Input: path to dir containing *counts.tsv files from running `umi-tools counts`
# TODO: allow for argument passing (get rid of hard-coded example)
# Returns combined qc_df data.frame for all scifi samples
# Note: in this instance, one scifi sample corresponds to one sample well from a 96-well plate

# load in relevant libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DropletUtils))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))

setwd("~/Desktop/santoro/scifi/JH207/emptyDrops_out/")


# get list of sample *counts.tsv files; ignore files from blank wells
# TODO: should we include files from blank wells for qc comparison?
files <- list.files(path = "~/Desktop/santoro/scifi/JH207/sample_counts", pattern = "*.tsv", full.names = TRUE, recursive = FALSE)
files <- files[!grepl("blank", files)]

qc_check <- function(file){
  res = tryCatch({ # `emptyDrops` fails on some files; TODO: move tryCatch?
    sample_name <- strsplit(file, "/")[[1]][9]
    sample_name <- gsub('_counts.tsv', '', sample_name)
    
    # read in counts.tsv file as data.frame df and convert to sparse matrix smtx
    suppressMessages(df <- read_tsv(file))
    dfw <- pivot_wider(df, names_from = "cell", values_from = "count", values_fill = 0L)
    dfw <- column_to_rownames(dfw, "gene")
    mtx <- as.matrix(dfw)
    smtx <- as(mtx, "dgCMatrix")
    
    val1 <- sample_name # Sample name
    val2 <- nrow(df) # Total # of entries
    val3 <- sum(df$count) # # of UMIs
    val4 <- n_distinct(df$cell) # # of barcodes
    val5 <- n_distinct(df$gene) # # of genes
    
    set.seed(100)
    e.out <- emptyDrops(smtx)
    
    is.cell <- e.out$FDR <= 0.01
    num_cells <- sum(is.cell, na.rm=TRUE)
    val6 <- num_cells # # of cells captured
    
    # create subset e.out data.frame to contain cells found by `emptyDrops`
    e.out.sub <- subset(e.out, FDR <= 0.01)
    
    total_umi <- sum(e.out$Total)
    sub_umi <- sum(e.out.sub$Total)
    per_umi <- (sub_umi/total_umi)*100
    val7 <- per_umi # % of UMIs in cell-associated droplet
    
    avg_umi <- sub_umi/num_cells
    val8 <- avg_umi # avg # of UMIs in cell-associated droplets
    
    # subset original data.frame df to contain cells found by `emptyDrops`
    cell_ids <- row.names(e.out.sub)
    df_sub <- df %>% filter(cell %in% cell_ids)
    gene_count <- count(df_sub, "cell")
    total_genes <- sum(gene_count$freq)
    avg_genes <- total_genes/num_cells
    val9 <- avg_genes # avg # of genes in cell-associated droplets
    
    # store qc metrics in data.frame
    res <- data.frame('Sample' = val1,
                      'Total # of entries' = val2,
                      '# of UMIs' = val3,
                      '# of barcodes' = val4,
                      '# of genes' = val5,
                      '# of cells captured' = val6,
                      '% of UMIs in cell-associated droplets'  = val7,
                      'avg # of UMIs in cell-associated droplets' = val8,
                      'avg # of genes in cell-associated droplets' = val9)
    res
  }, error = function(e){})
  res
}

dfs <- lapply(files, qc_check)

# combine data.frames from individual samples into one data.frame qc_df
qc_df <- do.call(rbind, dfs)

save(qc_df, file = "summaries/qc_df.Rda")