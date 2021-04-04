# make_big_mat.R is the same as make_mat.R, but combines all scifi samples into one output (to be used as a qc check, not for experimental analysis)
# Input: path to dir containing *counts.tsv files from running `umi-tools counts`
# Input: path to .gtf file used to run `STAR`
# TODO: allow for argument passing (get rid of hard-coded example)
# Returns combined matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz files for all scifi samples (CellRanger version 3.0 format)
# Note: in this instance, one scifi sample corresponds to one sample well from a 96-well plate

# load in relevant libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DropletUtils))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(BUSpaRse))

setwd("~/Desktop/santoro/scifi/JH207/emptyDrops_out/")
dir.create("cellranger") # output dir

# get list of sample *counts.tsv files; ignore files from blank wells
files <- list.files(path = "~/Desktop/santoro/scifi/JH207/sample_counts", pattern = "*.tsv", full.names = TRUE, recursive = FALSE)
files <- files[!grepl("blank", files)]

# read in .gtf file and grab gene transcript, id, and symbol info
gtf <- tr2g_gtf(file = "/Users/Caitlin/Desktop/santoro/scifi/JH207/data/genes.gtf",
                gene_biotype_use = "cellranger", get_transcriptome = FALSE)
# drop decimals from gene ids since those don't appear in our data
# TODO: will this always be the case, or do we need to add options here?
gtf$gene_adj <- gsub("[.][0-9]+$", "", gtf$gene)

make_big_mat <- function(file){
  res = tryCatch({ # `emptyDrops` fails on some files; TODO: move tryCatch?
    sample_name <- strsplit(file, "/")[[1]][9]
    sample_name <- gsub('_counts.tsv', '', sample_name)
    
    # read in counts.tsv file as data.frame df and convert to sparse matrix smtx
    suppressMessages(df <- read_tsv(file))
    dfw <- pivot_wider(df, names_from = "cell", values_from = "count", values_fill = 0L)
    dfw <- column_to_rownames(dfw, "gene")
    mtx <- as.matrix(dfw)
    smtx <- as(mtx, "dgCMatrix")
    
    set.seed(100)
    e.out <- emptyDrops(smtx)
    
    # create subset e.out data.frame to contain cells found by `emptyDrops`
    e.out.sub <- subset(e.out, FDR <= 0.01)
    
    # subset original data.frame df to contain cells found by `emptyDrops`
    cell_ids <- row.names(e.out.sub)
    df_sub <- df %>% filter(cell %in% cell_ids)
    
    # return subsetted data.frame df_sub
    res <- df_sub
    res
  }, error = function(e){})
  res
}

dfs <- lapply(files, make_big_mat)

# combine data.frames from individual samples into one data.frame big_df
big_df <- do.call(rbind, dfs)

# convert data.frame big_df to sparse matrix big_smtx
big_df <- pivot_wider(big_df, names_from = "cell", values_from = "count", values_fill = 0L)
big_df <- column_to_rownames(big_df, "gene")
big_mtx <- as.matrix(big_df)
big_smtx <- as(big_mtx, "dgCMatrix")

# save gene and cell barcodes
cell.ids <- colnames(big_smtx)
gene.ids <- rownames(big_smtx)

# get corresponding gene symbols for gene ids
gene.symbs <- vector()
for(gene in gene.ids){
  id <- match(gene, gtf$gene_adj)
  gene.symbs <- c(gene.symbs, gtf$gene_name[id])
}

# use `DropletUtils` function `write10xCounts` to output data in CellRanger format
tmpdir <- paste0("/Users/Caitlin/Desktop/santoro/scifi/JH207/emptyDrops_out/cellranger/combined")
write10xCounts(tmpdir, 
               big_smtx, 
               gene.id = gene.ids, 
               gene.symbol = gene.symbs, 
               barcodes = cell.ids,
               version = '3')