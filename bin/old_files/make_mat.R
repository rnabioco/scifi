# make_mat.R uses `DropletUtils` functions to return a directory containing a count matrix and cell/gene annotation from a spare matrix of UMI counts, in the format produced by CellRanger 
# This allows scifi data to be easily processed using popular downstream analysis tools such as Seurat and scanpy
# Input: path to dir containing *counts.tsv files from running `umi-tools counts`
# Input: path to .gtf file used to run `STAR`
# TODO: allow for argument passing (get rid of hard-coded example)
# Returns matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz files for each scifi sample (CellRanger version 3.0 format)
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

make_mat <- function(file){
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
    df_sub_w <- pivot_wider(df_sub, names_from = "cell", values_from = "count", values_fill = 0L)
    df_sub_w <- column_to_rownames(df_sub_w, "gene")
    mtx_sub <- as.matrix(df_sub_w)
    smtx_sub <- as(mtx_sub, "dgCMatrix")
    
    # save gene and cell barcodes
    cell.ids <- colnames(smtx_sub)
    gene.ids <- rownames(smtx_sub)
    
    # get corresponding gene symbols for gene ids
    gene.symbs <- vector()
    for(gene in gene.ids){
      id <- match(gene, gtf$gene_adj)
      gene.symbs <- c(gene.symbs, gtf$gene_name[id])
    }
    
    # use `DropletUtils` function `write10xCounts` to output data in CellRanger format
    tmpdir <- paste0("/Users/Caitlin/Desktop/santoro/scifi/JH207/emptyDrops_out/cellranger/", sample_name)
    write10xCounts(tmpdir, 
                   smtx_sub, 
                   gene.id = gene.ids, 
                   gene.symbol = gene.symbs, 
                   barcodes = cell.ids,
                   version = '3')
    
  }, error = function(e){})
  res
}

lapply(files, make_mat)