# qc_check.R checks the quality of scifi data using `emptyDrops` (`DropletUtils`)
# Input: path to dir containing *counts.tsv files from running `umi-tools counts`
# TODO: allow for argument passing (get rid of hard-coded example)
# Returns summary .txt file for each scifi sample
# TODO: merge summary files into one big file?
# Returns barcodeRanks plot of each scifi sample
# Note: in this instance, one scifi sample corresponds to one sample well from a 96-well plate

# load in relevant libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DropletUtils))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))

setwd("~/Desktop/santoro/scifi/JH207/emptyDrops_out/") # set home dir
dir.create("summaries") # output dir for summary .txt files
dir.create("plots") # output dir for barcodeRanks plots

# get list of sample *counts.tsv files; ignore files from blank wells
# TODO: should we include files from blank wells for qc comparison?
files <- list.files(path = "~/Desktop/santoro/scifi/JH207/sample_counts", pattern = "*.tsv", full.names = TRUE, recursive = FALSE)
files <- files[!grepl("blank", files)]

qc_check <- function(file){
  res = tryCatch({ # `emptyDrops` fails on some files; TODO: move tryCatch?
    sample_name <- strsplit(file, "/")[[1]][9]
    sample_name <- gsub('_counts.tsv', '', sample_name)
    output <- file(str_interp('summaries/${sample_name}_summary.txt'))
    
    # read in counts.tsv file as data.frame df and convert to sparse matrix smtx
    suppressMessages(df <- read_tsv(file))
    dfw <- pivot_wider(df, names_from = "cell", values_from = "count", values_fill = 0L)
    dfw <- column_to_rownames(dfw, "gene")
    mtx <- as.matrix(dfw)
    smtx <- as(mtx, "dgCMatrix")
      
    ln1 <- str_interp('Sample: ${sample_name}')
    ln2 <- str_interp('Total # of entries: ${nrow(df)}')
    ln3 <- str_interp('# of UMIs: ${sum(df$count)}')
    ln4 <- str_interp('# of barcodes: ${n_distinct(df$cell)}')
    ln5 <- str_interp('# of genes: ${n_distinct(df$gene)}')
  
    set.seed(100)
    e.out <- emptyDrops(smtx)
          
    is.cell <- e.out$FDR <= 0.01
    num_cells <- sum(is.cell, na.rm=TRUE)
    ln6 <- str_interp('# of cells captured: ${num_cells}')
    
    # create subset e.out data.frame to contain cells found by `emptyDrops`
    e.out.sub <- subset(e.out, FDR <= 0.01)
          
    total_umi <- sum(e.out$Total)
    sub_umi <- sum(e.out.sub$Total)
    per_umi <- (sub_umi/total_umi)*100
    ln7 <- str_interp('% of UMIs in cell-associated droplets: ${format(round(per_umi, 3), nsmall = 3)}%')
          
    avg_umi <- sub_umi/num_cells
    ln8 <- str_interp('average # of UMIs in cell-associated droplets: ${format(round(avg_umi, 3), nsmall = 3)}')
    
    # subset original data.frame df to contain cells found by `emptyDrops`     
    cell_ids <- row.names(e.out.sub)
    df_sub <- df %>% filter(cell %in% cell_ids)
    gene_count <- count(df_sub, "cell")
    total_genes <- sum(gene_count$freq)
    avg_genes <- total_genes/num_cells
    ln9 <- str_interp('average # of genes detected per cell-associated droplet: ${format(round(avg_genes, 3), nsmall = 3)}')
          
    writeLines(c(ln1, ln2, ln3, ln4, ln5, ln6, ln7, ln8, ln9), output)
    close(output)
          
    br.out <- barcodeRanks(smtx)
          
    # make the `barcodeRanks` plot (from `DropletUtils` tutorial)
    plot_name <- sub("%", "%%", str_interp('${sample_name}_barcodeRanks.png'))
    png(paste("plots/", plot_name))
    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col="red")
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))
    dev.off()
    
  }, error = function(e){}) # TODO: add more informative error function
  res
}

lapply(files, qc_check)