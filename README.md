#  scifi-scRNA-seq processing pipeline
This is a pipeline to process scifi-scRNA-seq ([Datlinger et al.](https://doi.org/10.1038/s41592-021-01153-z)) data from fastq files. It utilizes a custom python script to preprocess the fastq files, splitting them by well/sample and formatting them for downstream processing. We use [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) to quantify exons and introns from [salmon](https://github.com/COMBINE-lab/salmon) output, largely based on [this](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) tutorial.

## Dependencies
This pipeline requires the following exectuables in your PATH:
* [Snakemake](https://github.com/snakemake/snakemake) (developed with version 3.11.2)
* [Python 3](https://www.python.org/) (developed with version 3.6.3)
* [R](https://www.r-project.org/) (developed with version 4.0.3)
* [salmon](https://github.com/COMBINE-lab/salmon) (developed with version 1.5.2)
* [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) (developed with version 0.4.1)

## Usage
To run this pipeline, edit `config.yaml` to specify the following parameters:

1. `RAW_DATA`: This is the directory that contains the raw fastq data and experiment allowlist (a 3-column, tab-separated file that contains the well id, sample name, and sample index of each sample).
2. `OUT_DATA`: This is the directory that the output results will be placed.
3. `FASTA`: This is the path to a fasta reference file, appropriate to the experiment.
4. `GTF`: This is the path to a gene annotation file, appropriate to the experiment.
5. `SAMPLES`: This is the path to the allowlist contained in the `RAW_DATA` directory. Sample names consist of the well id followed by the sample name, i.e. `{well_id}_{sample_name}`.
6. `SRC`: This is the `src` directory that contains the `preprocessing.py` and `make_splici_txome.R` scripts.

## Output Files
If successful you should generate the following files per sample (located under `{OUT_DIR}/scifi/{sample}/quant_res/alevin`):

* quants_mat.mtx (the count matrix; [num_cells] x 3[num_genes])
* quants_mat_cols.txt (the gene names for each column of the matrix)
* quants_mat_rows.txt (the cell barcodes for each row of the matrix)

The counts matrix contains UMIs separately attributed to _spliced_ or _unspliced_ (intronic) gene sequences, or as _ambiguous_. For more information, see the [alevin-fry tutorial](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/).

## Downstream Processing
After running this pipeline subsequent processing can be done in R or Python using the `load_fry.py` or `load_fry.R` scripts available in `/bin/downstream_processing`.
<br />
<br />
<br />
See the [Snakemake](https://snakemake.readthedocs.io/en/stable/) documentation for general information in executing and manipulating snakemake pipelines.

https://github.com/rnabioco/scifi