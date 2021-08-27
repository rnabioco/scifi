#  scifi-scRNA-seq processing pipeline
This is a pipeline to process scifi-scRNA-seq ([Datlinger et al.](https://doi.org/10.1038/s41592-021-01153-z)) data from fastq files. It utilizes a custom python script to preprocess the fastq files, splitting them by well/sample and formatting them for downstream processing. We use [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) to quantify exons and introns from [salmon](https://github.com/COMBINE-lab/salmon) output, largely based on [this](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) tutorial.

## Dependencies
This pipeline requires the following exectuables in your PATH:
* [Snakemake](https://github.com/snakemake/snakemake) (developed with version 6.7.0)
* [Python 3](https://www.python.org/) (developed with version 3.6.3)
* [R](https://www.r-project.org/) (developed with version 4.0.3)
* [salmon](https://github.com/COMBINE-lab/salmon) (developed with version 1.5.2)
* [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) (developed with version 0.4.1)

## Example Usage
This pipeline requires the following as input (defined in config.yaml):
* scifi-RNA-seq FASTQs
* allowlist (tab-separated file) containing corresponding well ids, sample names, and sample index sequences
* 10X Genomics reference directory

See the [Snakemake](https://snakemake.readthedocs.io/en/stable/) documentation for general information in executing and manipulating snakemake pipelines.

https://github.com/rnabioco/scifi