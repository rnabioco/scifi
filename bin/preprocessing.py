# preprocessing.py
#
# Take the I1, R1, R2, and R3 fastq files from a scifi experiment and do the following:
#   - split files by well-specific sample index
#   - format for Alevin (10X v3 data)
#
# New R1 files will contain the cell barcodes and UMIs:
#   - cell barcode: 16 base index sequence from the 10X ATAC-seq oligos (from original R2)
#   - UMI: 8 base sequence ([0:8] from original R1)
# New R2 files will contain the mRNA sequence (from original R3)
#
# Expected fastq file input format:
#   - R1: 8 base UMI + 1 fixed base + 11 base round 1 barcode + 1 fixed base
#   - R2: 16 base round 2 barcode
#   - R3: 70 base mRNA sequence
#   - I1: 8 base i7 sample index
#
# TODO: do we need to include the sample index in individual file headers?
# TODO: would is be more useful to label files with the well ID vs sample index (i.e., A4 vs TGTAAGAATCA)?

import sys
import argparse
import os
import glob
from itertools import islice
import gzip

# Specifying command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("data_dir", help = "path to folder with fastq files from scifi scRNA-seq experiment (I1, R1, R2, R3)")
parser.add_argument("allowlist", help = "path to allowlist")

args = parser.parse_args()

# Helper functions
#   -getsamples; get well-specific sample indexes from provided allowlist
#   -getchunk; go through fastq files 4 lines at a time
#   -writechunk; store chunks and corresponding filenames in dictionary for efficient file writing

def getsamples(file):
    with open(file) as f:
        samples = f.readlines()
    samples = [sample.strip() for sample in samples]
    return samples

def getchunk(file): 
    return [x.strip() for x in islice(file, 4)]

def writechunk(file, chunk, file_dict): 
    if file not in file_dict.keys():
        file_dict[file] = gzip.open(file, "wt")
    file_dict[file].write("\n".join(chunk) + "\n")

# Main function
def main(argv):
    # get list of well-specific sample indexes from allowlist
    samples = getsamples(args.allowlist)

    # change dir to data_dir
    os.chdir(args.data_dir)

    # make outfile directory, "preprocessed_fastqs"
    outdir = "preprocessed_fastqs"
    os.makedirs(f"./{outdir}", exist_ok = True)

    # initialize output dictionary
    output_dict = {}

    try:
        # open I1, R1, R2, and R3 fastq read files
        i1 = gzip.open(glob.glob("*I1*.fastq.gz")[0], "rt")
        r1 = gzip.open(glob.glob("*R1*.fastq.gz")[0], "rt")
        r2 = gzip.open(glob.glob("*R2*.fastq.gz")[0], "rt")
        r3 = gzip.open(glob.glob("*R3*.fastq.gz")[0], "rt")

        while True:
            # get reads in 4 line chunks from corresponding infiles
            r1_chunk = getchunk(r1)
            if not r1_chunk:
                break
            r2_chunk = getchunk(r2)
            r3_chunk = getchunk(r3)
            i1_chunk = getchunk(i1)

            # get well-specific sample index from original R1
            sample_index = f"{r1_chunk[1][9:20]}"

            # get cell barcdoe from original R2 and UMI from original R1
            cell_barcode = f"{r2_chunk[1]}"
            cell_barcode_qc = f"{r2_chunk[3]}"
            umi = f"{r1_chunk[1][0:8]}"
            umi_qc = f"{r1_chunk[3][0:8]}"

            # write new R1 sequence info (cell barcode + UMI)
            r1_chunk[1] = f"{cell_barcode}{umi}"
            r1_chunk[3] = f"{cell_barcode_qc}{umi_qc}"

            # if sample_index is in the allowlist of well-specific sample indexes:
            # TODO: is there a way to allow for sequencing error (off by one)? Levenshtein? Custom helper function?
            if sample_index in samples:
                writechunk(f"preprocessed_fastqs/{sample_index}_I1.fastq.gz", i1_chunk, output_dict)
                writechunk(f"preprocessed_fastqs/{sample_index}_R1.fastq.gz", r1_chunk, output_dict)
                writechunk(f"preprocessed_fastqs/{sample_index}_R2.fastq.gz", r3_chunk, output_dict)

            # if sample_index is not in the allowlist of well-specific sample indexes it is an orphaned read:
            else:
                writechunk(f"preprocessed_fastqs/orphan_I1.fastq.gz", i1_chunk, output_dict)
                writechunk(f"preprocessed_fastqs/orphan_R1.fastq.gz", r1_chunk, output_dict)
                writechunk(f"preprocessed_fastqs/orphan_R2.fastq.gz", r3_chunk, output_dict)
            
    except IOError:
        print("Could not open file.")
        sys.exit(1)

    except Exception as e:
        print(e)
        sys.exit(1)

    finally: # close the files
        i1.close()
        r1.close()
        r2.close()
        r3.close()

        for _, file in output_dict.items():
            file.close()

if __name__ == '__main__':
    main(sys.argv[1:])