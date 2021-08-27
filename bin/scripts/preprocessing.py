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

import sys
import argparse
import os
import glob
from pysam import FastxFile
from itertools import islice
import gzip

# Specifying command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("data_dir", help = "path to folder with fastq files from scifi scRNA-seq experiment (I1, R1, R2, R3)")
parser.add_argument("allowlist", help = "path to allowlist (a tab-separated file that contains the well id, sample name, and sample index of each sample")
parser.add_argument("tolerance", help = "hamming distance for sequencing error allowance", nargs = "?", default = 1)

args = parser.parse_args()

# Helper functions
#   -get_samples; get well-specific ids, sample names, and sample indexes from provided allowlist and store in dictionary
#   -write_entry; store entries and corresponding filenames in dictionary for efficient file writing
#   -hamming_distance; calculates the number of bases in which two sequences of equal length are different

def get_samples(file):
    sample_dict = {}
    with open(file) as f:
        for line in f:
            splitline = line.split()
            well_id = splitline[0]
            sample_name = splitline[1]
            sample_seq = splitline[2]
            sample_dict[sample_seq] = (well_id, sample_name)
    return sample_dict

def write_entry(file, entry, file_dict): 
    if file not in file_dict.keys():
        file_dict[file] = gzip.open(file, "wt", compresslevel = 6)
    file_dict[file].write(str(entry) + "\n")

def hamming_distance(seq1, seq2):
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))

# Main function
def main(argv):
    # get list of well-specific sample indexes from allowlist
    sample_dict = get_samples(args.allowlist)
    samples = list(sample_dict.keys())

    # change dir to data_dir
    os.chdir(args.data_dir)

    # make outfile directory, "preprocessed_fastqs"
    outdir = "preprocessed_fastqs"
    os.makedirs(f"../{outdir}", exist_ok = True)

    # initialize output dictionary
    output_dict = {}

    try:
        # open I1, R1, R2, and R3 fastq read files
        i1 = FastxFile(glob.glob("*I1*.fastq.gz")[0])
        r1 = FastxFile(glob.glob("*R1*.fastq.gz")[0])
        r2 = FastxFile(glob.glob("*R2*.fastq.gz")[0])
        r3 = FastxFile(glob.glob("*R3*.fastq.gz")[0])

        while True:
            # get reads in 4 line entries from corresponding infiles
            r1_entry = next(r1)
            if not r1_entry:
                break
            r2_entry = next(r2)
            r3_entry = next(r3)
            i1_entry = next(i1)

            # get well-specific sample index from original R1
            sample_index = r1_entry.sequence[9:20]

            # get cell barcdoe from original R2 and UMI from original R1
            cell_barcode = r2_entry.sequence
            cell_barcode_qc = r2_entry.quality
            umi = r1_entry.sequence[0:8]
            umi_qc = r1_entry.quality[0:8]

            # write new R1 sequence info (cell barcode + UMI)
            r1_entry.sequence = cell_barcode + umi
            r1_entry.quality = cell_barcode_qc + umi_qc

            # if sample_index is in the allowlist of well-specific sample indexes (hamming distance <= 1 by default):
            match = next((sample for sample in samples if hamming_distance(sample, sample_index) <= args.tolerance), None)
            if match:
                write_entry(f"../preprocessed_fastqs/{sample_dict[match][0]}_{sample_dict[match][1]}_I1.fastq.gz", i1_entry, output_dict)
                write_entry(f"../preprocessed_fastqs/{sample_dict[match][0]}_{sample_dict[match][1]}_R1.fastq.gz", r1_entry, output_dict)
                write_entry(f"../preprocessed_fastqs/{sample_dict[match][0]}_{sample_dict[match][1]}_R2.fastq.gz", r3_entry, output_dict)

            # if sample_index is not in the allowlist of well-specific sample indexes it is an unassigned read:
            else:
                write_entry(f"../preprocessed_fastqs/unassigned_I1.fastq.gz", i1_entry, output_dict)
                write_entry(f"../preprocessed_fastqs/unassigned_R1.fastq.gz", r1_entry, output_dict)
                write_entry(f"../preprocessed_fastqs/unassigned_R2.fastq.gz", r3_entry, output_dict)
            
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