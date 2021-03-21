# fastq_merge.py
# Take the R1, R2, and R3 fastq files from a scifi experiment and convert them to one merged fastq file (merged.fastq)
# merged.fastq contains the mRNA reads with the cell barcode + UMI encoded in the header
# cell barcodes = 16bp index sequence from the 10X ATAC-seq oligos (from R2) + 11bp index sequence from the RT primers ([9:20] from R1)
# UMI = 8bp sequence ([0:8] from R1)

# Example header
# Before running fastq_merge.py: @A00405:365:H72K3DRXY:1:2101:1452:1063 3:N:0:CGACTTAG
# After running fastq_merge.py: @A00405:365:H72K3DRXY:1:2101:1452:1063 3:N:0:CGACTTAG_AACGTTAGCTTACAGCCCTGCCCCTGT_CNTTGAAT
# Format: @A00405:365:H72K3DRXY:1:2101:1452:1063 3:N:0:CGACTTAG_CELLBARCODE_UMI
# CELLBARCODE_UMI = (16bp_index_seq)(11bp_index_seq)_(8bp_seq)

import sys
import argparse
import os
import glob
from itertools import islice

# Specifying command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("data_dir", help="path to folder with fastq files from scifi scRNA-seq experiment (R1, R2, R3)")

args = parser.parse_args()

def getchunk(file): # helper function to go through fastq files 4 lines at a time
    return [x.strip() for x in islice(file, 4)]

def main(argv):
    # change dir to data_dir
    os.chdir(args.data_dir)

    try:
        # open R1, R2, and R3 fastq read files
        r1 = open(glob.glob("*R1*.fastq")[0])
        r2 = open(glob.glob("*R2*.fastq")[0])
        r3 = open(glob.glob("*R3*.fastq")[0])

        # open outfile as {sample}_merged.fastq
        sample = glob.glob("*R1*.fastq")[0].split("_")[0]
        outfile = open(f"{sample}_merged.fastq", "w")

        while True:
            # get mRNA reads
            r3_chunk = getchunk(r3)
            if not r3_chunk:
                break
            # add 16bp index sequence (from 10X ATAC-seq oligos) to header (from R2)
            r3_chunk[0] += f"_{getchunk(r2)[1]}"
            
            # add 11bp index sequence (from RT primers) and 8bp UMI to header (from R1)
            r1_chunk = getchunk(r1)
            r3_chunk[0] += f"{r1_chunk[1][9:20]}_{r1_chunk[1][0:8]}"

            # write mRNA reads and updated headers to outfile
            for line in r3_chunk:
                outfile.write(f"{line}\n")
            
    except IOError:
        print("Could not open file.")
        sys.exit(1)

    except:
        print("Error.")
        sys.exit(1)

    finally: # close the files
        r1.close()
        r2.close()
        r3.close()
        outfile.close()

if __name__ == '__main__':
    main(sys.argv[1:])