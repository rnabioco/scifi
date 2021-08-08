# split_counts.py
# Splits a given *_counts.tsv.gz file (output from umi_tools count) into separate samples as defined by an allowlist.tsv file
# allowlist.tsv is expected to contain three tab-deliminated columns: sample_index, sample_name, and plate_well
# sample_index is the 11bp index sequence from the RT primers, which is unqiue per plate well; sample_index is part of the cell barcode
# sample_name is the biologically informative name of samples distributed across a 96 well plate
# plate_well is the positional identifier of a well that corresponds to a specific sample_index

# Example allowlist format:
# sample_index    sample_name plate_well
# AGTGATTAGCA	P14-M1-L-1%	A01
# GGCCGGAACAT	P28-F2-L-4%	A11
# CCCTAAACGCT	P28-M1-L-1%	E05
# CGTACGCTTAT	blank	G09

import sys
import argparse
import os
import glob
import gzip

# Specifying command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("data_dir", help="path to directory containing umi_tools count output file and sample index allowlist")

args = parser.parse_args()

def main(argv):
    # change dir to data_dir
    os.chdir(args.data_dir)

    # grab count file and sample index allowlist from data_dir
    count_file = f"{args.data_dir}JH207_counts.tsv.gz" # TODO change file input so it isn't hard-coded
    allowlist_file = f"{args.data_dir}allowlist.tsv"

    # create and change to output directory, sample_counts
    os.mkdir("../sample_counts")
    os.chdir("../sample_counts")

    try:
        with open(allowlist_file, 'r') as allowlist, gzip.open(count_file, 'r') as counts:
            next(allowlist) # skip header row
            for sample in allowlist:
                index = sample.split()[0] # index to compare to cell barcode
                filename = f"{sample.split()[1]}_{sample.split()[2]}" # filename = sample_name + plate_well
                with open(f"{filename}_counts.tsv", 'w') as outfile:
                    outfile.write("gene\tcell\tcount\n")
                    for line in counts:
                        line = line.decode('utf8')
                        barcode = line.split()[1][-11:] # grab portion of cell barcode that contains sample index
                        if barcode == index:
                            outfile.write(line)
                    counts.seek(0) # rewind counts file for next sample iterations
                outfile.close()

    except Exception as e:
        print(e)
        sys.exit(1)

    finally:
        allowlist.close()
        counts.close()

if __name__ == '__main__':
    main(sys.argv[1:])