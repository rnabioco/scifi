# fastq_split.py
# Take demultiplexed, merged .bam files from Datlinger et al. (2021) (GSE168620) and reformat to Illumina NovaSeq format
# Using SRR13925098 for test data (PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells; sample wells B13, B15, F13, F15)
#
# Preprocessing: 
#   - `samtools view -b -r READ_GROUP PD213_scifi_2_CRISPR_TCR_77300_MeOH_cells.bam` to extract specific well samples
#       - READ_GROUP = ID value of @RG; ex. `BSF_0774_HNNGMDMXX_1#PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells_B13_01`
#   - `samtools fastq -T BC,RX,r2 PD213_scifi_2_CRISPR_TCR_77300_MeOH_cells.bam`
#       - r2: 16 base round 2 (10X ATAC oligos) barcode
#       - BC: 13 base round 1 (well-specific) barcode + 8 base i7 sample index
#       - RX: 8 base UMI
#
# Expected fastq format:
#   - R1: 8 base UMI (RX) + 13 base round 1 barcode (first 13 bases of BC)
#   - R2: 16 base round 2 barcode (r2)
#   - R3: 70 base mRNA sequence
#   - I1: 8 base i7 sample index (last 8 bases of BC)
#
# Note: Using dummy quality scores for R2, R3, and I1 .fastq files.
#
# Formatting example:
# Input .fastq:
#   @HWI-A00245_BSF_0774:1:1101:10013:6684#PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells_B13_01	BC:Z:ATGTAAGAATCAATCCTGAGC	RX:Z:ATTATGTA	r2:Z:GAGATTCGTGCTAGTT
#   ACTCTCTTCTCCCCAAAAATTGAGAACAAGGCAAGAATATCTACTTTCATCACTCCTATTCAGTTTTATA
#   +
#   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Output R1.fastq:
#   @HWI-A00245_BSF_0774:1:1101:10013:6684#PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells_B13_01
#   ATTATGTAATGTAAGAATCAA
#   +
#   FFFFFFFFFFFFFFFFFFFFF
# Output R2.fastq:
#   @HWI-A00245_BSF_0774:1:1101:10013:6684#PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells_B13_01
#   GAGATTCGTGCTAGTT
#   +
#   FFFFFFFFFFFFFFFF
# Output R3.fastq:
#   @HWI-A00245_BSF_0774:1:1101:10013:6684#PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells_B13_01
#   ACTCTCTTCTCCCCAAAAATTGAGAACAAGGCAAGAATATCTACTTTCATCACTCCTATTCAGTTTTATA
#   +
#   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Output I1.fastq:
#   @HWI-A00245_BSF_0774:1:1101:10013:6684#PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells_B13_01
#   TCCTGAGC
#   +
#   FFFFFFFF

import sys
import argparse
import os
import glob
from itertools import islice

# Specifying command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("data_dir", help="path to folder with preprocessed .fastq files")

args = parser.parse_args()

def getchunk(file): # helper function to go through .fastq file 4 lines at a time
    return [x.strip() for x in islice(file, 4)]

def main(argv):
    # change dir to data_dir
    os.chdir(args.data_dir)

    for current_inpath in glob.glob("*.fastq"):
        try:
            # open .fastq file
            infile = open(current_inpath)

            # open outfiles as {sample}_R1.fastq, {sample}_R2.fastq, etc.
            sample = current_inpath.split(".")[0]
            R1_file = open(f"{sample}_R1.fastq", "w")
            R2_file = open(f"{sample}_R2.fastq", "w")
            R3_file = open(f"{sample}_R3.fastq", "w")
            I1_file = open(f"{sample}_I1.fastq", "w")

            while True:
                chunk = getchunk(infile)
                if not chunk:
                    break

                # get header; this is the same for all outfiles    
                header = chunk[0].split()[0]

                # get R1 sequence (8 base UMI (RX) + 13 base round 1 barcode (first 13 bases of BC))
                R1_seq = f"{chunk[0].split()[2][5:]}{chunk[0].split()[1][5:18]}"

                # get R2 sequence (16 base round 2 barcode (r2))
                R2_seq = f"{chunk[0].split()[3][5:]}"

                # get I1 sequence (8 base i7 sample index (last 8 bases of BC))
                I1_seq = f"{chunk[0].split()[1][18:]}"

                R1_file.write(f"{header}\n")
                R1_file.write(f"{R1_seq}\n")
                R1_file.write(f"+\n")
                R1_file.write(f"FFFFFFFFFFFFFFFFFFFFF\n")

                R2_file.write(f"{header}\n")
                R2_file.write(f"{R2_seq}\n")
                R2_file.write(f"+\n")
                R2_file.write(f"FFFFFFFFFFFFFFFF\n")

                # the R3 sequence is the mRNA sequence included in the infile
                R3_file.write(f"{header}\n")
                R3_file.write(f"{chunk[1]}\n")
                R3_file.write(f"{chunk[2]}\n")
                R3_file.write(f"{chunk[3]}\n")

                I1_file.write(f"{header}\n")
                I1_file.write(f"{I1_seq}\n")
                I1_file.write(f"+\n")
                I1_file.write(f"FFFFFFFF\n")
                
        except IOError:
            print("Could not open file.")
            sys.exit(1)

        except:
            print("Error.")
            sys.exit(1)

        finally: # close the files
            infile.close()
            R1_file.close()
            R2_file.close()
            R3_file.close()
            I1_file.close()

if __name__ == '__main__':
    main(sys.argv[1:])