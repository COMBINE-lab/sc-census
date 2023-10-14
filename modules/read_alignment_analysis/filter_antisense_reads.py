#! /usr/bin/env python3
import pysam
import os
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("in_bam", help="input bam file", type=str)
parser.add_argument("in_pkl", help="in pickle file", type=str)
parser.add_argument("out_bam", help="out bam file", type=str)
args = parser.parse_args()

# load the CB-UMI list dict
cb_umi_dict = pickle.load(open(args.in_pkl, "rb"))

# open bam
if not os.path.exists(args.in_bam+".bai"):
    pysam.index(args.in_bam)
in_bam = pysam.AlignmentFile(args.in_bam, "rb")

# open out_bam
out_bam = pysam.AlignmentFile(args.out_bam, "wb", template=in_bam)

# Iterate over each record in the BAM file
for read in in_bam:
    # if the cb-umi pair is in the dict, write it to the output bam
    if read.get_tag("UB") not in cb_umi_dict[read.get_tag("CB")]:
        # add the read to the output bam
        out_bam.write(read)

# Close the BAM files
in_bam.close()
