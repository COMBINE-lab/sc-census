#! /usr/bin/env python3
import pysam
import os
import pickle

def create_dict_from_bam(in_bam_file: str, out_pickle_file: str):
    # initialize an empty dictionary
    cb_umi_dict = {}

    # open bam
    if not os.path.exists(in_bam_file+".bai"):
        pysam.index(in_bam_file)
    in_bam = pysam.AlignmentFile(in_bam_file, "rb")

    # Iterate over each record in the BAM file
    for read in in_bam:
        # we process the read if it has a valid CB and UB tag
        if read.has_tag("CB") and read.has_tag("UB"):
            cb = read.get_tag("CB")
            if cb not in cb_umi_dict:
                cb_umi_dict[cb] = set()
            cb_umi_dict.get(cb).add(read.get_tag("UB"))
    # Close the BAM files
    in_bam.close()
    pickle.dump(cb_umi_dict, open(out_pickle_file, "wb"))
    return cb_umi_dict

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("in_bam", help="input bam file", type=str)
parser.add_argument("out_pickle", help="output pickle file", default="cb-umi_dict.pickle", type=str)
args = parser.parse_args()
create_dict_from_bam(args.in_bam, args.out_pickle)