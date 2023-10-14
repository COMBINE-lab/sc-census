#!/usr/bin/env python3

import pysam
import os
import shutil
import scipy
import gzip
import pandas as pd
import pickle
import argparse

##    1. load intragenic BAM file and get the intragenic UMIs of each cell
##    2. load the intergenic BED file as a dict(Read name: (the closest peak, distance))
##    3. load the intergenic BAM file, for each alignment, if the read is in the dict, and its CB-UMI is not an intragenic UMI, put it in the dict((CB, UMI): (peak,distance))
##    4. for each record in dict((CB, UMI): (peak,distance)), add 1 to the peak of the CB that has the minimum distance 

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("in_bam", help="input bam file", type=str)
parser.add_argument("in_tsv", help="input read_name-closest_peak-distance tsv file", type=str)
parser.add_argument("in_pkl", help="input CB-UMI list dict in a pickle file", type=str)
parser.add_argument("out_dir", help="output count matrix dir", type=str)
parser.add_argument("--intragenic", help="include UMIs with intragenic alignments in the count matrix", action="store_true")
args = parser.parse_args()
# args = parser.parse_args(["bam.bam","read_peak_distance.tsv","cb-umi_dict.pickle","ocr_count_intergenic"])

# if the output directory exists, we delete it
if os.path.exists(args.out_dir):
    shutil.rmtree(args.out_dir)

# make a new one    
os.mkdir(args.out_dir)

# load the CB-UMI list dict
cb_umi_dict = pickle.load(open(args.in_pkl, "rb"))

# check if the dict is empty
if len(cb_umi_dict) == 0:
    raise Exception("The input CB-UMI list dict is empty!")

# load tsv file
read_peak_distance_dict = {}
with open(args.in_tsv, "r") as f:
    for line in f:
        line = line.strip().split("\t")
        read_peak_distance_dict[line[0]] = (line[1], int(line[2]))

# check if the dict is empty
if len(read_peak_distance_dict) == 0:
    raise Exception("The input tsv file is empty!")

# parse the bam file
if not os.path.exists(args.in_bam+".bai"):
    pysam.index(args.in_bam)
    
in_bam = pysam.AlignmentFile(args.in_bam, "rb")

peak_set = set()
cb_umi_peak_dict = {cb: {} for cb in cb_umi_dict.keys()}
for read in in_bam:
    # the read has to have a peak assigned
    peak_distance_tuple = read_peak_distance_dict.get(read.query_name) 
    if peak_distance_tuple is not None:
        cb = read.get_tag("CB")
        umi = read.get_tag("UB")
            
        # we only process UMIs that are not intragenic
        if args.intragenic or umi not in cb_umi_dict.get(cb):
            # if the CB-UMI pair is not in the dict, we add it
            # else, we update it if the distance is smaller
            if cb_umi_peak_dict[cb].get(umi) is None:
                cb_umi_peak_dict[cb][umi] = peak_distance_tuple
            else:
                if cb_umi_peak_dict[cb][umi][1] > peak_distance_tuple[1]:
                    cb_umi_peak_dict[cb][umi] = peak_distance_tuple
            
            # add the peak to the set
            peak_set.add(peak_distance_tuple[0])
in_bam.close()

# for (cb, umi_peak_dict) in cb_umi_peak_dict.items():
#     print(cb, len(umi_peak_dict))

# generate the count matrix
# first, we need to create an empty pandas dataframe
cb_list = list(cb_umi_dict.keys())
df = pd.DataFrame(0, index=list(peak_set), columns=cb_list)

# then, we traverse the dict to fill in the count matrix
for (cb, umi_peak_dict) in cb_umi_peak_dict.items():
    for (umi, peak_distance_tuple) in umi_peak_dict.items():
        df.at[peak_distance_tuple[0],cb] += 1

index_df = df.index.to_frame(name="feature_name")
index_df["feature_id"] = index_df.index
index_df["feature_type"] = "Peaks"
index_df["chr"] = index_df["feature_name"].apply(lambda x: x.split(":")[0])
index_df[["start", "end"]] = index_df["feature_name"].apply(lambda x: x.split(":")[1]).str.split("-", expand=True)

index_df.to_csv(os.path.join(args.out_dir, "features.tsv.gz"), header=False, index=False, sep="\t", compression="gzip")

df.columns.to_frame().to_csv(os.path.join(args.out_dir, "barcodes.tsv.gz"), header=False, index=False, sep="\t", compression="gzip")

with gzip.open(os.path.join(args.out_dir, "matrix.mtx.gz"), 'wb') as f:
    scipy.io.mmwrite(f, scipy.sparse.csr_matrix(df.values), comment='\n intergenic peak count using intergenic reads.\n')
    
