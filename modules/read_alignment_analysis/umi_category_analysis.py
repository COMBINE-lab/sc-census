# mamba create -n steve -c bioconda -c conda-forge -c default pysam scanpy pickle argparse upsetplot matplotlib ipython -y

import pysam
import os
import pickle
import gzip
import itertools
import pandas as pd
import argparse
from upsetplot import from_contents
from upsetplot import UpSet
from matplotlib import pyplot

def create_dict_from_bam(bam_filename_dict: list):    
    # initialize an empty dictionary
    cb_umi_dict = {}
    cat_cbumi_dict = {}
    cat_read_dict = {}
    # iterate over each BAM file
    for (cat_name,in_bam_file) in bam_filename_dict.items():
        cbumi_set = set()
        read_set = set()
        # Open the BAM file for reading
        if not os.path.exists(in_bam_file+".bai"):
            pysam.index(in_bam_file)
        in_bam = pysam.AlignmentFile(in_bam_file, "rb")
        # Iterate over each record in the BAM file
        for read in in_bam:
            # Get the name of the read
            read_name = read.query_name
            cb = read.get_tag('CB')
            umi = read.get_tag('UB')
            if cb not in cb_umi_dict:
                cb_umi_dict[cb] = set()
            cb_umi_dict[cb].add(umi)
            
            cbumi_set.add((cb,umi))
            read_set.add(read_name)
        # Close the BAM files
        in_bam.close()
        cat_cbumi_dict[cat_name] = cbumi_set
        cat_read_dict[cat_name] = read_set
    return (cat_cbumi_dict, cat_read_dict, cb_umi_dict)


parser = argparse.ArgumentParser()
parser.add_argument("species", type = str, help="sample species")
parser.add_argument("sample_type", type = str, help="sample type")
parser.add_argument("sample_name", type = str, help="sample name")
parser.add_argument("out_dir", default="read_alignment_analysis", type = str, help="Path to the output directory")
parser.add_argument("exon_exon_junction_bam", default="exon_exon_junctions.bam", type = str, help="Path to the exon-exon junction bam file")
parser.add_argument("intron_exon_junction_bam", default="intron_exon_junctions.bam", type = str, help="Path to the intron-exon junction bam file")
parser.add_argument("transcripts_terminal_exon_bam", default="transcripts_terminal_exons.bam", type = str, help="Path to the transcript terminal exon bam file")
parser.add_argument("spliced_txp_terminal_kilobase_bam", default="spliced_txp_three_prime_terminal_kilobase.bam", type = str, help="Path to the spliced transcript terminal kilobase bam file")
parser.add_argument("spliced_txp_terminal_kilobase_exonic_bam", default="spliced_txp_three_prime_terminal_kilobase_exonic.bam", type = str, help="Path to the spliced transcript terminal kilobase exonic bam file")
parser.add_argument("unspliced_txp_terminal_kilobase_bam", default="unspliced_txp_three_prime_terminal_kilobase.bam", type = str, help="Path to the unspliced transcript terminal kilobase bam file")
parser.add_argument("exon_bam", default="exons.bam", type = str, help="Path to the exon bam file")
parser.add_argument("introns_bam", default="introns.bam", type = str, help="Path to the intron bam file")
args = parser.parse_args()

bam_filename_dict = {
    "E-E-junc": args.exon_exon_junction_bam,
    "I-E-junc": args.intron_exon_junction_bam,
    "tx-last-exon": args.transcripts_terminal_exon_bam,
    "S-tx-last-kb": args.spliced_txp_terminal_kilobase_bam,
    "U-tx-last-kb": args.unspliced_txp_terminal_kilobase_bam,
    "E": args.exon_bam,
    "I": args.introns_bam
}

fig_dir = os.path.join(args.out_dir, "read_and_umi_feature_category_upset_plot")

if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

(cat_cbumi_dict, cat_read_dict, cb_umi_dict) = create_dict_from_bam(bam_filename_dict)

# plot the upset plot
# UMI

num_umis = 0
for (cat_name, read_set) in cat_cbumi_dict.items():
    num_umis += len(read_set)

in_data = from_contents(cat_cbumi_dict)
UpSet(in_data, subset_size='count').plot()
# pyplot.show()
pyplot.title("".join(["The compatible feature categories of UMIs (Total # of UMIs: ", str(num_umis), ")"]))
pyplot.savefig(os.path.join(fig_dir, "UMI_category_upsetplot.pdf"))

# read

num_reads = 0
for (cat_name, read_set) in cat_read_dict.items():
    num_reads += len(read_set)

in_data = from_contents(cat_read_dict)
UpSet(in_data, subset_size='count',show_percentages = True).plot()
# pyplot.show()
pyplot.title("".join(["The compatible feature categories of Reads (Total # of reads: ", str(num_reads), ")"]))
pyplot.savefig(os.path.join(fig_dir , "read_category_upsetplot.pdf"))

pickle.dump(cb_umi_dict, open(os.path.join(args.out_dir, "cb-umi_dict.pickle"), "wb"))
pickle.dump(cat_cbumi_dict, open(os.path.join(args.out_dir, "cat-cbumi_dict.pickle"), "wb"))
pickle.dump(cat_read_dict, open(os.path.join(args.out_dir, "cat-read_dict.pickle"), "wb"))









