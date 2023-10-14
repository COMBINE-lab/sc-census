#! /usr/bin/env python3
import scanpy as sc
import pandas as pd
import pyroe
import celltypist
from celltypist import models
import os
import shutil

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("sample_type", type=str)
parser.add_argument("cells_or_nuclei", type=str)
parser.add_argument("quant_dir", help="alevin-fry quant dir", type=str)
parser.add_argument("output_dir", help="output dir", type=str)
args = parser.parse_args()
# args = parser.parse_args(["human_pbmc", "cell", ".", "."])

# if the output directory exists, we delete it
if os.path.exists(args.out_dir):
    shutil.rmtree(args.out_dir)

# make a new one    
os.mkdir(args.out_dir)

output_format = 'snrna'
if args.cells_or_nuclei == "cell":
    output_format = "scrna"

# we read in the id to name mapping tsv file (no header), and convert it to a dict
gid2name = {id: n for (id, n) in pd.read_csv("gene_id_to_name.tsv", sep="\t",header=None).values}

adata = pyroe.load_fry(args.quant_dir, output_format=output_format)

sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=3)


adata.var_names = [gid2name.get(id) for id in adata.var_names]

adata.var_names_make_unique() 

sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sample_type_to_model_name = {"human_pbmc": "Immune_All_Low.pkl","human_bmmc": "Immune_All_Low.pkl", "human_brain": "Developing_Human_Brain.pkl", "mouse_brain": "Developing_Mouse_Brain.pkl"}

model_name = sample_type_to_model_name[args.sample_type]

models.download_models(model=model_name)

model = models.Model.load(model = model_name)

predictions = celltypist.annotate(adata, model = model_name, majority_voting = True)

predictions.predicted_labels.to_csv(os.path.join(args.out_dir,"celltypist_predictions.csv"))




