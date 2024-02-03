# sc-census
scCensus is an analysis pipeline focusing on establishing the interpretability of off-target priming reads in scRNA-seq. You can find the preprint at https://www.biorxiv.org/content/10.1101/2024.01.29.577807v1. 

To run the pipeline, you need to have Nextflow and conda installed. If not, please follow the instructions [here](https://www.nextflow.io/docs/latest/getstarted.html) to install Nextflow and [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install conda.

As the pipeline involves analyzing more than 10 single-cell datasets and requires a large amount of time, computational resources, and disk space (3 TB), we recommend running the pipeline on a high-performance computing (HPC) cluster. We also provide the analysis results of the pipeline on [Zenodo](https://doi.org/10.5281/zenodo.10520670), which can be downloaded and visualized directly.


## Running the pipeline
The pipeline can be devided into the following steps:
1. **Prepare the input data**. In this step, we download the reference genome, the selected single cell datasets, required software, to generate the input data for the pipeline. This step is done by running the bash script `data_prep.sh` in the `data` folder.
2. **Transcript compostion analysis**: In this analysis, we explore some of the interesting characteristics of human and mouse transcripts. For example, we wnat to know how many transcripts have a terminal exon that is longer than 1000 base pairs. For those transcripts, priming the polyA tail will be almost impossible to generate reads that have a definitive splicing status, such as exon-exon junctional reads or intronic reads. 
3. **Read alignment analtysis**： In this analysis, we analyze the alignment of the reads. We want to know how many reads map to exons, introns, transcripts' terminal exons, etc. We want to know how many reads are off-target priming reads.
4. **Quantification**: In this analysis, we quantify the expression of the defined features in intragenic regions and intergenic regions. We consider sense, antisense, and intergenic reads. 
5. **Downstream analysis**: In this analysis, we perform the downstream analysis, such as differential expression analysis, gene set enrichment analysis, for sense, antisense, and intergenic reads.

The actual code for running the pipeline is the following:
1. we need to execute the bash script `data_prep.sh` in the `data` folder to prepare the input data for the pipeline. **Throughout this section, we suppose that the root working directory is the `sc-census` folder.**

```bash
cd data
# !! We have to stay in the data folder to run the bash script !!
bash data_prep.sh
cd ..
```

2. Then, we want to setup the config file `nextflow.config` in the root working directory. If we want to use the default configuration, we can just modify the `process` to fit our HPC cluster. The parameter that has to be modified is `clusterOptions`, where you have to specify the slurm setting of the clsuter we want to use for executing the pipeline. 

3. Finally, we call the pipeline with the following command:

```bash
nextflow main.nf --resume
```

We specify the `--resume` option to run the pipeline from the last checkpoint. It doesn't hurt if this is the first time we run the pipeline, but if the pipeline fails, specifying it will use catched results instead of computing them one more time. If we want to run the pipeline from the beginning, we can just remove the `--resume` option.








