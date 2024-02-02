#!/usr/bin/env bash

THREADS=10

# download celltypist markers
# mkdir -p celltypist
# cd celltypist
# wget https://raw.githubusercontent.com/Teichlab/celltypist_wiki/main/atlases/Pan_Immune_CellTypist/v2/encyclopedia/encyclopedia_table.xlsx
# cd ..

# download cellranger
mkdir -p cellranger
cd cellranger

# cellranger 
curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1697253251&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=EIpcT-x-N3~kc51tcZqqAar4wOaJLsKTprE6DwsuKP1-kPr486rCRsWloiLp5jCVnMeeCmmNQ8tI6B5NVsI0D8JLKozXZrM2NFDbcaptbZCTx8zkv2yzmdRghZpg-kJj~avtClEckYc0ceDU6h4GqW1cy1EPCdK6napJgutmfx9N0Ft0IVX6FK351h~tMRNQyzWYim2hYc0N-N5MYbFyzDVw7FgLYu4fBpkCRBepgLXpaBsw~-eptWxByjAW9FPsyDMh0AnGnGDpXS0JCZvimL1HGrLInkaSrwD3q9W4uoMI1EA8QecBX~u7u3DWfJPEGH5BaNDn~FrsYEkS8NVs2w__"

tar -xzf cellranger-7.2.0.tar.gz

CR="$PWD/cellranger-7.2.0/cellranger"

# cellranger arc
curl -o cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1697253284&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjk3MjUzMjg0fX19XX0_&Signature=hYtypbgPWCLZatF2lKxhtxwydM7oK5eyDA9ugypQ4Mes-5JjHXScsnvoUmMgaYBF0mhhEDcyi5XcXi8gtxilx19f4Au10thRjHupVlBMg0xy7BeNIJPIk6pEqITYaOaqwi-UwkNfHv1M-fhfpBxk4UkoSQMiUBf708yhDR4FHmomN3l5H4zwfirPFebldWe~k7rXsW0IeAWIoSt7QE3CgoGOXO2vUwlcvGZOubHu1GlCji89VnpqPdL7W1jsuCRpWcs4G1aO0cj7N3dpaq3uKE8CZrIA2eFpsrbNsyzSmx~nuIFc9X8oCqJx8usnkKLuVSsA-~g0PVYOg9mBo8J1gw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -xzf cellranger-arc-2.0.2.tar.gz

CRarc="$PWD/cellranger-arc-2.0.2/cellranger-arc"

cd ..

# download reference sets
mkdir refs
cd refs

# human ref for cellranger arc
curl -O "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz"
tar -xzf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
REFHUMANARC="$PWD/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

# human ref for cellranger
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
REFHUMANGEX="$PWD/refdata-gex-GRCh38-2020-A"

# mouse ref 
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
tar -xzf refdata-gex-mm10-2020-A.tar.gz
REFMOUSEGEX="$PWD/refdata-gex-mm10-2020-A"

# Pool lab optimized mouse ref
# mouse
curl -o mouse_mm10_optimized_annotation_v2.gtf.gz -L "https://utsw.box.com/shared/static/0tjx19vc07n4fp2r7h2gal0m2tjcxbh2.gz"
gunzip mouse_mm10_optimized_annotation_v2.gtf.gz

# human
curl -o human_GRCh38_optimized_annotation_v2.gtf.gz -L "https://utsw.box.com/shared/static/8k86egl524r9l0nhjj33ja7w4274p45k.gz"
gunzip human_GRCh38_optimized_annotation_v2.gtf.gz

cd ..

# peaks
mkdir peaks
cd peaks
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz
curl -O https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_E18_brain_fresh_5k/atac_v1_E18_brain_fresh_5k_peaks.bed

cd ..

# GEX samples
mkdir datasets
cd datasets

# 10k_PBMC_Multiome_nextgem_Chromium_X
mkdir -p 10k_PBMC_Multiome_nextgem_Chromium_X
cd 10k_PBMC_Multiome_nextgem_Chromium_X

# download files
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_gex_possorted_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.tar.gz
tar -xzf 10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.tar.gz

cd ..

# 10k_PBMC_3p_nextgem_Chromium_X
mkdir -p 10k_PBMC_3p_nextgem_Chromium_X
cd 10k_PBMC_3p_nextgem_Chromium_X

# download files
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_possorted_genome_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.tar.gz
tar -xzf 10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.tar.gz

cd ..

# pbmc_10k_v3
mkdir -p pbmc_10k_v3
cd pbmc_10k_v3

# download files
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz
tar -xzf pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz

cd ..

# human_brain_3k_multiome
mkdir -p human_brain_3k_multiome
cd human_brain_3k_multiome

# curl -o bmmc_site4_donor08_cite.possorted_genome_bam.bam -L
# curl -o bmmc_site4_donor08_cite.possorted_genome_bam.bam -L
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_gex_possorted_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.tar.gz
tar -xzf human_brain_3k_filtered_feature_bc_matrix.tar.gz

cd ..

# e18_mouse_brain_fresh_5k_multiome
mkdir -p e18_mouse_brain_fresh_5k_multiome
cd e18_mouse_brain_fresh_5k_multiome

curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_gex_possorted_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.tar.gz
tar -xzf e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.tar.gz

cd ..

# SC3_v3_NextGem_DI_Nuclei_5K_Multiplex
mkdir -p SC3_v3_NextGem_DI_Nuclei_5K_Multiplex
cd SC3_v3_NextGem_DI_Nuclei_5K_Multiplex

# download and process fastqs
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_Nuclei_5K_Multiplex/SC3_v3_NextGem_DI_Nuclei_5K_Multiplex_fastqs.tar
tar -xf SC3_v3_NextGem_DI_Nuclei_5K_Multiplex_fastqs.tar

$CR count --id=SC3_v3_NextGem_DI_Nuclei_5K_Multiplex --fastqs SC3_v3_NextGem_DI_Nuclei_5K --transcriptome=$REFMOUSEGEX

mv SC3_v3_NextGem_DI_Nuclei_5K_Multiplex/outs/filtered_feature_bc_matrix .
mv SC3_v3_NextGem_DI_Nuclei_5K_Multiplex/outs/possorted_genome_bam.bam .

cd ..

# site4_donor08_cite
mkdir -p bmmc_site4_donor08_cite
cd bmmc_site4_donor08_cite

curl -o bmmc_site4_donor08_cite.possorted_genome_bam.bam -L "https://sra-pub-src-1.s3.amazonaws.com/SRR17693280/site4_donor08_cite.possorted_genome_bam.bam.1"

$CR/lib/bin/bamtofastq --nthreads=$THREADS bmmc_site4_donor08_cite.possorted_genome_bam.bam bmmc_site4_donor08_cite_fastqs
$CR count --id=site4_donor08_cite --fastqs bmmc_site4_donor08_cite_fastqs --transcriptome=$REFHUMANGEX

mv bmmc_site4_donor08_cite/outs/filtered_feature_bc_matrix .
mv bmmc_site4_donor08_cite/outs/possorted_genome_bam.bam .

cd ..

# site4_donor08_multiome
mkdir -p bmmc_site4_donor08_multiome
cd bmmc_site4_donor08_multiome

curl -o bmmc_site4_donor08_multiome_gex.possorted_genome_bam.bam -L https://sra-pub-src-2.s3.amazonaws.com/SRR17693267/site4_donor08_multiome_gex.possorted_genome_bam.bam.1
curl -o bmmc_site4_donor08_multiome_atac.possorted_genome_bam.bam -L https://sra-pub-src-2.s3.amazonaws.com/SRR17693254/site4_donor08_multiome_atac.possorted_genome_bam.bam.1

$CR/lib/bin/bamtofastq --nthreads=$THREADS bmmc_site4_donor08_multiome_gex.possorted_genome_bam.bam bmmc_site4_donor08_multiome_gex_fastqs
$CR/lib/bin/bamtofastq --nthreads=$THREADS bmmc_site4_donor08_multiome_atac.possorted_genome_bam.bam bmmc_site4_donor08_multiome_atac_fastqs


cat << EOF > libraries.csv
fastqs,sample,library_type
bmmc_site4_donor08_multiome_atac_fastqs,bamtofastq,Chromatin Accessibility
bmmc_site4_donor08_multiome_gex_fastqs,bamtofastq,Gene Expression
EOF

$CRarc count --id=bmmc_site4_donor08_multiome --reference=$REFHUMANARC --libraries=libraries.csv --localcores=$THREADS

mv bmmc_site4_donor08_multiome/outs/filtered_feature_bc_matrix .
mv bmmc_site4_donor08_multiome/outs/atac_peaks.bed ../peaks/bmmc_site4_donor08_multiome_peaks.bed
mv bmmc_site4_donor08_multiome/outs/gex_possorted_bam.bam .
mv bmmc_site4_donor08_multiome/outs/atac_fragments.tsv.gz .
mv bmmc_site4_donor08_multiome/outs/atac_fragments.tsv.gz.tbi .



