SAMPLE_NAME="10k_PBMC_Multiome_nextgem_Chromium_X"
DATA_DIR="/mnt/scratch4/dongze/read_lives_matter/data/10k_PBMC_Multiome_nextgem_Chromium_X_fastqs"
READ_ALIGNMENT_ANALYSIS_DIR="/mnt/scratch4/dongze/read_lives_matter/read_alignment_analysis/$SAMPLE_NAME/out_bam"
READ_ALIGNMENT_ANALYSIS_TEMP_DIR="$READ_ALIGNMENT_ANALYSIS_DIR/temp_files"

THREADS=10
GTF_PATH="/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
GENOME_PATH="/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa"

PEAK_BED_PATH="$DATA_DIR/10k_PBMC_Multiome_nextgem_Chromium_X_atac_peaks.bed"
INTERGENIC_BAM_PATH="$READ_ALIGNMENT_ANALYSIS_DIR/${SAMPLE_NAME}_intergenic.bam"
ANTISENSE_BAM_PATH="$READ_ALIGNMENT_ANALYSIS_DIR/${SAMPLE_NAME}_antisense.bam"
TXOME_BAM_PATH="$READ_ALIGNMENT_ANALYSIS_TEMP_DIR/${SAMPLE_NAME}_high_quality_cells_sense_transcriptome.bam"
ALL_BAM_PATH="$READ_ALIGNMENT_ANALYSIS_TEMP_DIR/${SAMPLE_NAME}_high_quality_cells.bam"
FASTQ_DIR="$DATA_DIR//mnt/scratch4/dongze/read_lives_matter/data/10k_PBMC_Multiome_nextgem_Chromium_X_fastqs/10k_PBMC_Multiome_nextgem_Chromium_X_gex"
FCB_DIR="/mnt/scratch4/dongze/read_lives_matter/feature_category_bed"

MAKE_CB_UMI_LIST_SCRIPT="/mnt/scratch4/dongze/read_lives_matter/code/make_cb-umi_list.py"
GENERATE_COUNT_MATRIX_SCRIPT="/mnt/scratch4/dongze/read_lives_matter/code/generate_count_matrix.py"

STAR_INDEX_DIR="/mnt/scratch3/alevin_fry_submission/indices/human-2020A/star_index"
OUT_DIR="/mnt/scratch4/dongze/read_lives_matter/quant_results/$SAMPLE_NAME"
TEMP_DIR="$OUT_DIR/temp"

FANTOM5_BED_PATH="/mnt/scratch4/dongze/read_lives_matter/data/fantom5/hg38_liftover+new_CAGE_peaks_phase1and2.bed"

OCR_DIR="$OUT_DIR/open_chomatin_regions"
mkdir $OCR_DIR

mkdir -p $TEMP_DIR

#########
# spliceu
#########
SPLICEU_DIR="$OUT_DIR/spliceu_simpleaf"
mkdir -p $SPLICEU_DIR

# index
simpleaf index \
--output $SPLICEU_DIR/simpleaf_index \
--threads $THREADS \
--fasta $GENOME_PATH \
--gtf $GTF_PATH \
--ref-type spliceu \
--use-piscem

# quant
# we first convert BAM to FASTQs
# sense transcriptomic reads
bamtofastq \
--nthreads=$THREADS \
$TXOME_BAM_PATH \
$TEMP_DIR/sense_txome_fastqs

# Define filename pattern
reads1_pat="_R1_"
reads2_pat="_R2_"

# Obtain and sort filenames
reads1="$(find -L $TEMP_DIR/sense_txome_fastqs -name "*$reads1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
reads2="$(find -L $TEMP_DIR/sense_txome_fastqs -name "*$reads2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"

simpleaf quant \
--chemistry 10xv3 \
--output $SPLICEU_DIR/simpleaf_quant \
--threads $THREADS \
--index $SPLICEU_DIR/simpleaf_index/index \
--reads1 $reads1 \
--reads2 $reads2 \
--unfiltered-pl \
--resolution cr-like

###################
# antisense spliceu 
###################
ANTISENSE_DIR="$OUT_DIR/antisense_simpleaf"
mkdir -p $ANTISENSE_DIR

# first convert BAM to FASTQs
bamtofastq \
--nthreads=$THREADS \
$ANTISENSE_BAM_PATH \
$TEMP_DIR/antisense_fastqs

# Define filename pattern
reads1_pat="_R1_"
reads2_pat="_R2_"

# Obtain and sort filenames
reads1="$(find -L $TEMP_DIR/antisense_fastqs -name "*$reads1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
reads2="$(find -L $TEMP_DIR/antisense_fastqs -name "*$reads2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"

#quant
simpleaf quant \
--chemistry 10xv3 \
--output $ANTISENSE_DIR/simpleaf_quant \
--threads $THREADS \
--index $SPLICEU_DIR/simpleaf_index/index \
--expected-ori rc \
--reads1 $reads1 \
--reads2 $reads2 \
--unfiltered-pl \
--resolution cr-like


#########################################################
# open chromatin regions
# non-coding RNAs
#########################################################

# get intragenic UMIs
$MAKE_CB_UMI_LIST_SCRIPT \
$TXOME_BAM_PATH \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl

# we get protein coding genes' regions
awk '/gene_type "protein_coding";/ && $3 == "gene" {print $1"\t"$4"\t"$5"\t"".""\t"$6"\t"$7}' $GTF_PATH > $OCR_DIR/protein_coding_genes.bed

#-----------------------------------#
# alignment set #1: intergenic reads
#-----------------------------------#
# 1. get intergenic peaks
bedtools intersect \
-a $PEAK_BED_PATH \
-b $FCB_DIR/genes.bed \
-wa -v > $OCR_DIR/${SAMPLE_NAME}_intergenic_peaks.bed

## 2. get intergenic alingments
bedtools intersect \
-a $ALL_BAM_PATH \
-b $FCB_DIR/genes.bed \
-wa -v > $OCR_DIR/${SAMPLE_NAME}_intergenic.bam

## 3. find intergenic alignments' closest intergenic peak
bedtools closest \
-d -t first \
-a $OCR_DIR/${SAMPLE_NAME}_intergenic.bam \
-b $OCR_DIR/${SAMPLE_NAME}_intergenic_peaks.bed | awk '$16 <= 9000 && $16 >= 0 {print $4"\t"$13":"$14"-"$15"\t"$16}' - > $OCR_DIR/${SAMPLE_NAME}_intergenic_read_closest_intergenic_peak_distance.tsv

## 4. generate count matrix
# generate_count_matrix.py in_bam in_tsv in_pkl out_dir
$GENERATE_COUNT_MATRIX_SCRIPT \
$OCR_DIR/${SAMPLE_NAME}_intergenic.bam \
$OCR_DIR/${SAMPLE_NAME}_intergenic_read_closest_intergenic_peak_distance.tsv \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl \
$OCR_DIR/${SAMPLE_NAME}_intergenic_count_matrix

#-------------------------------------#
# alignment set #2: not coding regions
#-------------------------------------#

# 1. get non-protein coding peaks
bedtools intersect \
-a $PEAK_BED_PATH \
-b $OCR_DIR/protein_coding_genes.bed \
-wa -v > $OCR_DIR/${SAMPLE_NAME}_non_coding_peaks.bed

# 2. get non-coding alignments
bedtools intersect \
-a $ALL_BAM_PATH \
-b $OCR_DIR/protein_coding_genes.bed \
-wa -v > $OCR_DIR/${SAMPLE_NAME}_non_coding.bam

# 3. find non-coding alignments' closest non-coding peak
bedtools closest \
-d -t first \
-a $OCR_DIR/${SAMPLE_NAME}_non_coding.bam \
-b $OCR_DIR/${SAMPLE_NAME}_non_coding_peaks.bed | awk '$16 < 8000 && $16 >= 0 {print $4"\t"$13":"$14"-"$15"\t"$16}' - > $OCR_DIR/${SAMPLE_NAME}_non_coding_read_closest_non_coding_peak_distance.tsv

## 4. generate count matrix
# generate_count_matrix.py in_bam in_tsv in_pkl out_dir
$GENERATE_COUNT_MATRIX_SCRIPT \
$OCR_DIR/${SAMPLE_NAME}_non_coding.bam \
$OCR_DIR/${SAMPLE_NAME}_non_coding_read_closest_non_coding_peak_distance.tsv \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl \
$OCR_DIR/${SAMPLE_NAME}_non_coding_count_matrix

#-------------------------------------------#
# alignment set #3: not sense-coding regions
#-------------------------------------------#
# 1. get non-coding alignments
bedtools intersect \
-a $ALL_BAM_PATH \
-b $OCR_DIR/protein_coding_genes.bed \
-wa -s -v > $OCR_DIR/${SAMPLE_NAME}_non_sense_coding.bam

# 2. find non-coding alignments' closest non-coding peak
bedtools closest \
-d -t first \
-a $OCR_DIR/${SAMPLE_NAME}_non_sense_coding.bam \
-b $PEAK_BED_PATH | awk '$16 < 8000 && $16 >= 0 {print $4"\t"$13":"$14"-"$15"\t"$16}' - > $OCR_DIR/${SAMPLE_NAME}_non_sense_coding_peak_distance.tsv

## 3. generate count matrix
# generate_count_matrix.py in_bam in_tsv in_pkl out_dir
$GENERATE_COUNT_MATRIX_SCRIPT \
$OCR_DIR/${SAMPLE_NAME}_non_sense_coding.bam \
$OCR_DIR/${SAMPLE_NAME}_non_sense_coding_peak_distance.tsv \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl \
$OCR_DIR/${SAMPLE_NAME}_non_sense_coding_count_matrix


#-------------------------------------------#
# alignment set #4: everything
#-------------------------------------------#
# 1. find alignments' closest non-coding peak
bedtools closest \
-d -t first \
-a $ALL_BAM_PATH \
-b $PEAK_BED_PATH | awk '$16 < 8000 && $16 >= 0 {print $4"\t"$13":"$14"-"$15"\t"$16}' - > $OCR_DIR/${SAMPLE_NAME}_peak_distance.tsv

## 2. generate count matrix
# generate_count_matrix.py in_bam in_tsv in_pkl out_dir
$GENERATE_COUNT_MATRIX_SCRIPT \
$ALL_BAM_PATH \
$OCR_DIR/${SAMPLE_NAME}_peak_distance.tsv \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl \
$OCR_DIR/${SAMPLE_NAME}_non_sense_coding_count_matrix \
--intragenic

#---------------------------------------------------#
# alignment set #5: intergenic in FANTOM5 enhancers
#---------------------------------------------------#

# 1. find alignments' closest non-coding peak

bedtools sort -i $FANTOM5_BED_PATH > $TEMP_DIR/fantom5_sorted.bed

bedtools closest \
-d -t first \
-a $OCR_DIR/${SAMPLE_NAME}_intergenic.bam \
-b $TEMP_DIR/fantom5_sorted.bed | awk '$16 < 8000 && $16 >= 0 {print $4"\t"$13":"$14"-"$15"\t"$16}' - > $OCR_DIR/${SAMPLE_NAME}_intergenic_fantom_distance.tsv

## 2. generate count matrix
# generate_count_matrix.py in_bam in_tsv in_pkl out_dir
$GENERATE_COUNT_MATRIX_SCRIPT \
$ALL_BAM_PATH \
$OCR_DIR/${SAMPLE_NAME}_intergenic_fantom_enhancer_distance.tsv \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl \
$OCR_DIR/${SAMPLE_NAME}_fantom_enhancer_matrix \
--intragenic

#---------------------------------------------------#
# alignment set #6: everything but FANTOM5 enhancers
#---------------------------------------------------#

# 1. find alignments' closest non-coding peak
bedtools closest \
-d -t first \
-a $ALL_BAM_PATH \
-b $TEMP_DIR/fantom5_sorted.bed | awk '$16 < 8000 && $16 >= 0 {print $4"\t"$13":"$14"-"$15"\t"$16}' - > $OCR_DIR/${SAMPLE_NAME}_fantom_distance.tsv

## 2. generate count matrix
# generate_count_matrix.py in_bam in_tsv in_pkl out_dir
$GENERATE_COUNT_MATRIX_SCRIPT \
$ALL_BAM_PATH \
$OCR_DIR/${SAMPLE_NAME}_fantom_enhancer_distance.tsv \
$OCR_DIR/${SAMPLE_NAME}_sense_txome_CB_UMI_list.pkl \
$OCR_DIR/${SAMPLE_NAME}_fantom_enhancer_matrix \
--intragenic




