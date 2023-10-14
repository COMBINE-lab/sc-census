# sample_name="10k_PBMC_Multiome_nextgem_Chromium_X"
# DATA_DIR="/mnt/scratch4/dongze/read_lives_matter/data/10k_PBMC_Multiome_nextgem_Chromium_X_fastqs"
# fcb_dir="/mnt/scratch4/dongze/read_lives_matter/feature_category_bed"
# OUT_DIR="/mnt/scratch4/dongze/read_lives_matter/read_alignment_analysis/$sample_name/out_bam"
# TEMP_DIR="$OUT_DIR/temp_files"
# in_bam="$DATA_DIR/10k_PBMC_Multiome_nextgem_Chromium_X_gex_possorted_bam.bam"
# CB_LIST_PATH="$DATA_DIR/filtered_feature_bc_matrix/barcodes.tsv.gz"
# MAKE_CB_UMI_LIST_SCRIPT="/mnt/scratch4/dongze/read_lives_matter/code/make_cb-umi_list.py"
# GENERATE_COUNT_MATRIX_SCRIPT="/mnt/scratch4/dongze/read_lives_matter/code/generate_count_matrix.py"
# num_threads=10

mkdir -p $OUT_DIR
mkdir -p $TEMP_DIR

gunzip -c $CB_LIST_PATH > $TEMP_DIR/barcodes.txt

####################################################
# Preprocess: split splice and contiguous alignments
####################################################
# filter only high quality cell's reads and mapped reads
samtools view -@ $num_threads -bh -F 4 -D CB:$TEMP_DIR/barcodes.txt $in_bam | samtools view -@ $num_threads -bh -d UB -o $TEMP_DIR/${sample_name}_high_quality_cells.bam -

# split sense and antisense alignments
## sense
bedtools intersect \
-a $TEMP_DIR/${sample_name}_high_quality_cells.bam \
-b $fcb_dir/genes.bed \
-wa \
-f 1 \
-s \
> $TEMP_DIR/${sample_name}_high_quality_cells_sense_transcriptome.bam

# split splice and contiguous alignments
samtools view -@ $num_threads -h $TEMP_DIR/${sample_name}_high_quality_cells_sense_transcriptome.bam | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -@ $num_threads -b > $TEMP_DIR/${sample_name}_spliced_alignments.bam
samtools index $TEMP_DIR/${sample_name}_spliced_alignments.bam

samtools view -@ $num_threads -h $TEMP_DIR/${sample_name}_high_quality_cells_sense_transcriptome.bam | awk '$0 ~ /^@/ || $6 !~ /N/' | samtools view -@ $num_threads -b > $TEMP_DIR/${sample_name}_contiguous_alignments.bam
samtools index $TEMP_DIR/${sample_name}_contiguous_alignments.bam

##########################
# exon-exon junction reads
##########################

# we first filtered out the exonic alignments from splice alignments

bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_alignments.bam \
-b $fcb_dir/exons.bed \
-wa \
-f 1 \
-s \
-v \
> $TEMP_DIR/${sample_name}_spliced_alignments_not_in_exons.bam

# we then filter the alignments who across an exon start site
bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_alignments_not_in_exons.bam \
-b $fcb_dir/exons_start_site.bed \
-wa \
-F 1 \
-s \
> $TEMP_DIR/${sample_name}_spliced_alignments_not_in_exons_span_start_site.bam

# we then filter the alignments who across an exon end site
bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_alignments_not_in_exons_span_start_site.bam \
-b $fcb_dir/exons_end_site.bed \
-wa \
-F 1 \
-s \
> $TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bam 

samtools index $TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bam

n_header=$(samtools view -H $TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bam | wc -l)

## we then make sure the start and end position of each read are both in exons 
### we first convert BAM to BED
bedtools bamtobed \
-i $TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bam \
> $TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bed

### we convert the BED file to represent only the first base of the read
awk -v n_header=$n_header '{print $1"\t"$2"\t"$2+1"\t"NR+n_header"\t"$5"\t"$6}' \
$TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bed \
> $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_first_base.bed

### we convert the BED file to represent only the last base of the read
awk -v n_header=$n_header '{print $1"\t"$3-1"\t"$3"\t"NR+n_header"\t"$5"\t"$6}' \
$TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bed \
> $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_last_base.bed

### we find the alignments whose start site is in an exon
bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_first_base.bed \
-b $fcb_dir/exons.bed \
-wa \
-f 1 \
-s | \
awk '{print $4}' | uniq > $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_first_base_in_exons_row_id.txt

### we find the alignments whose end site is in an exon
bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_last_base.bed \
-b $fcb_dir/exons.bed \
-wa \
-f 1 \
-s \
| awk '{print $4}' | uniq  > $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_last_base_in_exons_row_id.txt

### find the common row ids
comm -12 $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_first_base_in_exons_row_id.txt $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_last_base_in_exons_row_id.txt > $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_row_id.txt

### add header rowids to the file
seq 1 $n_header > $TEMP_DIR/header_row_id.txt

cat $TEMP_DIR/header_row_id.txt $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_row_id.txt > $TEMP_DIR/${sample_name}_spliced_alignments_span_exons_row_id_with_header_row_id.txt

### we filter the BAM file
awk 'NR==FNR{a[$1]; next} FNR in a' \
$TEMP_DIR/${sample_name}_spliced_alignments_span_exons_row_id_with_header_row_id.txt \
<(samtools view -h $TEMP_DIR/${sample_name}_spliced_alignments_span_exons.bam) \
| samtools view -@ $num_threads -b > $OUT_DIR/${sample_name}_exon_exon_junctions.bam

samtools index $OUT_DIR/${sample_name}_exon_exon_junctions.bam

############################
# intron exon junction reads
############################

# we first filtered the alignment who is across an 
bedtools intersect \
-a $TEMP_DIR/${sample_name}_contiguous_alignments.bam \
-b $fcb_dir/exon_intron_junctions.bed \
-wa \
-F 1 \
-s \
> $OUT_DIR/${sample_name}_intron_exon_junctions.bam

samtools index $OUT_DIR/${sample_name}_intron_exon_junctions.bam

########################
# 3' terminal exon reads
########################
bedtools intersect \
-a $TEMP_DIR/${sample_name}_contiguous_alignments.bam \
-b $fcb_dir/transcripts_terminal_exon.bed \
-wa \
-f 1 \
-s \
> $OUT_DIR/${sample_name}_transcripts_terminal_exon.bam

samtools index $OUT_DIR/${sample_name}_transcripts_terminal_exon.bam

#############################################
# 3' terminal kilobase of spliced transcripts
#############################################

## we first filter reads in the contiguous regions from the 1000th base to the last base of spliced transcripts in genome
bedtools intersect \
-a $TEMP_DIR/${sample_name}_high_quality_cells_sense_transcriptome.bam \
-b $fcb_dir/spliced_transcripts_terminal_kilobase_exons_range.bed \
-wa \
-f 1 \
-s \
> $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bam

n_header=$(samtools view -H $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bam | wc -l)

## we then make sure the start and end position of each read are both in exons 
### we first convert BAM to BED
bedtools bamtobed \
-i $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bam \
> $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bed

### we convert the BED file to represent only the first base of the read
awk -v n_header=$n_header '{print $1"\t"$2"\t"$2+1"\t"NR+n_header"\t"$5"\t"$6}' \
$TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bed \
> $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_first_base.bed

### we convert the BED file to represent only the last base of the read
awk -v n_header=$n_header '{print $1"\t"$3-1"\t"$3"\t"NR+n_header"\t"$5"\t"$6}' \
$TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bed \
> $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_last_base.bed
### we find the alignments whose start site is in an exon
bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_first_base.bed \
-b $fcb_dir/spliced_transcripts_terminal_kilobase_exons.bed \
-wa \
-f 1 \
-s | \
awk '{print $4}' | uniq > $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_first_base_in_exons_row_id.txt

### we find the alignments whose end site is in an exon
bedtools intersect \
-a $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_last_base.bed \
-b $fcb_dir/spliced_transcripts_terminal_kilobase_exons.bed \
-wa \
-f 1 \
-s \
| awk '{print $4}' | uniq  > $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_last_base_in_exons_row_id.txt

### find the common row ids
comm -12 $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_first_base_in_exons_row_id.txt $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range_last_base_in_exons_row_id.txt > $TEMP_DIR/500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id.txt

### add header rowids to the file
seq 1 $n_header > $TEMP_DIR/header_row_id.txt

cat $TEMP_DIR/header_row_id.txt $TEMP_DIR/500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id.txt > $TEMP_DIR/500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id_with_header_row_id.txt

### we filter the BAM file
awk 'NR==FNR{a[$1]; next} FNR in a' \
$TEMP_DIR/500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id_with_header_row_id.txt \
<(samtools view -h $TEMP_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_range.bam) \
| samtools view -@ $num_threads -b > $OUT_DIR/${sample_name}_spliced_transcripts_terminal_kilobase.bam

samtools index $OUT_DIR/${sample_name}_spliced_transcripts_terminal_kilobase.bam

###############################################################
# exons in spliced transcripts' 3' terminal kilobase
###############################################################
bedtools intersect \
-a $TEMP_DIR/${sample_name}_contiguous_alignments.bam \
-b $fcb_dir/spliced_transcripts_terminal_kilobase_exons.bed \
-wa \
-f 1 \
-s \
> $OUT_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_exonic.bam

samtools index $OUT_DIR/${sample_name}_spliced_transcripts_terminal_kilobase_exonic.bam

###############################################
# 3' terminal kilobase of unspliced transcripts
###############################################
bedtools intersect \
-a $TEMP_DIR/${sample_name}_contiguous_alignments.bam \
-b $fcb_dir/unspliced_transcripts_terminal_kilobase.bed \
-wa \
-f 1 \
-s \
> $OUT_DIR/${sample_name}_unspliced_transcripts_terminal_kilobase.bam

samtools index $OUT_DIR/${sample_name}_unspliced_transcripts_terminal_kilobase.bam

########
# exonic
########
bedtools intersect \
-a $TEMP_DIR/${sample_name}_contiguous_alignments.bam \
-b $fcb_dir/exons.bed \
-wa \
-f 0.95 \
-s \
> $OUT_DIR/${sample_name}_exons.bam

samtools index $OUT_DIR/${sample_name}_exons.bam

#########
# introns
#########
bedtools intersect \
-a $TEMP_DIR/${sample_name}_contiguous_alignments.bam \
-b $fcb_dir/introns.bed \
-wa \
-f 0.95 \
-s > $OUT_DIR/${sample_name}_introns.bam

samtools index $OUT_DIR/${sample_name}_introns.bam


