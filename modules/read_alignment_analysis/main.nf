
// at this moment, we have species, sample_name, bam_file, filtered_mtx_dir, ref_name, genome, gtf, fcb_dir

process high_quality_cells_bam {
    label "cmd"
    afterScript 'rm -rf barcodes.txt'

    input:
        tuple val(species),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(in_bam),
            path(filtered_mtx_dir),
            val(ref_name),
            val(genome),
            val(gtf),
            path(fcb_dir)
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            path(fcb_dir),
            path("high_quality_cells.bam")
    script:
    // run script
    """
    gunzip -c $filtered_mtx_dir/barcodes.tsv.gz > barcodes.txt

    # filter only high quality cell's reads and mapped reads
    samtools view -@ ${params.num_threads} -bh -F 4 -D CB:barcodes.txt $in_bam | samtools view -@ ${params.num_threads} -bh -d UB -o high_quality_cells.bam -
    """
}

process txome_sense_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(high_quality_cells_bam)
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path("txome_sense.bam")
    
    script:
    // run script
    """
    bedtools intersect \
    -a $high_quality_cells_bam \
    -b $fcb_dir/genes.bed \
    -wa \
    -f 1 \
    -s \
    > txome_sense.bam
    """
}

process txome_antisense_bam {
    label "cmd"
    afterScript "rm txome_antisense_all.bam"

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(high_quality_cells_bam),
            path(sense_txome_umis)

    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path("txome_antisense.bam")
    
    script:
    // run script
    """
    bedtools intersect \
    -a $high_quality_cells_bam \
    -b $fcb_dir/genes.bed \
    -wa \
    -f 1 \
    -S \
    > txome_antisense_all.bam

    cp $moduleDir/filter_antisense_reads.py .
    python filter_antisense_reads.py \
    txome_antisense_all.bam \
    $sense_txome_umis \
    txome_antisense.bam
    """
}

process txome_splice_bam {
    label "cmd"

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(high_quality_cells_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path("txome_splice.bam")
    
    script:
    """
    samtools view -@ ${params.num_threads} -h $high_quality_cells_bam | awk '\$0 ~ /^@/ || \$6 ~ /N/' | samtools view -@ ${params.num_threads} -b > txome_splice.bam
    """
}

process txome_contiguous_bam {
    label "cmd"

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(high_quality_cells_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path("txome_contiguous.bam")
    
    script:
    """
    samtools view -@ ${params.num_threads} -h $high_quality_cells_bam | awk '\$0 ~ /^@/ || \$6 !~ /N/' | samtools view -@ ${params.num_threads} -b > txome_contiguous.bam
    """
}

process exon_exon_junction_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'
    afterScript 'rm -rf spliced_alignments_not_in_exons.bam spliced_alignments_not_in_exons_span_start_site.bam  spliced_alignments_span_exons.bam '


    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(txome_splice_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("exon_exon_junctions.bam")
    
    script:
    """
    # we first filtered out the exonic alignments from splice alignments
    bedtools intersect \
    -a $txome_splice_bam \
    -b $fcb_dir/exons.bed \
    -wa \
    -f 1 \
    -s \
    -v \
    > spliced_alignments_not_in_exons.bam

    # we then filter the alignments who across an exon start site
    bedtools intersect \
    -a spliced_alignments_not_in_exons.bam \
    -b $fcb_dir/exons_start_site.bed \
    -wa \
    -F 1 \
    -s \
    > spliced_alignments_not_in_exons_span_start_site.bam

    # we then filter the alignments who across an exon end site
    bedtools intersect \
    -a spliced_alignments_not_in_exons_span_start_site.bam \
    -b $fcb_dir/exons_end_site.bed \
    -wa \
    -F 1 \
    -s \
    > spliced_alignments_span_exons.bam 

    samtools index spliced_alignments_span_exons.bam

    n_header=\$(samtools view -H spliced_alignments_span_exons.bam | wc -l)

    ## we then make sure the start and end position of each read are both in exons 
    ### we first convert BAM to BED
    bedtools bamtobed \
    -i spliced_alignments_span_exons.bam \
    > spliced_alignments_span_exons.bed

    ### we convert the BED file to represent only the first base of the read
    awk -v n_header=\$n_header '{print \$1"\\t"\$2"\\t"\$2+1"\\t"NR+n_header"\\t"\$5"\\t"\$6}' \
    spliced_alignments_span_exons.bed \
    > spliced_alignments_span_exons_first_base.bed

    ### we convert the BED file to represent only the last base of the read
    awk -v n_header=\$n_header '{print \$1"\\t"\$3-1"\\t"\$3"\\t"NR+n_header"\\t"\$5"\\t"\$6}' \
    spliced_alignments_span_exons.bed \
    > spliced_alignments_span_exons_last_base.bed

    ### we find the alignments whose start site is in an exon
    bedtools intersect \
    -a spliced_alignments_span_exons_first_base.bed \
    -b $fcb_dir/exons.bed \
    -wa \
    -f 1 \
    -s | \
    awk '{print \$4}' | uniq > spliced_alignments_span_exons_first_base_in_exons_row_id.txt

    ### we find the alignments whose end site is in an exon
    bedtools intersect \
    -a spliced_alignments_span_exons_last_base.bed \
    -b $fcb_dir/exons.bed \
    -wa \
    -f 1 \
    -s \
    | awk '{print \$4}' | uniq  > spliced_alignments_span_exons_last_base_in_exons_row_id.txt

    ### find the common row ids
    comm --nocheck-order -12 --nocheck-order spliced_alignments_span_exons_first_base_in_exons_row_id.txt spliced_alignments_span_exons_last_base_in_exons_row_id.txt > spliced_alignments_span_exons_row_id.txt

    ### add header rowids to the file
    seq 1 \$n_header > header_row_id.txt

    cat header_row_id.txt spliced_alignments_span_exons_row_id.txt > spliced_alignments_span_exons_row_id_with_header_row_id.txt

    ### we filter the BAM file
    awk 'NR==FNR{a[\$1]; next} FNR in a' \
    spliced_alignments_span_exons_row_id_with_header_row_id.txt \
    <(samtools view -h spliced_alignments_span_exons.bam) \
    | samtools view -@ ${params.num_threads} -b > exon_exon_junctions.bam
    """
}

process intron_exon_junction_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(contiguous_alignments_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("intron_exon_junctions.bam")
    
    script:
    """
    # we first filtered the alignment who is across an 
    bedtools intersect \
    -a $contiguous_alignments_bam \
    -b $fcb_dir/exon_intron_junctions.bed \
    -wa \
    -F 1 \
    -s \
    > intron_exon_junctions.bam
    """
}

process txp_terminal_exon_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'
    
    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(contiguous_alignments_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("transcripts_terminal_exons.bam")
    
    script:
    """
    bedtools intersect \
    -a $contiguous_alignments_bam \
    -b $fcb_dir/transcripts_terminal_exon.bed \
    -wa \
    -f 1 \
    -s \
    > transcripts_terminal_exons.bam
    """
}


process spliced_txp_terminal_kilobase_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'
    afterScript 'rm -rf spliced_transcripts_terminal_kilobase_range.bam '

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(txome_sense_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path("spliced_txp_three_prime_terminal_kilobase.bam")
    
    script:
    """
    ## we first filter reads in the contiguous regions from the 1000th base to the last base of spliced transcripts in genome
    bedtools intersect \
    -a $txome_sense_bam \
    -b $fcb_dir/spliced_transcripts_terminal_kilobase_exons_range.bed \
    -wa \
    -f 1 \
    -s \
    > spliced_transcripts_terminal_kilobase_range.bam

    n_header=\$(samtools view -H spliced_transcripts_terminal_kilobase_range.bam | wc -l)

    ## we then make sure the start and end position of each read are both in exons 
    ### we first convert BAM to BED
    bedtools bamtobed \
    -i spliced_transcripts_terminal_kilobase_range.bam \
    > spliced_transcripts_terminal_kilobase_range.bed

    ### we convert the BED file to represent only the first base of the read
    awk -v n_header=\$n_header '{print \$1"\\t"\$2"\\t"\$2+1"\\t"NR+n_header"\\t"\$5"\\t"\$6}' \
    spliced_transcripts_terminal_kilobase_range.bed \
    > spliced_transcripts_terminal_kilobase_range_first_base.bed

    ### we convert the BED file to represent only the last base of the read
    awk -v n_header=\$n_header '{print \$1"\\t"\$3-1"\\t"\$3"\\t"NR+n_header"\\t"\$5"\\t"\$6}' \
    spliced_transcripts_terminal_kilobase_range.bed \
    > spliced_transcripts_terminal_kilobase_range_last_base.bed
    ### we find the alignments whose start site is in an exon
    bedtools intersect \
    -a spliced_transcripts_terminal_kilobase_range_first_base.bed \
    -b $fcb_dir/spliced_transcripts_terminal_kilobase_exons.bed \
    -wa \
    -f 1 \
    -s | \
    awk '{print \$4}' | uniq > spliced_transcripts_terminal_kilobase_range_first_base_in_exons_row_id.txt

    ### we find the alignments whose end site is in an exon
    bedtools intersect \
    -a spliced_transcripts_terminal_kilobase_range_last_base.bed \
    -b $fcb_dir/spliced_transcripts_terminal_kilobase_exons.bed \
    -wa \
    -f 1 \
    -s \
    | awk '{print \$4}' | uniq  > spliced_transcripts_terminal_kilobase_range_last_base_in_exons_row_id.txt

    ### find the common row ids
    comm -12 --nocheck-order spliced_transcripts_terminal_kilobase_range_first_base_in_exons_row_id.txt spliced_transcripts_terminal_kilobase_range_last_base_in_exons_row_id.txt > 500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id.txt

    ### add header rowids to the file
    seq 1 \$n_header > header_row_id.txt

    cat header_row_id.txt 500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id.txt > 500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id_with_header_row_id.txt

    ### we filter the BAM file
    awk 'NR==FNR{a[\$1]; next} FNR in a' \
    500_PBMC_3p_LT_Chromium_X_possorted_genome_spliced_transcripts_terminal_kilobase_alignments_row_id_with_header_row_id.txt \
    <(samtools view -h spliced_transcripts_terminal_kilobase_range.bam) \
    | samtools view -@ ${params.num_threads} -b > spliced_txp_three_prime_terminal_kilobase.bam
    """
}


process spliced_txp_terminal_kilobase_exonic_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(txome_sense_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("spliced_txp_three_prime_terminal_kilobase_exonic.bam")
    
    script:
    """
    bedtools intersect \
    -a $txome_sense_bam \
    -b $fcb_dir/spliced_transcripts_terminal_kilobase_exons.bed \
    -wa \
    -f 1 \
    -s \
    > spliced_txp_three_prime_terminal_kilobase_exonic.bam
    """
}


process unspliced_txp_terminal_kilobase_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(contiguous_alignments_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("unspliced_txp_three_prime_terminal_kilobase.bam")
    
    script:
    """
    bedtools intersect \
    -a $contiguous_alignments_bam \
    -b $fcb_dir/unspliced_transcripts_terminal_kilobase.bed \
    -wa \
    -f 1 \
    -s \
    > unspliced_txp_three_prime_terminal_kilobase.bam
    """
}

process exon_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(contiguous_alignments_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("exons.bam")
    
    script:
    """
    bedtools intersect \
    -a $contiguous_alignments_bam \
    -b $fcb_dir/exons.bed \
    -wa \
    -f 1 \
    -s \
    > exons.bam
    """
}


process introns_bam {
    label "cmd"
    publishDir "${params.output_dir}/read_alignment_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(contiguous_alignments_bam)
    
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("introns.bam")
    
    script:
    """
    bedtools intersect \
    -a $contiguous_alignments_bam \
    -b $fcb_dir/introns.bed \
    -wa \
    -f 1 \
    -s \
    > introns.bam
    """
}

// This process take the feature_category_bams as input
process umi_category_analysis {
    label "py"
    publishDir "${params.output_dir}/downstream_analysis/${sample_name}/${ref_name}", mode: 'symlink'
    input:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(exon_exon_junction_bam),
            path(intron_exon_junction_bam),
            path(transcripts_terminal_exon_bam),
            path(fcb_dir),
            path(spliced_txp_terminal_kilobase_bam),
            path(spliced_txp_terminal_kilobase_exonic_bam),
            path(unspliced_txp_terminal_kilobase_bam),
            path(exon_bam),
            path(introns_bam)
    output:
        tuple val(species), 
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("read_alignment_analysis/cb-umi_dict.pickle"), emit: sense_txome_umis
        path("read_alignment_analysis")
    """
    cp $moduleDir/umi_category_analysis.py .

    python umi_category_analysis.py \
    $species \
    $sample_type \
    $sample_name \
    read_alignment_analysis \
    $exon_exon_junction_bam \
    $intron_exon_junction_bam \
    $transcripts_terminal_exon_bam \
    $spliced_txp_terminal_kilobase_bam \
    $spliced_txp_terminal_kilobase_exonic_bam \
    $unspliced_txp_terminal_kilobase_bam \
    $exon_bam \
    $introns_bam
    """
}


process sense_txome_umis {
    label "cmd"
    input:
        tuple val(species), 
            val(sample_type),
            val(ref_name),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(txome_sense_bam)
    output:
        tuple val(species),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("sense_txome_umis.pkl")

    """
    cp $moduleDir/make_cb-umi_list.py .
    python make_cb-umi_list.py \
    $txome_sense_bam \
    sense_txome_umis.pkl
    """
}


workflow read_alignment_analysis {
    take: species_refname_genome_gtf_fcb

    main:
        // The sample sheet records the information of the RNA-seq samples 
        samples = Channel
            .fromPath(params.input_files.sample_sheet)
            .splitCsv(header:true, sep: ",", strip: true)
            .map{ row-> tuple(row.species,
                            row.sample_type,
                            row.sample_name,
                            row.cells_or_nuclei,
                            row.bam_file,
                            row.filtered_mtx_dir)
            }

        workflow_input = samples
            .combine(species_refname_genome_gtf_fcb, by: 0)
        // at this moment, we have species, ref_name, sample_type, sample_name, bam_file, filtered_mtx_dir, genome, gtf, fcb_dir
        // ##########################################
        // 1. we first filter the high quality cells
        // ##########################################
        high_quality_cells_bam(workflow_input)
        // out: species, ref_name, sample_type, sample_name, fcb_dir, high_quality_cells_bam

        // ##############################################
        // 2. split txome sense and antisense alignments
        // ##############################################
        txome_sense_bam(high_quality_cells_bam.out)
        // out: species, ref_name, sample_type, sample_name, fcb_dir, txome_sense_bam

        // ################################################
        // 3. split txome splice and contiguous alignments
        // ################################################
        txome_splice_bam(high_quality_cells_bam.out)
        // out: species, ref_name, sample_type, sample_name, fcb_dir, txome_splice_bam

        txome_contiguous_bam(high_quality_cells_bam.out)
        // out: species, ref_name, sample_type, sample_name, fcb_dir, txome_contiguous_bam

        // ################################################################
        // 4. split exon-exon junction and intron-exon junction alignments
        // ################################################################
        exon_exon_junction_bam(txome_splice_bam.out)
        // out: species, ref_name, sample_type, sample_name, exon_exon_junction_bam

        intron_exon_junction_bam(txome_contiguous_bam.out)
        // out: species, ref_name, sample_type, sample_name, intron_exon_junction_bam

        // ##############################################
        // 5. split transcripts terminal exon alignments
        // ##############################################
        txp_terminal_exon_bam(txome_contiguous_bam.out)
        // out: species, ref_name, sample_type, sample_name, transcripts_terminal_exon_bam

        // ##########################################################
        // 6. split spliced and unspliced transcripts terminal kilobase alignments
        // ##########################################################
        spliced_txp_terminal_kilobase_bam(txome_sense_bam.out)
        // out: species, ref_name, sample_type, sample_name, spliced_txp_three_prime_terminal_kilobase_bam

        spliced_txp_terminal_kilobase_exonic_bam(txome_sense_bam.out)
        // out: species, ref_name, sample_type, sample_name, spliced_txp_three_prime_terminal_kilobase_exonic_bam

        unspliced_txp_terminal_kilobase_bam(txome_contiguous_bam.out)
        // out: species, ref_name, sample_type, sample_name, unspliced_txp_three_prime_terminal_kilobase_bam

        // ####################################
        // 7. split exon and intron alignments
        // ####################################
        exon_bam(txome_contiguous_bam.out)
        // out: species, ref_name, sample_type, sample_name, exon_bam

        introns_bam(txome_contiguous_bam.out)
        // out: species, ref_name, sample_type, sample_name, introns_bam

        // ##################
        // 8. post-processing
        // ##################
        // To analyze feature categories' alignments, we need their BAM files
        feature_category_bams = exon_exon_junction_bam.out.combine(intron_exon_junction_bam.out, by: [0,1,2,3,4])
                .combine(txp_terminal_exon_bam.out, by: [0,1,2,3,4])
                .combine(spliced_txp_terminal_kilobase_bam.out, by: [0,1,2,3,4])
                .combine(spliced_txp_terminal_kilobase_exonic_bam.out, by: [0,1,2,3,4])
                .combine(unspliced_txp_terminal_kilobase_bam.out, by: [0,1,2,3,4])
                .combine(exon_bam.out, by: [0,1,2,3,4])
                .combine(introns_bam.out, by: [0,1,2,3,4])

        umi_category_analysis(feature_category_bams)

        // ##########################################
        // 9. txome antisense alignments
        // ##########################################
        txome_antisense_bam_input = high_quality_cells_bam.out.combine(umi_category_analysis.out.sense_txome_umis, by: [0,1,2,3,4])

        txome_antisense_bam(txome_antisense_bam_input)

        emit:
            txome_sense_bam = txome_sense_bam.out
            txome_antisense_bam = txome_antisense_bam.out
            high_quality_cells_bam = high_quality_cells_bam.out
            sense_txome_umis = umi_category_analysis.out.sense_txome_umis
}   






