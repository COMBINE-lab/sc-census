workflow ocr_count {
    take:
        txome_sense_bam
        high_quality_cells_bam
        sense_txome_umis

    main:
        // load peak sheet
        // species,sample_type,peak_file_path
        peaks = Channel
        .fromPath(params.input_files.peak_sheet)
        .splitCsv(header:true, sep: ",", strip: true)
        .map{ row-> tuple(row.species,
                            row.peak_name,
                            row.sample_type,
                            "${projectDir}/${row.peak_file_path}")
        }
        
        // get sense txome UMIs
        // out: species, sample_type, ref_name, sample_name, sense_txome_umis
        // make input
        // combine by species, sample_type
        input = high_quality_cells_bam
                    .combine(sense_txome_umis, by:[0,1,2,3,4])
                    // out: species, ref_name, sample_type, sample_name, fcb_dir, high_quality_cells_bam, sense_txome_umis
                    .combine(peaks, by:[0,2])
                    // out: species, sample_type, ref_name, sample_name, fcb_dir, high_quality_cells_bam, sense_txome_umis,peak_name, peak_file_path

        // #########
        // ocr count
        // #########
        ocr_count_intergenic(input)
        ocr_count_non_coding(input)
        ocr_count_not_sense_coding(input)
        ocr_count_all(input)

        emit:
            intergenic = ocr_count_intergenic.out
            non_coding = ocr_count_non_coding.out
            not_sense_coding = ocr_count_not_sense_coding.out
            all = ocr_count_all.out
}

process ocr_count_intergenic {
    label "cmd"
    publishDir "${params.output_dir}/quantification/${sample_name}/${ref_name}", mode: 'symlink'
    afterScript "rm -rf read_peak_distance.tsv bam.bam peaks.bed"

    input:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(sample_name),
            val(cells_or_nuclei),
            path(filtered_mtx_dir),
            path(fcb_dir),
            path(high_quality_cells_bam),
            path(sense_txome_umis),
            val(peak_name),
            path(peak_file_path)

    output:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(peak_name),
            val(sample_name),
            path("ocr_count_intergenic")

    """
    # 1. get peaks
    bedtools intersect \
    -a ${peak_file_path} \
    -b ${fcb_dir}/genes.bed \
    -wa -v > peaks.bed

    ## 2. get alingments
    bedtools intersect \
    -a ${high_quality_cells_bam} \
    -b ${fcb_dir}/genes.bed \
    -wa -v > bam.bam

    # 3. find non-coding alignments' closest non-coding peak
    bedtools closest \
    -d -t first \
    -a bam.bam \
    -b peaks.bed | awk '\$NF ~ /^[0-9]+\$/ && \$NF < ${params.read_peak_dist_threshold} && \$NF >= 0 {print \$4"\\t"\$13":"\$14"-"\$15"\\t"\$NF}' - > read_peak_distance.tsv

    cp $moduleDir/generate_count_matrix.py .
    python generate_count_matrix.py \
    bam.bam \
    read_peak_distance.tsv \
    $sense_txome_umis \
    $filtered_mtx_dir \
    ocr_count_intergenic

    """
}

process ocr_count_non_coding {
    label "cmd"
    publishDir "${params.output_dir}/quantification/${sample_name}/${ref_name}", mode: 'symlink'
    afterScript "rm -rf read_peak_distance.tsv bam.bam peaks.bed"

    input:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(sample_name),
            val(cells_or_nuclei),
            path(filtered_mtx_dir),
            path(fcb_dir),
            path(high_quality_cells_bam),
            path(sense_txome_umis),
            val(peak_name),
            path(peak_file_path)

    output:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(peak_name),
            val(sample_name),
            path("ocr_count_non_coding")

    """
    # 1. get peaks
    bedtools intersect \
    -a ${peak_file_path} \
    -b ${fcb_dir}/protein_coding_genes.bed \
    -wa -v > peaks.bed

    ## 2. get alingments
    bedtools intersect \
    -a ${high_quality_cells_bam} \
    -b ${fcb_dir}/protein_coding_genes.bed \
    -wa -v > bam.bam

    # 3. find non-coding alignments' closest non-coding peak
    bedtools closest \
    -d -t first \
    -a bam.bam \
    -b peaks.bed | awk '\$NF ~ /^[0-9]+\$/ && \$NF < ${params.read_peak_dist_threshold} && \$NF >= 0 {print \$4"\\t"\$13":"\$14"-"\$15"\\t"\$NF}' - > read_peak_distance.tsv


    cp $moduleDir/generate_count_matrix.py .
    python generate_count_matrix.py \
    bam.bam \
    read_peak_distance.tsv \
    $sense_txome_umis \
    $filtered_mtx_dir \
    ocr_count_non_coding
    """
}


process ocr_count_not_sense_coding {
    label "cmd"
    publishDir "${params.output_dir}/quantification/${sample_name}/${ref_name}", mode: 'symlink'
    afterScript "rm -rf read_peak_distance.tsv bam.bam"
    input:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(sample_name),
            val(cells_or_nuclei),
            path(filtered_mtx_dir),
            path(fcb_dir),
            path(high_quality_cells_bam),
            path(sense_txome_umis),
            val(peak_name),
            path(peak_file_path)

    output:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(peak_name),
            val(sample_name),
            path("ocr_count_not_sense_coding")

    """

    ## 2. get alingments
    bedtools intersect \
    -a ${high_quality_cells_bam} \
    -b ${fcb_dir}/protein_coding_genes.bed \
    -wa -v -s > bam.bam

    # 3. find non-coding alignments' closest non-coding peak
    bedtools closest \
    -d -t first \
    -a bam.bam \
    -b $peak_file_path | awk '\$NF ~ /^[0-9]+\$/ && \$NF < ${params.read_peak_dist_threshold} && \$NF >= 0 {print \$4"\\t"\$13":"\$14"-"\$15"\\t"\$NF}' - > read_peak_distance.tsv

    cp $moduleDir/generate_count_matrix.py .
    python generate_count_matrix.py \
    bam.bam \
    read_peak_distance.tsv \
    $sense_txome_umis \
    $filtered_mtx_dir \
    ocr_count_not_sense_coding
    """
}

process ocr_count_all {
    label "cmd"
    publishDir "${params.output_dir}/quantification/${sample_name}/${ref_name}", mode: 'symlink'
    afterScript "rm -rf read_peak_distance.tsv"
    input:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(sample_name),
            val(cells_or_nuclei),
            path(filtered_mtx_dir),
            path(fcb_dir),
            path(high_quality_cells_bam),
            path(sense_txome_umis),
            val(peak_name),
            path(peak_file_path)

    output:
        tuple val(species),
            val(sample_type),
            val(ref_name),
            val(peak_name),
            val(sample_name),
            path("ocr_count_all")

    """
    # 3. find non-coding alignments' closest non-coding peak
    bedtools closest \
    -d -t first \
    -a $high_quality_cells_bam \
    -b $peak_file_path | awk '\$NF ~ /^[0-9]+\$/ && \$NF < ${params.read_peak_dist_threshold} && \$NF >= 0 {print \$4"\\t"\$13":"\$14"-"\$15"\\t"\$NF}' - > read_peak_distance.tsv

    cp $moduleDir/generate_count_matrix.py .
    python generate_count_matrix.py \
    $high_quality_cells_bam \
    read_peak_distance.tsv \
    $sense_txome_umis \
    $filtered_mtx_dir \
    ocr_count_all \
    --intragenic
    """
}
