process simpleaf_index() {
    label "cmd"
    input:
        tuple val(species),
            val(ref_name),
            path(genome_path),
            path(gtf_path)
    output:
        tuple val(species),
            val(ref_name),
            path("${ref_name}_simpleaf_index/index")
    
    """
    export ALEVIN_FRY_HOME="af_home"
    # ulimit -n 2048
    
    simpleaf set-paths

    simpleaf index \
        --output ${ref_name}_simpleaf_index \
        --threads $params.num_threads \
        --ref-type spliceu \
        --fasta ${genome_path} \
        --gtf ${gtf_path}

    cp ${ref_name}_simpleaf_index/ref/gene_id_to_name.tsv ${ref_name}_simpleaf_index/index
    """
}

process simpleaf_quant() {
    label "cmd"
    publishDir "${params.output_dir}/quantification/${sample_name}/${ref_name}/simpleaf_quant/${orientation}", mode: 'symlink'
    input:
        tuple val(species),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fcb_dir),
            path(in_bam),
            path(index_dir)
        val orientation
    output:
        tuple val(species),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            val(orientation),
            path("simpleaf_quant/af_quant/*")
    
    """
    export ALEVIN_FRY_HOME="af_home"
    simpleaf set-paths

    #####################################
    # quantify sense transcriptomic reads
    #####################################
    bamtofastq \
    --nthreads=${params.num_threads} \
    $in_bam \
    fastqs

    reads1="\$(find -L fastqs -name "*_R1_*" -type f | sort | awk -v OFS=, '{\$1=\$1;print}' | paste -sd, -)"
    reads2="\$(find -L fastqs -name "*_R2_*" -type f | sort | awk -v OFS=, '{\$1=\$1;print}' | paste -sd, -)"

    simpleaf quant \
        --output simpleaf_quant \
        --threads ${params.num_threads} \
        --expected-ori ${orientation} \
        --chemistry 10xv3 \
        --index ${index_dir} \
        --reads1 \$reads1 \
        --reads2 \$reads2 \
        --unfiltered-pl \
        --resolution cr-like

    cp ${index_dir}/gene_id_to_name.tsv simpleaf_quant/af_quant
    """
}

