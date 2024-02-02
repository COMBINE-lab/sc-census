// 1. read alignment analysis - Python
//      In this analysis, we want to know how many reads and UMIs fall into each category. An upset plot will be plotted for reads and UMIs separately for each dataset. This analysis was integrated into the read_alignment_analysis pipeline.
// 2. clustering analysis - R
//      In this analysis, we want to
//      1. show that when fixing the number of clusters as a small number, the clustering results of all count matrices and their combinations are similar
//      2. show that when increasing the number of clusters, the clustering results of count matrices and their combinations are diverging
//      3. Two UMAP/tSNE plots will be plotted for each final count matrix, one using the cell types from "standard" clustering, the other using its own clustering results 
// 3. Antisense as the imputation of sense count matrix - R
//     In this analysis, we want to show that the 0s in sense count matrix can be imputed by antisense count matrix
//     And some of the antisense only genes are important for each cell type
// 4. DE analysis - R
//     In this analysis, we want to show that 
//      1. the DE genes are not the same but similar between spliced, unspliced and ambiguous count matrices
//      2. the DE genes are not the same but similar between sense and antisense count matrices
//      3. the DEGs in those count matrices are biologically meaningful
// 5. open chromatin region-associated count matrix - R
//      In this analysis, we want to show that 
//      1. the ARI of the clustering results using RNA peak count matrices are high
//      2. There are some existing explainations, but we want to come up with a new one, which is enhancer RNAs and other non-coding RNAs
//      3. show that the similarity between the same cells' RNA peak count and ATAC count matrix is significantly higher than other cells
//      4. show that when using bipartite graph to connect RNA peak and ATAC peaks, we can cluster the bipartite graph like what Rob suggested

workflow downstream_analysis {
    take:
        ocr
        simpleaf

    main:
        sctype(simpleaf)
        simpleaf_analysis_input = simpleaf.combine(sctype.out.out, by : [0,1,2,3,4])
        clustering_analysis(simpleaf_analysis_input)
        antisense_count_analysis(simpleaf_analysis_input)

        ocr_count_analysis_input = ocr.combine(sctype.out.out, by : [0,1,2,3,4])
        ocr_count_analysis(ocr_count_analysis_input)
        atac_ocr_correlation(ocr_count_analysis_input)
        
}

process celltypist() {
    publishDir "${params.output_dir}/downstream_analysis/celltypist_predictions/${sample_name}/${ref_name}", mode: 'symlink'
    label "py"
    input:
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fw_af_quant),
            path(rc_af_quant)

    output:
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fw_af_quant),
            path("celltypist_out/celltypist_predictions.csv")

    """
    cp $moduleDir/run_celltypist.py .
    python run_celltypist.py \
    $sample_type \
    $cells_or_nuclei \
    $fw_af_quant \
    celltypist_out
    """
}


process sctype() {
    publishDir "${params.output_dir}/downstream_analysis/sctype_predictions/${sample_name}/${ref_name}", mode: 'symlink'
    label "r"
    input:
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fw_af_quant),
            path(rc_af_quant)

    output:
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path(fw_af_quant),
            path("sctype_out/sctype_predictions.csv"), emit: out
        path("sctype_out/run_sctype.html")
    """

    cp $moduleDir/run_sctype.rmd .
    cp $moduleDir/r_utils.R .
    mkdir sctype_out

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "run_sctype.rmd", 
        output_format = "html_document", 
        output_file = "sctype_out/run_sctype.html",
        clean = TRUE,
        params = list(
            utils_path = "r_utils.R",
            fw_quant_dir = "${fw_af_quant}", 
            out_dir = "sctype_out",
            cells_or_nuclei = "${cells_or_nuclei}",
            sample_type= "${sample_type}",
            num_threads= ${params.large_num_threads}
        )
    )
    EOF

    Rscript run_r.R
    
    """
}

process clustering_analysis() {
    label "r"
    publishDir "${params.output_dir}/downstream_analysis/clustering_analysis/${sample_name}/${ref_name}", mode: 'symlink'

    input:
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("fw"),
            path("rc"),
            path("nevermind"),
            path(sctype_predictions)
    
    output:
        path("clustering_analysis")

    script:

    """

    cp $moduleDir/clustering_analysis.rmd .
    cp $moduleDir/r_utils.R .
    mkdir clustering_analysis

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "clustering_analysis.rmd", 
        output_format = "html_document", 
        output_file = "clustering_analysis/clustering_analysis.html",
        clean = TRUE,
        params = list(
            utils_path = "r_utils.R",
            sctype_csv = "${sctype_predictions}",
            fw_quant_dir = "fw", 
            out_dir = "clustering_analysis",
            sample_type = "${sample_type}",
            cells_or_nuclei = "${cells_or_nuclei}",
            num_threads= ${params.large_num_threads}
        )
    )
    EOF

    Rscript run_r.R
    exit 0
    """
}


process antisense_count_analysis() {
    label "r"
    publishDir "${params.output_dir}/downstream_analysis/antisense_count_analysis/${sample_name}/${ref_name}", mode: 'symlink', pattern: "antisense_count_analysis"

    input:
        tuple val(species_for_rerun_three),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            path("fw"),
            path("rc"),
            path("nevermind"),
            path(sctype_predictions)
    
    output:
        path("antisense_count_analysis")

    """
    cp $moduleDir/antisense_count_analysis.rmd .
    cp $moduleDir/r_utils.R .
    mkdir -p antisense_count_analysis

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "antisense_count_analysis.rmd", 
        output_format = "html_document", 
        output_file = "antisense_count_analysis/antisense_count_analysis.html",
        clean = TRUE,
        params = list(
            utils_path = "r_utils.R",
            sctype_csv = "${sctype_predictions}",
            fw_quant_dir = "fw", 
            rc_quant_dir = "rc", 
            out_dir = "antisense_count_analysis",
            sample_type= "${sample_type}",
            num_threads= ${params.large_num_threads}
        )
    )
    quit()
    EOF

    Rscript run_r.R
    exit 0
    """
}

process ocr_count_analysis() {
    label "r"
    publishDir "${params.output_dir}/downstream_analysis/ocr_count_analysis/${sample_name}/${ref_name}/${peak_name}", mode: 'symlink', pattern: "ocr_count_analysis"

    input:
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            val(peak_name),
            path(intergenic_count),
            path(ocr_count_intergenic),
            path(ccres_count_intergenic),
            path(ocr_count_non_coding),
            path(ccres_count_non_coding),
            path(ocr_count_not_sense_coding),
            path(ccres_count_not_sense_coding),
            path(ocr_count_all),
            path(ccres_count_all),
            path(fw_quant_dir),
            path(sctype_predictions)
    
    output:
        path("ocr_count_analysis")

    """

    cp $moduleDir/ocr_count_analysis.rmd .
    cp $moduleDir/r_utils.R .
    mkdir -p ocr_count_analysis

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "ocr_count_analysis.rmd", 
        output_format = "html_document", 
        output_file = "ocr_count_analysis/ocr_count_analysis.html",
        clean = TRUE,
        params = list(
            utils_path = "r_utils.R",
            sctype_csv = "${sctype_predictions}",
            simpleaf_count_dir = "${fw_quant_dir}",
            intergenic_count_dir = "${intergenic_count}",
            ocr_count_intergenic_dir = "${ocr_count_intergenic}",
            ocr_count_non_coding_dir = "${ocr_count_non_coding}",
            ocr_count_not_sense_coding_dir = "${ocr_count_not_sense_coding}",
            ocr_count_all_dir = "${ocr_count_all}",
            ccres_count_intergenic_dir = "${ccres_count_intergenic}",
            ccres_count_non_coding_dir = "${ccres_count_non_coding}",
            ccres_count_not_sense_coding_dir = "${ccres_count_not_sense_coding}",
            ccres_count_all_dir = "${ccres_count_all}",
            out_dir = "ocr_count_analysis", 
            sample_type = "${cells_or_nuclei}",
            num_threads = ${params.large_num_threads},
            read_peak_dist_threshold = ${params.read_peak_dist_threshold}
        )
    )
    quit()
    EOF

    Rscript run_r.R
    exit 0
    """
}

process atac_ocr_correlation {
    label "r"
    publishDir "${params.output_dir}/downstream_analysis/atac_ocr_correlation/${sample_name}/${ref_name}", mode: 'symlink', pattern: "atac_ocr_correlation"

    input: 
        tuple val(species_for_rerun_two),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei),
            val(peak_name),
            path(intergenic_count),
            path(ocr_count_intergenic),
            path(ccres_count_intergenic),
            path(ocr_count_non_coding),
            path(ccres_count_non_coding),
            path(ocr_count_not_sense_coding),
            path(ccres_count_not_sense_coding),
            path(ocr_count_all),
            path(ccres_count_all),
            path(fw_quant_dir),
            path(sctype_predictions)

    output:
        path("atac_ocr_correlation")

    when:
    sample_name =~ /.*ultiome.*/ && peak_name == 'multiomics_peaks'
    
    script:

    atac_count_dir=file("${projectDir}/data/datasets/${sample_name}/filtered_feature_bc_matrix")
    """
    cp $moduleDir/atac_ocr_correlation.rmd .
    cp $moduleDir/r_utils.R .
    mkdir -p atac_ocr_correlation

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "atac_ocr_correlation.rmd", 
        output_format = "html_document", 
        output_file = "atac_ocr_correlation/atac_ocr_correlation.html",
        clean = TRUE, 
        params = list(
            utils_path = "r_utils.R",
            atac_count_dir = "${atac_count_dir}",
            ocr_count_intergenic_dir = "${ocr_count_intergenic}",
            ocr_count_non_coding_dir = "${ocr_count_non_coding}",
            ocr_count_not_sense_coding_dir = "${ocr_count_not_sense_coding}",
            ocr_count_all_dir = "${ocr_count_all}",
            ccres_count_intergenic_dir = "${ccres_count_intergenic}",
            ccres_count_non_coding_dir = "${ccres_count_non_coding}",
            ccres_count_not_sense_coding_dir = "${ccres_count_not_sense_coding}",
            ccres_count_all_dir = "${ccres_count_all}",
            out_dir = "atac_ocr_correlation", 
            sample_type = "${cells_or_nuclei}",
            num_threads = ${params.large_num_threads}
        )
    )
    quit()
    EOF

    Rscript run_r.R
    exit 0
    """

}

