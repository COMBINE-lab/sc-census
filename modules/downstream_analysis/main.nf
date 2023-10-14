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

process celltypist() {
    label "py"
    input:
        tuple val(species),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei)
            val(orientation),
            path(quant_dir)
    output:
        val(species),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei)
            val(orientation),
            path(quant_dir),
            path("celltypist_predictions.csv")

    """
    copy $moduleDir/run_celltypist.py .
    python run_celltypist.py \
    $sample_type \
    $cells_or_nuclei \
    $quant_dir \
    .
    """
}

process clustering_analysis() {
    label "r"
    input:
        tuple val(species),
            val(ref_name),
            val(sample_type),
            val(sample_name),
            val(cells_or_nuclei)
            val(orientation),
            path(quant_dir),
            path(celltypist_predictions)
    
    output:
        path("clustering_analysis")

    """

    cp $moduleDir/clustering_analysis.rmd .
    cp $moduleDir/r_utils.rmd .
    mkdir clustering_analysis

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "clustering_analysis.rmd", 
        output_format = "html_document", 
        output_file = "${ref_name}/clustering_analysis.html",
        clean = TRUE,
        params = list(
            celltypist_csv = "${celltypist_predictions}",
            quant_dir = "${quant_dir}", 
            out_dir = "${clustering_analysis}",
            sample_type= "${sample_type}",
            num_threads= ${params.num_threads}
            )
        )
    EOF

    Rscript run_r.R
    """
}


