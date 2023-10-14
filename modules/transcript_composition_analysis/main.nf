workflow transcript_composition_analysis {
  // read in the reference sheet
    refs = Channel
      .fromPath(params.input_files.ref_sheet)
      .splitCsv(header:true, sep: ",", strip: true)
      .map{ row-> tuple(row.species,
                        row.ref_name,
                        row.genome_path,
                        row.gtf_path)
      }
    r_transcript_composition_analysis(refs)

    emit:
        r_transcript_composition_analysis.out.species_refname_genome_gtf_fcb
}

process r_transcript_composition_analysis {
    label 'r'
    publishDir "${params.output_dir}/transcript_composition_analysis", mode: 'copy', pattern: "${ref_name}"
    input:
        tuple val(species),
            val(ref_name),
            path(genome_path),
            path(gtf_path)
    output:
        path "${ref_name}"
        tuple val(species),
            val(ref_name),
            path(genome_path),
            path(gtf_path),
            path("${ref_name}/feature_category_bed"), emit: species_refname_genome_gtf_fcb

    script:
    p = params.transcript_composition_analysis
    
    """

    cp $moduleDir/transcript_composition_analysis.rmd .
    mkdir ${ref_name}

    cat << EOF > run_r.R
    rmarkdown::render(
        input = "transcript_composition_analysis.rmd", 
        output_format = "html_document", 
        output_file = "${ref_name}/transcript_composition_analysis.html",
        clean = TRUE,
        params = list(
            gtf_path = "${gtf_path}",
            genome_path = "${genome_path}", 
            out_dir = "${ref_name}",
            terminal_length= ${p.terminal_length},
            threads= ${params.num_threads},
            min_length= ${p.min_length},
            max_mismatch= ${p.max_mismatch}
            )
        )
    EOF

    Rscript run_r.R
    """
}

