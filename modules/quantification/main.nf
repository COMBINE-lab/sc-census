include {ocr_count} from "./ocr"
include {simpleaf_index} from "./simpleaf"
include {simpleaf_quant as simpleaf_quant_fw} from "./simpleaf"
include {simpleaf_quant as simpleaf_quant_rc} from "./simpleaf"


workflow quantification {
    take:
        txome_sense_bam
        txome_antisense_bam
        high_quality_cells_bam
        sense_txome_umis
    
    main:
    ocr_count(txome_sense_bam, high_quality_cells_bam, sense_txome_umis)
    simpleaf(txome_sense_bam, txome_antisense_bam)    

    
    emit:
        ocr = ocr_count.out
        simpleaf = simpleaf.out
}

workflow simpleaf {
    take:
        txome_sense_bam
        txome_antisense_bam

    main:
        // read in the reference sheet
        refs = Channel
        .fromPath(params.input_files.ref_sheet)
        .splitCsv(header:true, sep: ",", strip: true)
        .map{ row-> tuple(row.species,
                          row.ref_name,
                          "${projectDir}/${row.genome_path}",
                          "${projectDir}/${row.gtf_path}"
                        )
        }
        // ##############
        // simpleaf_index
        // ##############
        simpleaf_index(refs)
        // out: species, ref_name, index_dir

        // ##############
        // simpleaf_quant
        // ##############
        // txome sense quant
        
        fw_input = txome_sense_bam.combine(simpleaf_index.out, by:[0,1])
        simpleaf_quant_fw(fw_input, "fw")

        rc_input = txome_antisense_bam.combine(simpleaf_index.out, by:[0,1])
        simpleaf_quant_rc(rc_input, "rc")

        out = simpleaf_quant_fw.out.combine(simpleaf_quant_rc.out, by:[0,1,2,3,4])

        emit:
            out
        
}
