// We have the following analyses:
// 1. transcript composition analysis
// 2. read alignment analysis
// 3. quantification
// 4. downstream analysis:
    // 1. count matrix processing
    // 2. differential expression analysis
    // 3. cluster similarity analysis

include {transcript_composition_analysis} from "./modules/transcript_composition_analysis"
include {read_alignment_analysis} from "./modules/read_alignment_analysis"
include {quantification} from "./modules/quantification"
include {downstream_analysis} from "./modules/downstream_analysis"

workflow {
  // ##################################
  // 1. transcript composition analysis
  // ##################################
  // process the reference spreadsheet and perform transcript composition analysis
  transcript_composition_analysis()
  // out: tuple[species, ref_name, genome_path, gtf_path, fcb_dir]

  // ##########################
  // 2. read alignment analysis
  // ##########################
  read_alignment_analysis(transcript_composition_analysis.out)

  // #################
  // 3. quantification
  // #################
  quantification(read_alignment_analysis.out)

  downstream_analysis(quantification.out)
  
}
