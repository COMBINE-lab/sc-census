
count_percentage <- function(p = ".") {
    cols = c('10k_PBMC_Multiome_nextgem_Chromium_X','10k_PBMC_3p_nextgem_Chromium_X','bmmc_site4_donor08_multiome','bmmc_site4_donor08_cite','human_brain_3k_multiome','e18_mouse_brain_fresh_5k_multiome','SC3_v3_NextGem_DI_Nuclei_5K_Multiplex')
    rows = c("Txome Count", "% Spliced Txome", "% Unspliced Txome", "% Ambiguous Txome",
        "Antisense Count", "% Spliced Antisense", "% Unspliced Antisense", "% Ambiguous Antisense",
        "Intergenic Count", 
        "cCREs Intergenic Count", "% cCREs Intergenic",  "cCREs Non-Coding Count", "cCREs All Count",
        "OCRs Intergenic Count", "% OCRs Intergenic", "OCRs Non-Coding Count", "OCR All Count",
        " OCRs Intergenic Count", " % OCRs Intergenic", " OCRs Non-Coding Count", " OCR All Count"
    )
    df = data.frame(matrix(NA, ncol = length(cols), nrow = length(rows)))
    rownames(df) = rows
    colnames(df) = cols
    df_list = list(ref_pool_df = df, ref_10x_df = df)

    for(dataset_dir in list.dirs(p, recursive = FALSE)) {
        dataset_name = basename(dataset_dir)
        for(ref_dir in list.dirs(dataset_dir, recursive = FALSE)) {
            df_name = if(grepl("pool_lab", basename(ref_dir))) {"ref_pool_df"} else {"ref_10x_df"}
            # sense     
            sce = fishpond::loadFry(file.path(ref_dir, "simpleaf_quant", "fw_af_quant"), outputFormat = "raw")
            s = sum(sce@assays@data$spliced)
            u = sum(sce@assays@data$unspliced)
            a = sum(sce@assays@data$ambiguous)
            total = s+u+a
            df_list[[df_name]]["Txome Count", dataset_name] = formatC(total, format="d", big.mark=",")
            df_list[[df_name]]["% Spliced Txome", dataset_name] = paste0(round(s/total*100, 2),"%")
            df_list[[df_name]]["% Unspliced Txome", dataset_name] = paste0(round(u/total*100, 2),"%")
            df_list[[df_name]]["% Ambiguous Txome", dataset_name] = paste0(round(a/total*100, 2),"%")

            # antisense
            sce = fishpond::loadFry(file.path(ref_dir, "simpleaf_quant", "rc_af_quant"), outputFormat = "raw")
            s = sum(sce@assays@data$spliced)
            u = sum(sce@assays@data$unspliced)
            a = sum(sce@assays@data$ambiguous)
            total = s+u+a
            df_list[[df_name]]["Antisense Count", dataset_name] = formatC(total, format="d", big.mark=",")
            df_list[[df_name]]["% Spliced Antisense", dataset_name] = paste0(round(s/total*100, 2),"%")
            df_list[[df_name]]["% Unspliced Antisense", dataset_name] = paste0(round(u/total*100, 2),"%")
            df_list[[df_name]]["% Ambiguous Antisense", dataset_name] = paste0(round(a/total*100, 2),"%")

            peak_name = "multiomics_peaks"
            intergenic_count = Seurat::Read10X(file.path(ref_dir, peak_name, "intergenic_count"))
            ocr_count_intergenic = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_intergenic"))
            ocr_count_non_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_non_coding"))
            ocr_count_not_sense_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_not_sense_coding"))
            ocr_count_all = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_all"))
            ccres_count_intergenic = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_intergenic"))  
            ccres_count_non_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_non_coding"))
            ccres_count_not_sense_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_not_sense_coding"))
            ccres_count_all = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_all"))

            df_list[[df_name]]["Intergenic Count", dataset_name] = formatC(sum(intergenic_count), format="d", big.mark=",")
            df_list[[df_name]]["cCREs Intergenic Count", dataset_name] = formatC(sum(ccres_count_intergenic), format="d", big.mark=",")
            df_list[[df_name]]["% cCREs Intergenic", dataset_name] = paste0(round(sum(ccres_count_intergenic)/sum(intergenic_count)*100, 2),"%")
            df_list[[df_name]]["cCREs Non-Coding Count", dataset_name] = formatC(sum(ccres_count_non_coding), format="d", big.mark=",")
            df_list[[df_name]]["cCREs All Count", dataset_name] = formatC(sum(ccres_count_all), format="d", big.mark=",")
            df_list[[df_name]]["OCRs Intergenic Count", dataset_name] = formatC(sum(ocr_count_intergenic), format="d", big.mark=",")
            df_list[[df_name]]["% OCRs Intergenic", dataset_name] = paste0(round(sum(ocr_count_intergenic)/sum(intergenic_count)*100, 2),"%")
            df_list[[df_name]]["OCRs Non-Coding Count", dataset_name] = formatC(sum(ocr_count_non_coding), format="d", big.mark=",")
            df_list[[df_name]]["OCR All Count", dataset_name] = formatC(sum(ocr_count_all), format="d", big.mark=",")

            peak_name = "atac_peaks"
            if(dir.exists(file.path(ref_dir, peak_name))) {
                ocr_count_intergenic = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_intergenic"))
                ocr_count_non_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_non_coding"))
                ocr_count_not_sense_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_not_sense_coding"))
                ocr_count_all = Seurat::Read10X(file.path(ref_dir, peak_name, "ocr_count_all"))
                ccres_count_intergenic = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_intergenic"))  
                ccres_count_non_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_non_coding"))
                ccres_count_not_sense_coding = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_not_sense_coding"))
                ccres_count_all = Seurat::Read10X(file.path(ref_dir, peak_name, "cCREs_count_all"))
                df_list[[df_name]][" cCREs Intergenic Count", dataset_name] = formatC(sum(ccres_count_intergenic), format="d", big.mark=",")
                df_list[[df_name]][" % cCREs Intergenic", dataset_name] = paste0(round(sum(ccres_count_intergenic)/sum(intergenic_count)*100, 2),"%")
                df_list[[df_name]][" cCREs Non-Coding Count", dataset_name] = formatC(sum(ccres_count_non_coding), format="d", big.mark=",")
                df_list[[df_name]][" cCREs All Count", dataset_name] = formatC(sum(ccres_count_all), format="d", big.mark=",")
                df_list[[df_name]][" OCRs Intergenic Count", dataset_name] = formatC(sum(ocr_count_intergenic), format="d", big.mark=",")
                df_list[[df_name]][" % OCRs Intergenic", dataset_name] = paste0(round(sum(ocr_count_intergenic)/sum(intergenic_count)*100, 2),"%")
                df_list[[df_name]][" OCRs Non-Coding Count", dataset_name] = formatC(sum(ocr_count_non_coding), format="d", big.mark=",")
                df_list[[df_name]][" OCR All Count", dataset_name] = formatC(sum(ocr_count_all), format="d", big.mark=",")

            }
        }
    }
    df_list
}


