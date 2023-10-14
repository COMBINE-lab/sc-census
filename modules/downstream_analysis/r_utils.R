packages = c(
    "Seurat",
    "SingleCellExperiment",
    "findPC",
    "pheatmap",
    "RColorBrewer",
    "ggpubr",
    "ggplot2",
    "igraph",
    "openxlsx",
    "HGNChelper"
)

filter_empty <- function(af_raw_sce) {
    unfiltered_counts <- af_raw_sce@assays@data$spliced +
        af_raw_sce@assays@data$unspliced +
        af_raw_sce@assays@data$ambiguous
    e.out <- emptyDrops(unfiltered_counts)
    is.cell <- e.out$FDR <= 0.01
    is.cell[is.na(is.cell)] <- FALSE
    is.cell
}

create_seurat_obj <-
    function(m,
             sce,
             min.cells = 0,
             min.features = 0,
             gid2name_df = NULL,
             ...) {
        colnames(m) <- colnames(sce)
        if (is.null(gid2name_df)) {
            rownames(m) <- rownames(sce)
        } else {
            gid2name = gid2name_df$V2
            names(gid2name) = gid2name_df$V1
            rownames(m) <- gid2name[rownames(sce)]
        }
        seurat_obj <-
            CreateSeuratObject(
                counts = m,
                min.cells = min.cells,
                min.features = min.features,
                ...
            )
        seurat_obj
    }

find_seurat_clusters <-
    function(seurat_obj,
             findPC_npcs = 100,
             pcs = NULL,
             clustering_resolution = 0.5,
             redo_knn = FALSE,
             verbose = FALSE) {
        # if scale.data is not there, find variable features
        if (IsMatrixEmpty(GetAssayData(seurat_obj, slot = "scale.data"))) {
            seurat_obj <-
                PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
            
            # From here https://satijalab.org/seurat/articles/sctransform_vignette.html
            seurat_obj <-
                SCTransform(
                    seurat_obj,
                    method = "glmGamPoi",
                    vars.to.regress = "percent.mt",
                    assay = seurat_obj@active.assay,
                    verbose = verbose
                )
        }
        # if pca hasn't been computed, do it
        if (!"pca" %in% names(seurat_obj@reductions)) {
            if (is.null(pcs)) {
                pcs <- max(findPC(
                    sdev =  RunPCA(seurat_obj,
                                   npcs = findPC_npcs,
                                   verbose = verbose)[["pca"]]@stdev,
                    number = findPC_npcs,
                    method = "all"
                ))
            }
            seurat_obj <- RunPCA(seurat_obj,
                                 npcs = max(pcs),
                                 verbose = verbose)
        } else {
            if (is.null(pcs)) {
                pcs <- length(seurat_obj[["pca"]])
            }
        }
        
        if (!"RNA_nn" %in% names(seurat_obj@graphs) || redo_knn) {
            if (length(pcs) == 1)
                pcs <- 1:pcs
            seurat_obj <-
                FindNeighbors(seurat_obj, dims = pcs, verbose = verbose)
        }
        seurat_obj <- FindClusters(seurat_obj,
                                   resolution = clustering_resolution,
                                   verbose = verbose)
        
        seurat_clusters = as.numeric(as.character(seurat_obj$seurat_clusters)) + 1
        Idents(seurat_obj) <- as.factor(seurat_clusters)
        
        list(
            seurat_clusters = seurat_clusters,
            n_clusters = max(seurat_clusters),
            clustering_resolution = clustering_resolution,
            pcs = pcs,
            seurat_obj = seurat_obj,
            knn_graph = seurat_obj@graphs$SCT_snn
        )
    }

seurat_clusters <-
    function(seurat_obj,
             n_clusters = NULL,
             pcs = NULL,
             start_resolution = 0.5,
             step_size = 0.1,
             redo_knn = FALSE,
             verbose = FALSE) {
        clusters_result <- find_seurat_clusters(
            seurat_obj,
            pcs = pcs,
            clustering_resolution = start_resolution,
            redo_knn = redo_knn,
            verbose = verbose
        )
        
        # if we don't need to tune clustering_resolution, just return
        if (is.null(n_clusters)) {
            return(clusters_result)
        }
        
        # if we need to have a fixed # of clusters, tune the clustering_resolution to achieve this
        tuned_resolution <- start_resolution
        last_nclusters <- clusters_result$n_clusters
        
        
        while (last_nclusters != n_clusters) {
            
            if (tuned_resolution <= step_size) {
                step_size = step_size/2
            }
            
            if (last_nclusters > n_clusters) {
                tuned_resolution <- round(tuned_resolution - step_size, 5)
            } else {
                tuned_resolution <- round(tuned_resolution + step_size, 5)
            }
            
            tuned_resolution = abs(tuned_resolution)
            
            if (verbose) {
                message(paste0("last_nclusters = ", last_nclusters, 
                            "\n tuned_resolution = ", tuned_resolution, 
                            "\n step_size = ", step_size))
            }
            
            clusters_result <-
                find_seurat_clusters(
                    seurat_obj = clusters_result$seurat_obj,
                    pcs = pcs,
                    clustering_resolution = tuned_resolution,
                    verbose = verbose
                )
            
            if ((last_nclusters > n_clusters) != (clusters_result$n_clusters > n_clusters)) {
                step_size <- round(step_size / 1.5, 5)
            }
            # avoid infinit loop
            last_nclusters <-
                length(unique(clusters_result$seurat_clusters))
        }
        clusters_result
    }

get_variable_features <- function(normalized_seurat_obj,
                                  nfeatures = 3000,
                                  verbose = FALSE) {
    seurat_obj <- FindVariableFeatures(normalized_seurat_obj,
                                       nfeatures = nfeatures,
                                       verbose = verbose)
    VariableFeatures(seurat_obj)
}

sctype <- function(clusters_result,
                   tissue = "Immune system") {
    require(dplyr)
    require(HGNChelper)
    require(Seurat)
    # get data
    seurat_obj <- clusters_result$seurat_obj
    seurat_clusters <- clusters_result$seurat_clusters
    
    # load gene set preparation function
    source(
        "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
    )
    # load cell type annotation function
    source(
        "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
    )
    db_ <-
        "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    
    
    # prepare gene sets
    suppressWarnings({
        suppressMessages({
            gs_list <- gene_sets_prepare(db_, tissue)
        })
    })
    
    es.max <-
        sctype_score(
            scRNAseqData = seurat_obj[[seurat_obj@active.assay]]@scale.data,
            scaled = TRUE,
            gs = gs_list$gs_positive,
            gs2 = gs_list$gs_negative
        )
    
    # merge by cluster
    cL_resutls <-
        do.call("rbind", lapply(unique(seurat_clusters), function(cl) {
            es.max.cl <-
                sort(rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_clusters == cl,])]), decreasing = TRUE)
            head(data.frame(
                cluster = cl,
                type = names(es.max.cl),
                scores = es.max.cl,
                ncells = sum(seurat_clusters == cl)
            ),
            10)
        }))
    sctype_scores <-
        cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells /
                           4] <- "Unknown"
    sctype_scores <-
        sctype_scores[order(as.numeric(sctype_scores$cluster)), ]
    seurat_obj@meta.data$sctype_clusters <- ""
    for (j in unique(sctype_scores$cluster)) {
        cl_type <- sctype_scores[sctype_scores$cluster == j, ]
        seurat_obj@meta.data$sctype_clusters[seurat_clusters == j] <-
            as.character(cl_type$type[1])
    }
    
    
    clusters_result$sctype_clusters <-
        seurat_obj@meta.data$sctype_clusters
    clusters_result$sctype_score <- sctype_scores
    clusters_result$seurat_obj <- seurat_obj
    
    # return updated clusters_result
    return(clusters_result)
}

get_clusters_from_count_type_multiple <-
    function(input_dataset_list,
             pcs = NULL,
             clustering_resolution = 0.5,
             fontsize = NULL,
             metrics = c("ari", "knn_spearman"),
             n_threads = 2,
             tissue = "Immune system",
             ...) {
        if (n_threads < 2) {
            registerDoSEQ()
        } else {
            registerDoParallel(cores = n_threads)
        }
        
        clusters_result_list <-
            foreach (input_dataset = input_dataset_list, .packages = packages) %dopar% {
                # apply standard pipeline
                seurat_obj_list <- input_dataset$seurat_obj_list
                dataset_id <- input_dataset$dataset_id
                n_clusters <- input_dataset$n_clusters
                
                # create clusters result list
                clusters_result_list <-
                    lapply(seurat_obj_list, function(seurat_obj_set) {
                        seurat_obj <- seurat_obj_set$seurat_obj
                        count_type <- seurat_obj_set$count_type
                        
                        # apply pipeline
                        clusters_result <- seurat_clusters(
                            seurat_obj = seurat_obj,
                            n_clusters = n_clusters,
                            start_resolution = clustering_resolution,
                            pcs = pcs,
                            verbose = verbose
                        )
                        
                        ####################################################################################################################################
                        # scType procedure
                        clusters_result <- sctype(clusters_result)
                        #############################
                        
                        list("name" = count_type,
                             "clusters_result" = clusters_result)
                    })
                
                heatmap_list <- make_heatmaps(
                    clusters_result_list = clusters_result_list,
                    dataset_id = dataset_id,
                    fontsize = fontsize,
                    metrics = metrics,
                    ...
                )
                
                list(
                    dataset_id = dataset_id,
                    clusters_result_list = clusters_result_list,
                    heatmap_list = heatmap_list
                )
            }
        
        names(clusters_result_list) <-
            sapply(clusters_result_list, function(x)
                x$dataset_id)
        clusters_result_list
        
    }



create_count_type_dataset_list <-
    function(pq_list,
             nfeatures = 3000,
             gid2name_df = NULL,
             fix_n_clusters = TRUE,
             clustering_resolution = 0.5) {
        input_dataset_list <- lapply(pq_list, function(pq) {
            dataset_id <- as.character(pq@dataset_id)
            sce <- pq@sce
            is.cell <- filter_empty(sce)
            sce <- sce[, is.cell]
            seurat_obj_list <- list(
                STD =  list(
                    count_type = "STD",
                    seurat_obj = create_seurat_obj(
                        m = sce@assays@data$spliced + sce@assays@data$ambiguous,
                        sce = sce,
                        gid2name_df = gid2name_df,
                        project = "STD"
                    )
                ),
                USA = list(
                    count_type = "USA",
                    seurat_obj = create_seurat_obj(
                        m = sce@assays@data$spliced + sce@assays@data$unspliced + sce@assays@data$ambiguous,
                        sce = sce,
                        gid2name_df = gid2name_df,
                        project = "USA"
                    )
                ),
                U = list(
                    count_type = "U",
                    seurat_obj = create_seurat_obj(
                        m = sce@assays@data$unspliced,
                        sce = sce,
                        gid2name_df = gid2name_df,
                        project = "U"
                    )
                )
            )
            
            # n_clusters_list
            if (fix_n_clusters) {
                seurat_obj <- seurat_obj_list$STD$seurat_obj
                clusters_result = seurat_clusters(seurat_obj = seurat_obj,
                                                  start_resolution = clustering_resolution)
                n_clusters = clusters_result$n_clusters
                
            } else {
                n_clusters = NULL
            }
            
            list(
                dataset_id = dataset_id,
                seurat_obj_list = seurat_obj_list,
                n_clusters = n_clusters
            )
            
        })
        
        names(input_dataset_list) <-
            sapply(input_dataset_list, function(x)
                x$dataset_id)
        input_dataset_list
    }


get_variable_features <- function(normalized_seurat_obj,
                                  nfeatures = 3000,
                                  verbose = FALSE) {
    seurat_obj <- FindVariableFeatures(normalized_seurat_obj,
                                       nfeatures = nfeatures,
                                       verbose = verbose)
    VariableFeatures(seurat_obj)
}


make_heatmaps_sub <- function(h_mtx,
                              main,
                              fontsize = NULL,
                              breaks_max = 1,
                              display_numbers = TRUE,
                              angle_col = 45,
                              ...) {
    my.breaks <- c(seq(0, breaks_max, by = breaks_max / 50))
    my.colors <-
        colorRampPalette(colors = c("white", "yellow", "orange", "red"))(length(my.breaks))
    
    if (is.null(fontsize))
        fontsize = round(50 / ncol(h_mtx))
    # if (is.null(fontsize)) fontsize = round(log(ncol(h_mtx),15),1)
    
    pheatmap(
        h_mtx,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        angle_col = angle_col,
        drop_levels = FALSE,
        color = my.colors,
        breaks = my.breaks,
        main = main,
        silent = TRUE,
        display_numbers = display_numbers,
        legend = FALSE,
        fontsize = fontsize,
        ...
    )[[4]]
    
}

compute_ari <- function(clusters_result_list) {
    # init matrix
    h_mtx <- matrix(
        0,
        nrow = length(clusters_result_list),
        ncol = length(clusters_result_list)
    )
    colnames(h_mtx) <-
        rownames(h_mtx) <- sapply(clusters_result_list, function(x)
            x$name)
    
    # compute similarity
    for (c_std in 1:length(clusters_result_list)) {
        for (c_query in c_std:length(clusters_result_list)) {
            # ARI
            ari <- round(
                mclust::adjustedRandIndex(
                    clusters_result_list[[c_std]]$clusters_result$seurat_clusters,
                    clusters_result_list[[c_query]]$clusters_result$seurat_clusters
                ),
                2
            )
            h_mtx[c_std, c_query] <- ari
            h_mtx[c_query, c_std] <- ari
        }
    }
    h_mtx
}

compute_metric <- function(clusters_result_list, fn) {
    # init matrix
    h_mtx <- matrix(
        0,
        nrow = length(clusters_result_list),
        ncol = length(clusters_result_list)
    )
    colnames(h_mtx) <-
        rownames(h_mtx) <- sapply(clusters_result_list, function(x)
            x$name)
    
    # compute similarity
    for (c_std in 1:length(clusters_result_list)) {
        for (c_query in c_std:length(clusters_result_list)) {
            # ARI
            ari <- round(
                fn(
                    as.numeric(clusters_result_list[[c_std]]$clusters_result$seurat_clusters),
                    as.numeric(clusters_result_list[[c_query]]$clusters_result$seurat_clusters)
                ),
                2
            )
            h_mtx[c_std, c_query] <- ari
            h_mtx[c_query, c_std] <- ari
        }
    }
    h_mtx
}

compute_knn_spearman_r <- function(clusters_result_list) {
    # init matrix
    h_mtx <- matrix(
        0,
        nrow = length(clusters_result_list),
        ncol = length(clusters_result_list)
    )
    colnames(h_mtx) <-
        rownames(h_mtx) <- sapply(clusters_result_list, function(x)
            x$name)
    
    
    # compute similarity
    for (c_std in 1:length(clusters_result_list)) {
        for (c_query in c_std:length(clusters_result_list)) {
            # kNN separman R
            spearman_cor = knn_spearman_cor(
                clusters_result_list[[c_std]]$clusters_result$knn_graph,
                clusters_result_list[[c_query]]$clusters_result$knn_graph
            )$spearman_cor
            
            h_mtx[c_std, c_query] <- spearman_cor
            h_mtx[c_query, c_std] <- spearman_cor
        }
    }
    h_mtx
}


make_heatmaps <- function(clusters_result_list,
                          dataset_id ,
                          main = NA,
                          metrics = c("ari", "knn_spearman", "nmi", "fmi"),
                          fontsize = NULL,
                          ...) {
    if (is.na(main)) {
        main <- paste0("dataset #", dataset_id)
    }
    
    result_list = list()
    # ari
    if ("ari" %in% metrics) {
        # ari_h_mtx <- compute_ari(clusters_result_list)
        ari_h_mtx <- compute_metric(clusters_result_list, aricode::ARI)
        result_list$ari_heatmap <- make_heatmaps_sub(h_mtx = ari_h_mtx,
                                   main = paste("ARI", main),
                                   fontsize = fontsize,
                                   ...)
    }
    
    # knn spearman
    if ("knn_spearman" %in% metrics) {
        knn_h_mtx <- compute_knn_spearman_r(clusters_result_list)
        result_list$knn_spearman_heatmap <- make_heatmaps_sub(h_mtx = knn_h_mtx,
                                        main = paste("KNN spearman", main),
                                        fontsize = fontsize,
                                        ...)
    } 

    # NMI
    if ("nmi" %in% metrics) {
        nmi_h_mtx <- compute_metric(clusters_result_list, aricode::NMI)
        result_list$nmi_heatmap <- make_heatmaps_sub(h_mtx = nmi_h_mtx,
                                        main = paste("NMI", main),
                                        fontsize = fontsize,
                                        ...)
    }

    # FMI
    if ("fmi" %in% metrics) {
        # clevr::fowlkes_mallows
        fmi_h_mtx <- compute_metric(clusters_result_list, dendextend::FM_index)
        result_list$fmi_heatmap <- make_heatmaps_sub(h_mtx = fmi_h_mtx,
                                        main = paste("FMI", main),
                                        fontsize = fontsize,
                                        ...)
    }
    result_list
}

# functions are adapted from https://github.com/davisidarta/pca-eval/blob/master/pcaeval/pca_eval.py
geodesic_distance <- function(knn_adj_mtx, weighted = TRUE) {
    require(igraph)
    ds_g <-
        distances(graph_from_adjacency_matrix(as.matrix(knn_adj_mtx), weighted = weighted))
}

knn_spearman_cor <- function(knn_1, knn_2, weighted = TRUE) {
    # calculate distance matrix
    dist_knn_1 <- geodesic_distance(knn_1, weighted = weighted)
    dist_knn_2 <- geodesic_distance(knn_2, weighted = weighted)
    
    # spearman cor of the flatten matrix
    # use only the upper triangle, excluding diagonal.
    spearman_cor <- cor(c(dist_knn_1[upper.tri(dist_knn_1)]),
                        c(dist_knn_2[upper.tri(dist_knn_2)]),
                        method = "spearman")
    
    list(
        dist_knn_1 = dist_knn_1,
        dist_knn_2 = dist_knn_2,
        spearman_cor = spearman_cor
    )
}

match_clusters <-
    function(seurat_obj_set_ref,
             # this is the standard result
             seurat_obj_set_query,
             dataset_id) {
        # get number of clusters and the cluster assignment in each clustering result
        name_ref <- seurat_obj_set_ref$name
        name_query <- seurat_obj_set_query$name
        nc_ref <- seurat_obj_set_ref$clusters_result$n_clusters
        c_ref <- seurat_obj_set_ref$clusters_result$seurat_clusters
        nc_query <- seurat_obj_set_query$clusters_result$n_clusters
        c_query <-
            seurat_obj_set_query$clusters_result$seurat_clusters
        
        # init matrix
        # this matrix is used to see if ref clusters are divided into sub-clusters in query
        # each column sums to one
        cell_assignment_mtx <-
            matrix(
                data = 0,
                nrow = nc_query,
                ncol = nc_ref,
                dimnames = list(paste0(name_query, "_", seq(nc_query)), paste0(name_ref, "_", seq(nc_ref)))
            )
        
        for (c_ref_id in seq(nc_ref)) {
            # grab cell in the current processing c_ref cluster
            c_ref_cells <- which(c_ref == c_ref_id)
            
            # see what c_query clusters do those cells fall into
            ref_in_query_table <- c(table(c_query[c_ref_cells]))
            
            # record the distribution
            cell_assignment_mtx[paste0(name_query, "_",  names(ref_in_query_table)), c_ref_id] <-
                ref_in_query_table / length(c_ref)
        }
        # ref_in_query_mtx = ref_in_query_mtx/max(ref_in_query_mtx)
        
        main <-
            paste(
                "Assignment difference in",
                name_query,
                "and",
                name_ref,
                "clusters",
                "in dataset",
                dataset_id
            )
        hm <-
            make_heatmaps_sub(
                h_mtx = cell_assignment_mtx,
                main = main,
                breaks_max = max(cell_assignment_mtx)
            )
        
        name <- paste0(name_ref, "-", name_query)
        
        matched_clusters_list = list(
            name = name,
            cell_assignment_mtx = cell_assignment_mtx,
            heatmap = hm,
            description = paste(
                "The whole matrix sums to 1.",
                "The color and value of each bin shows the portion of cells in the bin."
            )
        )
        
        list(
            name = name,
            dataset_id = dataset_id,
            matched_clusters_list = matched_clusters_list
        )
    }

match_clusters_for_all <-
    function(input_dataset_list,
             n_threads = 2,
             verbose = FALSE,
             ...) {
        if (n_threads < 2) {
            registerDoSEQ()
        } else {
            registerDoParallel(cores = n_threads)
        }
        
        # find the markers of each seurat object
        clusters_result_list <-
            foreach(input_dataset = input_dataset_list, .packages = packages) %dopar% {
                dataset_id = input_dataset$dataset_id
                
                matched_clusters_list <- list()
                for (assay1 in seq(length(input_dataset$clusters_result_list) -
                                   1)) {
                    seurat_obj_set_ref = input_dataset$clusters_result_list[[assay1]]
                    for (assay2 in seq(assay1 + 1,
                                       length(input_dataset$clusters_result_list))) {
                        seurat_obj_set_query = input_dataset$clusters_result_list[[assay2]]
                        matched_clusters <-
                            match_clusters(seurat_obj_set_ref,
                                           seurat_obj_set_query,
                                           dataset_id)
                        matched_clusters_list[[matched_clusters$name]] = matched_clusters
                    }
                }
                
                names(matched_clusters_list) <-
                    sapply(matched_clusters_list, function(x)
                        x$name)
                
                # input_dataset$cluster_id_df <- do.call(cbind.data.frame, lapply(input_dataset$clusters_result_list, FUN = function(seurat_obj_set){
                #     seurat_obj_set$jaccard_result$matched_clusters_df
                # }))
                input_dataset$heatmap_list$matched_clusters_list = matched_clusters_list
                input_dataset
            }
        names(clusters_result_list) = sapply(clusters_result_list, function(x)
            x$dataset_id)
        
        clusters_result_list
    }

add_dimnames <- function(m, ref) {
    rownames(m) = rownames(ref)
    colnames(m) = colnames(ref)
    m
}

find_clusters_for_all <-
    function(pq_list,
             nfeatures = 3000,
             clustering_resolution = 0.5,
             n_threads = 2,
             pcs = NULL,
             verbose = FALSE,
             ...) {
        if (n_threads < 2) {
            registerDoSEQ()
        } else {
            registerDoParallel(cores = n_threads)
        }
        
        # find the markers of each seurat object
        clusters_result_list <-
            foreach(pq = pq_list, .packages = packages) %dopar% {
                # get data
                dataset_id <- as.character(pq@dataset_id)
                sce <- pq@sce
                
                # preprocessing
                is.cell <- filter_empty(sce)
                sce <- sce[, is.cell]
                
                # create seurat object
                seurat_obj = create_seurat_obj(
                    m = sce@assays@data$spliced + sce@assays@data$ambiguous,
                    sce = sce,
                    gid2name_df = gid2name_df,
                    assay = "STD"
                )
                
                # add assay for USA
                seurat_obj[["USA"]] = CreateAssayObject(
                    add_dimnames(
                        sce@assays@data$spliced + sce@assays@data$unspliced + sce@assays@data$ambiguous,
                        seurat_obj
                    )
                )
                
                # add assay for U
                seurat_obj[["U"]] = CreateAssayObject(add_dimnames(sce@assays@data$unspliced, seurat_obj))
                
                # do clustering using STD count
                ## apply the standard pipeline
                clusters_result <- seurat_clusters(
                    seurat_obj = seurat_obj,
                    nfeatures = nfeatures,
                    start_resolution = clustering_resolution,
                    pcs = pcs,
                    verbose = verbose
                )
                
                ## sctype
                clusters_result <- sctype(clusters_result)
                
                ## seurat transfer learning
                clusters_result$seurat_obj <-
                    RunAzimuth(clusters_result$seurat_obj,
                               assay = "SCT",
                               reference = "pbmcref")
                
                clusters_result
            }
        
        names(clusters_result_list) = sapply(clusters_result_list, function(x)
            x$dataset_id)
        
        clusters_result_list
    }

create_multimodal_dataset_for_all <- function(pq_list,
                                              clustering_resolution = 0.5,
                                              n_threads = 2,
                                              pcs = NULL,
                                              verbose = FALSE,
                                              ...) {
    if (n_threads < 2) {
        registerDoSEQ()
    } else {
        registerDoParallel(cores = n_threads)
    }
    
    multimodal_dataset_list <-
        foreach(pq = pq_list, .packages = packages) %dopar% {
            # get data
            dataset_id <- as.character(pq@dataset_id)
            sce <- pq@sce
            
            # preprocessing
            is.cell <- filter_empty(sce)
            sce <- sce[, is.cell]
            
            # create seurat object
            seurat_obj = create_seurat_obj(
                m = sce@assays@data$spliced + sce@assays@data$ambiguous,
                sce = sce,
                gid2name_df = gid2name_df,
                assay = "STD"
            )
            
            # add assay for USA
            seurat_obj[["USA"]] = CreateAssayObject(
                add_dimnames(
                    sce@assays@data$spliced + sce@assays@data$unspliced + sce@assays@data$ambiguous,
                    seurat_obj
                )
            )
            
            # add assay for U
            seurat_obj[["U"]] = CreateAssayObject(add_dimnames(sce@assays@data$unspliced, seurat_obj))
            
            list(dataset_id = dataset_id,
                 seurat_obj = seurat_obj)
        }
    
    # output
    names(multimodal_dataset_list) = sapply(multimodal_dataset_list, function(x)
        x$dataset_id)
    multimodal_dataset_list
}

find_cell_types_for_all <-
    function(multimodal_dataset_list_all,
             clustering_resolution = 0.5,
             n_threads = 2,
             pcs = NULL,
             verbose = FALSE,
             ...) {
        if (n_threads < 2) {
            registerDoSEQ()
        } else {
            registerDoParallel(cores = n_threads)
        }
        
        # find the markers of each seurat object
        clusters_result_list_all <-
            foreach(multimodal_dataset_list = multimodal_dataset_list_all,
                    .packages = packages) %dopar% {
                        # apply the standard pipeline
                        clusters_result <- seurat_clusters(
                            seurat_obj = multimodal_dataset_list$seurat_obj,
                            start_resolution = clustering_resolution,
                            verbose = verbose
                        )
                        
                        # sctype
                        clusters_result <- sctype(clusters_result)
                        
                        # seurat transfer learning
                        clusters_result$seurat_obj <-
                            RunAzimuth(
                                clusters_result$seurat_obj,
                                assay = clusters_result$seurat_obj@active.assay,
                                reference = "pbmcref",
                                verbose = verbose
                            )
                        
                        # write results
                        clusters_result$predicted.celltype.l1 <-
                            clusters_result$seurat_obj$predicted.celltype.l1
                        clusters_result$predicted.celltype.l2 <-
                            clusters_result$seurat_obj$predicted.celltype.l2
                        
                        list(dataset_id = multimodal_dataset_list$dataset_id,
                             clusters_result = clusters_result)
                    }
        
        names(clusters_result_list_all) = sapply(clusters_result_list_all, function(x)
            x$dataset_id)
        
        clusters_result_list_all
    }

# TODO: try presto to compute auc as well

find_markers_for_all <-
    function(multimodal_dataset_list_all,
             n_threads = 2,
             ident_names = c(
                 "seurat_clusters",
                 "sctype_clusters",
                 "predicted.celltype.l1",
                 "predicted.celltype.l2"
             ),
             assay_names = c("STD", "USA", "U"),
             verbose = FALSE,
             ...) {
        if (n_threads < 2) {
            registerDoSEQ()
        } else {
            registerDoParallel(cores = n_threads)
        }
        
        multimodal_dataset_list_all <-
            foreach(multimodal_dataset_list = multimodal_dataset_list_all,
                    .packages = packages) %dopar% {
                        # get data
                        dataset_id <-
                            as.character(multimodal_dataset_list$dataset_id)
                        
                        if (verbose)
                            message(paste0("Processing dataset - ", dataset_id))
                        
                        
                        # we find markers using each of the count types for each of the identification sets
                        marker_df_list <-
                            apply(
                                expand.grid(assay_names, ident_names, stringsAsFactors = FALSE),
                                MARGIN = 1,
                                FUN = function(assay_ident)
                                {
                                    # get correct assay
                                    name = paste0(assay_ident[1], "-", assay_ident[2])
                                    
                                    if (verbose)
                                        message(paste0("  - Finding markers for ", name))
                                    
                                    seurat_obj <-
                                        multimodal_dataset_list$clusters_result$seurat_obj
                                    Idents(seurat_obj) <-
                                        assay_ident[2]
                                    
                                    # we find markers for each of the identification set
                                    all_markers_df <-
                                        FindAllMarkers(
                                            object = seurat_obj,
                                            assay = assay_ident[1],
                                            verbose = FALSE,
                                            ...
                                        )
                                    # FindAllMarkers(object = seurat_obj, verbose = verboase, ...)
                                    
                                    list(
                                        name = name,
                                        assay_name = assay_ident[1],
                                        ident_name = assay_ident[2],
                                        all_markers_df = all_markers_df
                                    )
                                }
                            )
                        names(marker_df_list) = sapply(marker_df_list, function(x)
                            x$name)
                        
                        multimodal_dataset_list$marker_df_list = marker_df_list
                        multimodal_dataset_list
                    }
        names(multimodal_dataset_list_all) = sapply(multimodal_dataset_list_all, function(x)
            x$dataset_id)
        
        multimodal_dataset_list_all
    }

adjust_report_size <- function(assay_name,
                               assays.use.report_size_factor,
                               report_size,
                               report_size_factor,
                               verbose = FALSE) {
    # if report size is not set, then return everything
    if (is.null(report_size)) {
        if (!is.null(assays.use.report_size_factor)) {
            if (assay_name %in% assays.use.report_size_factor)
                .say(verbose,
                     "report_size is NULL, unable to resize report_size")
        }
        return(report_size)
    }
    
    # now report size is a number
    # if need resize
    if (!is.null(assays.use.report_size_factor)) {
        if (assay_name %in% assays.use.report_size_factor) {
            report_size = report_size_factor * report_size
        }
    }
    
    report_size
}

filter_deg_df = function(deg_df,
                         assay_name,
                         only_positive = FALSE,
                         report_size = NULL,
                         report_size_factor = 2,
                         max.p_val_adj = NULL,
                         min.avg_log2FC = NULL,
                         assays.use.report_size_factor = NULL) {
    # filter positive degs if needed
    if (only_positive)
        deg_df = deg_df[deg_df$avg_log2FC > 0,]
    
    # filter large p value if needed
    if (!is.null(max.p_val_adj))
        deg_df = deg_df[deg_df$p_val_adj < max.p_val_adj,]
    
    # filter low FC if needed
    if (!is.null(min.avg_log2FC))
        deg_df = deg_df[abs(deg_df$avg_log2FC) > min.avg_log2FC,]
    
    # rank degs
    deg_df = deg_df[with(deg_df,
                         order(-abs(avg_log2FC), -pct.1)), ]
    
    
    # could be NULL
    adjusted_report_size = adjust_report_size(
        assay_name,
        assays.use.report_size_factor,
        report_size,
        report_size_factor,
        verbose = verbose
    )
    
    if (is.null(adjusted_report_size)) {
        adjusted_report_size = nrow(deg_df)
    }
    
    
    deg_df = deg_df[0:min(adjusted_report_size,
                          nrow(deg_df)),]
    
    deg_df
}

get_degs_for_all <-
    function(multimodal_dataset_list_all,
             parent_dir = "deg_of_datasets",
             only_positive = FALSE,
             report_size = NULL,
             report_size_factor = 2,
             assays.use.report_size_factor = NULL,
             max.p_val_adj = NULL,
             min.avg_log2FC = NULL,
             # return.deg.list = FALSE,
             write_to_files = FALSE,
             verbose = FALSE) {
        dir_create(write_to_files, parent_dir)
        
        deg_list_all = list()
        
        for (multimodal_dataset_list in multimodal_dataset_list_all) {
            # get dataset_id
            dataset_id = multimodal_dataset_list$dataset_id
            
            if (verbose)
                cat(paste0("Processing dataset - ", dataset_id, "\n"))
            
            dataset_dir <-
                dir_create(write_to_files,
                           parent_dir,
                           paste0("dataset_", dataset_id, "_degs"))
            
            dataset_deg_list <- list()
            for (marker_df_list in multimodal_dataset_list$marker_df_list) {
                # get filename
                name = marker_df_list$name
                
                assay_name = strsplit(name, split = "-")[[1]][1]
                ident_name = strsplit(name, split = "-")[[1]][2]
                # assay_name =  marker_df_list$assay_name
                # ident_name =  marker_df_list$ident_name
                
                all_markers_df <-
                    as.data.frame(matrix(
                        NA,
                        nrow = 0,
                        ncol = ncol(marker_df_list$all_markers_df)
                    ))
                
                for (cluster_name in as.character(unique(marker_df_list$all_markers_df$cluster))) {
                    log_msg =
                        paste0(
                            "Identity: ",
                            ident_name,
                            ", assay: ",
                            assay_name,
                            ", cluster: ",
                            cluster_name
                        )
                    if (verbose)
                        cat(paste0("  - Processing ", log_msg, "\n"))
                    
                    # find the degs for the current cluster
                    cluster_degs = marker_df_list$all_markers_df[marker_df_list$all_markers_df$cluster == cluster_name, ]
                    
                    cluster_degs = filter_deg_df(
                        cluster_degs,
                        assay_name = assay_name,
                        only_positive = only_positive,
                        report_size = report_size,
                        report_size_factor = report_size_factor,
                        max.p_val_adj = max.p_val_adj,
                        min.avg_log2FC = min.avg_log2FC,
                        assays.use.report_size_factor = assays.use.report_size_factor
                    )
                    
                    if (isEmpty(cluster_degs) & verbose) {
                        cat("    No DEG passed filtering in ",
                            log_msg,
                            "; Skip",
                            "\n")
                        next
                    }
                    
                    # record the filtered degs
                    all_markers_df <-
                        rbind.data.frame(all_markers_df, cluster_degs)
                    
                    dataset_deg_list[[ident_name]][[cluster_name]][[assay_name]] <-
                        list(
                            assay_name = assay_name,
                            ident_name = ident_name,
                            cluster_name = cluster_name,
                            deg_df = cluster_degs
                        )
                }
                
                # wrtite the filtered degs
                fn = file.path(dataset_dir, paste0(name, "_deg", ".csv"))
                
                if (!isEmpty(all_markers_df) & write_to_files) {
                    write.csv(all_markers_df,
                              file = fn,
                              quote = FALSE)
                }
            }
            
            deg_list_all[[dataset_id]] = list(dataset_id = dataset_id,
                                              dataset_deg_list = dataset_deg_list)
            
        }
        return(deg_list_all)
    }

dir_create <- function(write_to_files = FALSE, ...) {
    d <- file.path(...)
    # create dir if needed
    if (!dir.exists(d) & write_to_files)
        dir.create(d,
                   recursive = TRUE,
                   showWarnings = FALSE)
    
    d
    
}

write_overlap_degs_across_assays_venn_diagrams_for_all <-
    function(deg_list_all,
             parent_dir = "deg_of_datasets",
             only_positive = FALSE,
             report_size = NULL,
             report_size_factor = 2,
             assays.use.report_size_factor = NULL,
             max.p_val_adj = NULL,
             min.avg_log2FC = NULL,
             show_empty_plot = TRUE,
             fontsize  = NULL,
             verbose = FALSE) {
        # This function needs ggVenDiagram
        require(ggVennDiagram)
        
        dir_create(TRUE, parent_dir)
        # each datasets gets a file
        for (dataset_deg_list in deg_list_all) {
            # get dataset_id
            dataset_id = dataset_deg_list$dataset_id
            
            # be patient, I'm running
            if (verbose)
                message(paste0("Processing dataset - ", dataset_id))
            
            # create dir if needed
            dataset_dir <-
                dir_create(TRUE,
                           parent_dir,
                           paste0("dataset_", dataset_id, "_degs"))
            # get ven diagram folder name
            dataset_vd_dir <-
                dir_create(TRUE,
                           dataset_dir,
                           paste0("dataset_degs_venn_diagrams"))
            
            # create a venn diagram file for each ident name
            lapply(
                dataset_deg_list$dataset_deg_list,
                FUN = function(ident_deg_list) {
                    # Take the needed info from the first record
                    ident_name = ident_deg_list[[1]][[1]]$ident_name
                    assay_name = ident_deg_list[[1]][[1]]$assay_name
                    
                    # set font size to make the plots pretty
                    if (is.null(fontsize))
                        fontsize = round(30 / length(ident_deg_list))
                    # fontsize = round(log(length(ident_deg_list),15),1)
                    
                    # create a venn diagram for each cluster
                    # R trick, use print if want to write PDF in a function
                    
                    plotlist = list()
                    
                    for (cluster_deg_list in ident_deg_list) {
                        cluster_name = cluster_deg_list[[1]]$cluster_name
                        # get the dataframe
                        
                        deg_list = list()
                        for (assay_deg_list in cluster_deg_list) {
                            assay_name = assay_deg_list$assay_name
                            cluster_name = assay_deg_list$cluster_name
                            # find the degs for the current cluster
                            deg_df =  filter_deg_df(
                                assay_deg_list$deg_df,
                                assay_name = assay_name,
                                only_positive = only_positive,
                                report_size = report_size,
                                report_size_factor = report_size_factor,
                                max.p_val_adj = max.p_val_adj,
                                min.avg_log2FC = min.avg_log2FC,
                                assays.use.report_size_factor = assays.use.report_size_factor
                            )
                            if (nrow(deg_df) > 0) {
                                deg_list[[assay_name]] = deg_df$gene
                            }
                        }
                        
                        # generate plot only if there are more than one deg vector
                        if (length(deg_list) > 1) {
                            # names(deg_list) = lapply(cluster_deg_list, function(assay_deg_list)
                            #     assay_deg_list$assay_name)
                            
                            plotlist[[cluster_name]] = ggVennDiagram(
                                deg_list,
                                color = "black",
                                lwd = 0.8,
                                lty = 1,
                                set_size = fontsize,
                                label_size = fontsize
                            ) +
                                scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
                                labs(title = paste0("cluster - ", cluster_name)) +
                                theme(legend.position = "none")
                        } else {
                            .say(verbose,
                                 "Failed finding DEGs in enough assays; Skip")
                            if (show_empty_plot) {
                                plotlist[[cluster_name]] = ggplot(as.data.frame(matrix())) + geom_blank() +
                                    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
                                    labs(title = paste0(
                                        "cluster - ",
                                        cluster_name,
                                        "\nNo DEG Found"
                                    )) +
                                    theme_minimal()
                            }
                        }
                    }
                    # start a PDF file
                    pdf(
                        file = file.path(dataset_vd_dir, paste0(ident_name, ".pdf")),
                        # The directory you want to save the file in
                        width = 20,
                        # The width of the plot in inches
                        height = 10
                    )
                    
                    if (length(plotlist) > 0) {
                        print(ggarrange(plotlist = plotlist))
                    }
                    
                    dev.off()
                }
            )
        }
    }


write_overlap_degs_across_datasets_venn_diagrams_for_all <-
    function(deg_list_all,
             ident_names = c("predicted.celltype.l1",
                             "predicted.celltype.l2",
                             "sctype_clusters"),
             assay_names = c("USA", "U", "STD"),
             parent_dir = "DEG_overlap_venn_diagrams",
             fontsize = NULL,
             only_positive = FALSE,
             report_size = NULL,
             report_size_factor = 2,
             assays.use.report_size_factor = NULL,
             max.p_val_adj = NULL,
             min.avg_log2FC = NULL,
             show_empty_plot = TRUE,
             verbose = FALSE) {
        require(ggVennDiagram)
        
        dir_create(TRUE, parent_dir)
        # each ident - assay  pair has a file
        result = lapply(ident_names, function(ident_name) {
            ident_file_dir = file.path(parent_dir, ident_name)
            dir_create(TRUE, ident_file_dir)
            
            lapply(assay_names, function(assay_name) {
                # get candidate cell types
                candidate_cell_types = unique(unlist(lapply(deg_list_all, function(deg_list) {
                    dataset_id = deg_list$dataset_id
                    names(deg_list$dataset_deg_list[[ident_name]])
                })))
                
                # prepare plot list
                plotlist = list()
                
                # each cell type gets a plot
                for (curr_cell_type in candidate_cell_types) {
                    curr_deg_list = list()
                    
                    # gather deg list of all datasets
                    for (deg_list in deg_list_all) {
                        if (curr_cell_type %in% names(deg_list$dataset_deg_list[[ident_name]])) {
                            # find the degs for the current cluster
                            deg_df = deg_list$dataset_deg_list[[ident_name]][[curr_cell_type]][[assay_name]]$deg_df
                            
                            # filter if needed
                            deg_df = filter_deg_df(
                                deg_df,
                                assay_name = assay_name,
                                only_positive = only_positive,
                                report_size = report_size,
                                report_size_factor = report_size_factor,
                                max.p_val_adj = max.p_val_adj,
                                min.avg_log2FC = min.avg_log2FC,
                                assays.use.report_size_factor = assays.use.report_size_factor
                            )
                            
                            if (nrow(deg_df) > 0) {
                                curr_deg_list[[deg_list$dataset_id]] = deg_df$gene
                            }
                        }
                    }
                    
                    # set font size
                    if (is.null(fontsize))
                        fontsize = round(10 / length(curr_deg_list))
                    
                    if (length(curr_deg_list) > 1) {
                        # plot and record
                        plotlist[[curr_cell_type]] = ggVennDiagram(
                            curr_deg_list,
                            color = "black",
                            lwd = 0.8,
                            lty = 1
                            ,
                            set_size = fontsize,
                            label_size = fontsize
                        ) +
                            scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
                            labs(title = curr_cell_type) +
                            theme(legend.position = "none")
                    } else {
                        .say(verbose,
                             "Failed finding DEGs in enough assays, skip")
                        if (show_empty_plot) {
                            plotlist[[curr_cell_type]] = ggplot(as.data.frame(matrix())) + geom_blank() +
                                scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
                                labs(title = paste0(curr_cell_type, "\nNo DEG Found")) +
                                theme_minimal()
                            
                        }
                    }
                }
                
                pdf(
                    file = file.path(
                        ident_file_dir,
                        paste0(assay_name, "_", ident_name, ".pdf")
                    ),
                    # The directory you want to save the file in
                    width = 20,
                    # The width of the plot in inches
                    height = 10
                )
                
                if (length(plotlist) > 0) {
                    print(ggarrange(plotlist = plotlist))
                }
                
                dev.off()
            })
        })
    }

.say <- function(verbose = FALSE, ...)  {
    if (verbose) {
        cat(..., "\n")
    }
}


find_enriched_gset_for_all <- function(deg_list_all,
                                       fgsea_sets,
                                       ident_names = c("predicted.celltype.l1",
                                                       "predicted.celltype.l2",
                                                       "sctype_clusters"),
                                       assay_names = c("USA", "U", "STD"),
                                       # redo DEG and filter
                                       use_presto = FALSE,
                                       min.presto_auc = NULL,
                                       only_positive = FALSE,
                                       report_size = 2000,
                                       report_size_factor = 2,
                                       assays.use.report_size_factor = NULL,
                                       max.p_val_adj = NULL,
                                       min.avg_log2FC = NULL,
                                       max.fgsea_padj = 0.05,
                                       show_empty_plot = TRUE,
                                       write_to_files = TRUE,
                                       multimodal_dataset_list_all = NULL,
                                       verbose = FALSE) {
    # for each dataset
    gset_list_all <- list()
    for (dataset_deg_list in deg_list_all) {
        # get dataset_id
        dataset_id = dataset_deg_list$dataset_id
        
        # be patient, I'm running
        .say(verbose, "Processing dataset -", dataset_id)
        
        # ident_gset_list <- list()
        
        for (ident_name in ident_names) {
            ident_deg_list = dataset_deg_list$dataset_deg_list[[ident_name]]
            # dataset_plotlist <- list()
            for (cluster_deg_list in ident_deg_list) {
                # assay_gset_list <- list()
                for (assay_name in assay_names) {
                    assay_deg_list = cluster_deg_list[[assay_name]]
                    
                    # get all info
                    # assay_name = assay_deg_list$assay_name
                    # ident_name = assay_deg_list$ident_name
                    cluster_name = assay_deg_list$cluster_name
                    
                    .say(
                        verbose,
                        "  - Processing assay:",
                        assay_name,
                        ", cluster:",
                        cluster_name
                    )
                    
                    
                    if (use_presto &
                        !is.null(multimodal_dataset_list_all)) {
                        seurat_obj <-
                            multimodal_dataset_list_all[[dataset_id]]$clusters_result$seurat_obj
                        
                        deg_df = wilcoxauc(X = as.matrix(
                            GetAssayData(
                                seurat_obj,
                                assay = assay_name,
                                slot = "data"
                            )
                        ),
                        y = seurat_obj[[ident_name]][, 1])
                        deg_df = deg_df %>%
                            dplyr::filter(padj <= max.fgsea_padj) %>%
                            arrange(desc(NES))
                        
                    } else {
                        deg_df = assay_deg_list$deg_df
                        
                        
                        # filter if needed
                        deg_df = filter_deg_df(
                            deg_df,
                            assay_name = assay_name,
                            only_positive = only_positive,
                            report_size = report_size,
                            report_size_factor = report_size_factor,
                            max.p_val_adj = max.p_val_adj,
                            min.avg_log2FC = min.avg_log2FC,
                            assays.use.report_size_factor = assays.use.report_size_factor)
                        # get gene rank
                        gene_stats = deg_df$avg_log2FC
                        names(gene_stats) = deg_df$gene
                    }
                    
                    # run fgsea
                    suppressWarnings({
                        fgseaRes <-
                            fgsea(
                                fgsea_sets,
                                stats = gene_stats,
                                scoreType = "std"
                            ) %>%
                            dplyr::filter(padj <= max.fgsea_padj) %>%
                            arrange(desc(abs(NES)))
                    })
                    
                    if (!isEmpty(fgseaRes)) {
                        fgseaRes$leadingEdge = sapply(fgseaRes$leadingEdge, function(x) {
                            paste(x, collapse = " ")
                        }, simplify = TRUE)
                        
                        g = ggplot(
                            fgseaRes %>% head(n = 20),
                            aes(reorder(pathway, NES), NES)
                        ) +
                            geom_col(aes(fill = NES)) +
                            coord_flip() +
                            labs(
                                x = "Pathway",
                                y = "Normalized Enrichment Score",
                                title = paste0(
                                    "dataset: ",
                                    dataset_id,
                                    ", assay: ",
                                    assay_name
                                )
                            ) +
                            theme_minimal()
                    } else {
                        .say(verbose,
                             "    No significantly enriched gene set found; Skip")
                        if (show_empty_plot) {
                            g = ggplot(as.data.frame(matrix())) +
                                geom_blank() +
                                scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
                                labs(title = paste0(
                                    paste0(
                                        "dataset: ",
                                        dataset_id,
                                        ", assay: ",
                                        assay_name
                                    ),
                                    "\nNo Enriched Gene Set Found"
                                )) +
                                theme_minimal()
                        } # show_empty_plot
                    } # generat plot
                    
                    gset_list_all[[ident_name]][[cluster_name]][[assay_name]][[dataset_id]] <-
                        list(
                            assay_name = assay_name,
                            ident_name = ident_name,
                            cluster_name = cluster_name,
                            dataset_id = dataset_id,
                            g = g,
                            deg_df = deg_df,
                            gset_df = fgseaRes
                        )
                } # assay_deg_list
            } # cluster_deg_list
        } # dataset_deg_list
    } # deg_list_all
    
    # return
    gset_list_all
}

write_gset_to_file_for_all <- function(gest_list_all,
                                       parent_dir = "enriched_gsets_results") {
    # create parent dir
    dir_create(TRUE, parent_dir)
    
    # write gsets and plots to file
    for (ident_name in names(gset_list_all)) {
        # ident_gset_list = gset_list_all[[ident_name]]
        
        # One folder for each identity type
        ident_dir = dir_create(TRUE, parent_dir, ident_name)
        
        for (cluster_name in names(gset_list_all[[ident_name]])) {
            # cluster_gset_list <- ident_gset_list[[cluster_name]]
            # each cluster gets a plot file and a gset folder
            plotlist <- list()
            
            cluster_gset_dir <-
                dir_create(TRUE,
                           ident_dir,
                           paste0("cluster_", cluster_name))
            
            # Now, add plot to file
            for (assay_name in names(gset_list_all[[ident_name]][[cluster_name]])) {
                # assay_gset_list <- cluster_gset_list[[assay_name]]
                
                for (dataset_id in names(gset_list_all[[ident_name]][[cluster_name]][[assay_name]])) {
                    dataset_gset_list = gset_list_all[[ident_name]][[cluster_name]][[assay_name]][[dataset_id]]
                    gset_df = dataset_gset_list$gset_df
                    
                    g = dataset_gset_list$g
                    # write enriched gene set to file
                    fn = file.path(
                        cluster_gset_dir,
                        paste0("assay_",
                            assay_name,
                            "_dataset_",
                            dataset_id,
                            ".csv"
                        )
                    )
                    if (!isEmpty(gset_df)) {
                        write.csv(gset_df,
                                  file = fn,
                                  quote = FALSE)
                        plotlist[[paste0(assay_name, "-", dataset_id)]] = g
                    }
                }
            }
            
            plot_size_factor = round(length(plotlist) / 10) + 1
            # plot
            pdf(
                file = file.path(
                    ident_dir,
                    paste0("cluster_", cluster_name, "_top_enriched_gsets.pdf")
                ),
                # The directory you want to save the file in
                width = 20 + 10 *plot_size_factor,
                # The width of the plot in inches
                height = 10 + 5 * plot_size_factor
            )
            
            if (length(plotlist) > 0) {
                print(ggarrange(plotlist = plotlist))
            }
            
            dev.off()
        } # cluster_name
    } # ident_name
}
