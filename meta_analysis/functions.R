## ----parameter setting----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
columns.metadata <- c("orig.ident", "sample", "disease", "disease_duration",
                      "condition", "age", "race", "sex", 
                      "immunosuppressant", "immunosuppressant_duration",
                      "batch")


extract_cells <- function(query, reference, prefix){
  ## ----preprocessing--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query <- PercentageFeatureSet(query,
                                pattern = "^MT-",
                                col.name = "percent.mt")
  head(query@meta.data, 5)
  
  
  ## ----QC_ViolinPlot--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  options(repr.plot.width = 8, repr.plot.height = 5)
  
  VlnPlot(query,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3,
          pt.size = 0)
  ggsave(paste0(prefix, "_check_QC_vlnPlot_cell_extraction.pdf"),
         width = 15,
         height = 5)
  
  ## ----QC_subset------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query <- subset(query,
                  subset = nFeature_RNA > 200 & percent.mt < 10)
  
  ## ----sctransform_normalization--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query <- SCTransform(object = query,
                       assay = "RNA",
                       new.assay.name = "refAssay",
                       residual.features = rownames(x = reference$map),
                       reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
                       method = "glmGamPoi",
                       ncells = 2000,
                       n_genes = 2000,
                       do.correct.umi = FALSE,
                       do.scale = FALSE,
                       do.center = TRUE)
  
  
  ## ----Find anchors between query and reference-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  anchors <- FindTransferAnchors(reference = reference$map,
                                 query = query,
                                 k.filter = NA,
                                 reference.neighbors = "refdr.annoy.neighbors",
                                 reference.assay = "refAssay",
                                 query.assay = "refAssay",
                                 reference.reduction = "refDR",
                                 normalization.method = "SCT",
                                 features = intersect(rownames(x = reference$map),
                                                      VariableFeatures(object = query)),
                                 dims = 1:50,
                                 n.trees = 20,
                                 mapping.score.k = 100)
  
  
  ## ----Transfer cell type labels and impute protein expression--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  refdata <- lapply(X = "celltype.l1",
                    function(x) {
                      reference$map[[x, drop = TRUE]]
                      })
  names(x = refdata) <- "celltype.l1"
  if (TRUE) {
    refdata[["impADT"]] <- GetAssayData(
      object = reference$map[['ADT']],
      slot = 'data'
    )
  }
  
  query <- TransferData(reference = reference$map,
                        query = query,
                        dims = 1:50,
                        anchorset = anchors,
                        refdata = refdata,
                        n.trees = 20,
                        store.weights = TRUE)
  
  
  ## ----Calculate the embeddings of the query data on the reference SPCA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query <- IntegrateEmbeddings(anchorset = anchors,
                               reference = reference$map,
                               query = query,
                               reductions = "pcaproject",
                               reuse.weights.matrix = TRUE)
  
  
  ## ----Calculate the query neighbors----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query[["query_ref.nn"]] <- FindNeighbors(object = Embeddings(reference$map[["refDR"]]),
                                           query = Embeddings(query[["integrated_dr"]]),
                                           return.neighbor = TRUE,
                                           l2.norm = TRUE)
  
  
  ## ----NNTransform----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query <- Azimuth:::NNTransform(object = query,
                                 meta.data = reference$map[[]])
  
  
  ## ----Project the query to the reference UMAP------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query[["proj.umap"]] <- RunUMAP(object = query[["query_ref.nn"]],
                                  reduction.model = reference$map[["refUMAP"]],
                                  reduction.key = 'UMAP_')
  
  
  
  ## ----Calculate mapping score and add to metadata--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  query <- AddMetaData(object = query,
                       metadata = MappingScore(anchors = anchors),
                       col.name = "mapping.score")
  
  
  ## ----predicted metadata field---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  id <- "celltype.l1"
  predicted.id <- paste0("predicted.", id)
  
  # DimPlot of the reference
  #DimPlot(object = reference$plot,
  #        reduction = "refUMAP",
  #        group.by = id,
  #        label = TRUE) + NoLegend()
  
  # DimPlot of the query, colored by predicted cell type
  DimPlot(object = query,
          reduction = "proj.umap",
          group.by = predicted.id,
          label = TRUE) + NoLegend()
  ggsave(paste0(prefix, "_check_dimplot_cell_extraction.pdf"),
         width = 5,
         height = 5)
  
  # Plot the score for the predicted cell type of the query
  #FeaturePlot(object = query,
  #            features = paste0(predicted.id, ".score"),
  #            reduction = "proj.umap")
  #VlnPlot(object = query,
  #        features = paste0(predicted.id, ".score"),
  #        group.by = predicted.id,
  #        pt.size = 0) + NoLegend()
  
  # Plot the mapping score
  #FeaturePlot(object = query,
  #            features = "mapping.score",
  #            reduction = "proj.umap")
  #VlnPlot(object = query,
  #        features = "mapping.score",
  #        group.by = predicted.id,
  #        pt.size = 0) + NoLegend()
  
  genes <- c("CD3E", "CD4", "CD8A", "FOXP3", "IL7R", "CD74")
  
  ## ----DotPlot for finding CD4+Tcell----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  options(repr.plot.width = 10, repr.plot.height = 8)
  
  DotPlot(object = query,
          features = genes,
          group.by = predicted.id) + RotatedAxis()
  ggsave(paste0(prefix, "_check_dotplot_CD4T_extraction.pdf"),
         width = 10,
         height = 10)
  
  ## ----FeaturePlot for finding CD4+Tcell------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  options(repr.plot.width = 10, repr.plot.height = 8)
  
  FeaturePlot(object = query,
              features = genes,
              reduction = "proj.umap",
              ncol = 2,
              raster = TRUE)
  
  ggsave(paste0(prefix, "_check_featureplot_CD4T_extraction.pdf"),
         width = 10,
         height = 15)
  
  
  genes <- c("CD79A", "MS4A1")
  ## ----DotPlot for finding Bcell----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  options(repr.plot.width = 10, repr.plot.height = 8)
  
  DotPlot(object = query,
          features = genes,
          group.by = predicted.id) + RotatedAxis()
  ggsave(paste0(prefix, "_check_dotplot_B_extraction.pdf"),
         width = 10,
         height = 10)
  
  ## ----FeaturePlot for finding Bcell------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  options(repr.plot.width = 10, repr.plot.height = 8)
  
  FeaturePlot(object = query,
              features = genes,
              reduction = "proj.umap",
              ncol = 2,
              raster = TRUE)
  
  ggsave(paste0(prefix, "_check_featureplot_B_extraction.pdf"),
         width = 10,
         height = 5)
  
  
  ## ----Extraction and save of CD4+Tcell, Bcell ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (length(query@meta.data[query@meta.data['predicted.celltype.l1'] == "CD4 T"]) > 0 ) {
    pbmc_CD4T <- subset(query,
                        predicted.celltype.l1 == "CD4 T")
    
    pbmc_CD4T.meta.data <- pbmc_CD4T@meta.data
    pbmc_CD4T.meta.data <- pbmc_CD4T.meta.data %>%
      select(all_of(columns.metadata))
    
    pbmc_CD4T <- GetAssayData(object = pbmc_CD4T,
                              slot = "counts")
    
    saveRDS(pbmc_CD4T.meta.data,
            file = paste0(prefix, "_CD4T_MetaData.rds"))
    
    saveRDS(pbmc_CD4T,
            file = paste0(prefix, "_CD4T_AssayData.rds"))
  } else {
    print("No CD4T cell was detected in the sample.")
  }
  
  if (length(query@meta.data[query@meta.data['predicted.celltype.l1'] == "B"]) > 0 ) {
    pbmc_B <- subset(query,
                        predicted.celltype.l1 == "B")
    
    pbmc_B.meta.data <- pbmc_B@meta.data
    pbmc_B.meta.data <- pbmc_B.meta.data %>%
      select(all_of(columns.metadata))
    
    pbmc_B <- GetAssayData(object = pbmc_B,
                           slot = "counts")
    
    saveRDS(pbmc_B.meta.data,
            file = paste0(prefix, "_B_MetaData.rds"))
    
    saveRDS(pbmc_B,
            file = paste0(prefix, "_B_AssayData.rds"))
  } else {
    print("No B cell was detected in the sample.")
  }

}





reference_mapping <- function(ref, query_obj, prefix){
  ## ----parameter setting----------------------------------------------------------------------------
  vars.obj <- c("sample")
  vars.query1 <- c("batch")
  vars.query2 <- c("disease")
  
  genes <- c("CD3E", "MS4A1", "CST3", "CD4", "CD8A",
             "TBX21", "GATA3", "RORC", "FOXP3", "CCR8",
             "IL2RA", "S1PR1", "CCR7", "IL7R", "FAS",
             "CD28", "CTLA4", "IKZF4", "NRP1", "BCL6",
             "CXCR5", "TOX", "CSF2", "MX1")
  
  c1 <- c("sample")
  c2 <- c("sex")
  c3 <- c("project_name", "sample", "disease", "disease_duration",
          "condition", "age", "race", "sex", "immunosuppressant", "immunosuppressant_duration",
          "batch")
  treg <- c("Treg")
  
  
  ## ----query_obj--------------------------------------------------------------------------
  query_obj <- CreateSeuratObject(query_obj) %>%
    SCTransform(assay = "RNA",
                method = "glmGamPoi",
                return.only.var.genes = FALSE)
  
  query_obj.meta.data <- readRDS(paste0(prefix, "_CD4T_MetaData.rds"))
  
  query_obj <- AddMetaData(query_obj,
                           metadata = query_obj.meta.data)
  
  
  ## ----Map Query------------------------------------------------------------------------------------
  query <- mapQuery(exp_query = query_obj$SCT@scale.data,
                    metadata_query = query_obj@meta.data,
                    ref_obj = ref,
                    vars = vars.query1,
                    do_normalize = FALSE,
                    return_type = "Seurat")
  
  
  ## ----UMAP-----------------------------------------------------------------------------------------
  options(repr.plot.height = 4, repr.plot.width = 10)
  
    (
      DimPlot(query,
              reduction = "umap",
              group.by = vars.query1,
              shuffle = TRUE) +
        labs(title = paste0("Mapped Query (", vars.query1, ")"))
      ) +
    (
      DimPlot(query,
              reduction = "umap",
              group.by = vars.query2,
              shuffle = TRUE) +
        labs(title = paste0("Mapped Query (", vars.query2, ")"))
      )
  
  ggsave(paste0(prefix, "_mapped_query_vars_batch_disease_Reference_Mapping.pdf"),
         width = 10,
         height = 5)
  
  
  ## ----Predict clusters-----------------------------------------------------------------------------
  queryL1 <- knnPredict.Seurat(query_obj = query,
                               ref_obj = ref,
                               label_transfer = "clusterL1")
  
  queryL2 <- knnPredict.Seurat(query_obj = query,
                               ref_obj = ref,
                               label_transfer = "clusterL2")
  
  
  ## ----UMAP_predict_clusters------------------------------------------------------------------------
  options(repr.plot.height = 4, repr.plot.width = 10)
  
  (
    DimPlot(queryL1,
            reduction = "umap",
            group.by = "clusterL1",
            shuffle = TRUE) +
      labs(title = "Mapped Query (Predicted clusterL1)")
    ) +
    (
      DimPlot(queryL2,
              reduction = "umap",
              group.by = "clusterL2",
              shuffle = TRUE) +
        labs(title = "Mapped Query (Predicted clusterL2)")
    )
  
  ggsave(paste0(prefix, "_predict_clusterL1L2_Reference_Mapping.pdf"),
         width = 10,
         height = 5)
  
  
  ## ----FeaturePlot for Reference Mapping------------------------------------------------------------
  options(repr.plot.width = 10, repr.plot.height = 8)
  
  FeaturePlot(object = queryL1,
              features = genes,
              reduction = "umap",
              ncol = 2,
              raster = TRUE)
  
  ggsave(paste0(prefix, "_featureplot_queryL1_Reference_Mapping.pdf"),
         width = 10,
         height = 48)
  
  
  FeaturePlot(object = queryL2,
              features = genes,
              reduction = "umap",
              ncol = 2,
              raster = TRUE)
  
  ggsave(paste0(prefix, "_featureplot_queryL2_Reference_Mapping.pdf"),
         width = 10,
         height = 48)
  
  
  ## ----output files---------------------------------------------------------------------------------

  save(queryL1,
       file = paste0(prefix, "_queryL1_Reference_Mapping.RData"))
  save(queryL2,
       file = paste0(prefix, "_queryL2_Reference_Mapping.RData"))
  
  
  
  queryL1@meta.data <- queryL1@meta.data %>%
    rename(project_name = orig.ident)
  
  queryL1.meta.data <- queryL1@meta.data
  
  queryL1_clusterL1 <- queryL1.meta.data %>%
    group_by(.data[[c1]], clusterL1) %>%
    count()
  queryL1_clusterL1_prob <- queryL1.meta.data %>%
    group_by(.data[[c1]], clusterL1) %>%
    summarise(clusterL1_prob.mean = mean(clusterL1_prob),
              clusterL1_prob.max = max(clusterL1_prob),
              clusterL1_prob.min = min(clusterL1_prob))
  queryL1_clusterL1 <- left_join(queryL1_clusterL1,
                                 queryL1_clusterL1_prob,
                                 by = c(c1, "clusterL1"))  
  write.csv(queryL1_clusterL1,
            file = paste0(prefix, "_queryL1_Reference_Mapping.csv"),
            row.names = FALSE,
            quote = FALSE)
  
  queryL1.meta.data <- queryL1.meta.data %>%
    select(all_of(c3)) %>%
    distinct(.data[[c1]],
             .keep_all = TRUE)
  
  write.csv(queryL1.meta.data,
            file = paste0(prefix, "_queryL1_metadata_Reference_Mapping.csv"),
            row.names = FALSE,
            quote = FALSE)
  
  


  queryL2@meta.data <- queryL2@meta.data %>%
    rename(project_name = orig.ident)
  
  queryL2.meta.data <- queryL2@meta.data
  
  queryL2_clusterL2 <- queryL2.meta.data %>%
    group_by(.data[[c1]], clusterL2) %>%
    count()
  queryL2_clusterL2_prob <- queryL2.meta.data %>%
    group_by(.data[[c1]], clusterL2) %>%
    summarise(clusterL2_prob.mean = mean(clusterL2_prob),
              clusterL2_prob.max = max(clusterL2_prob),
              clusterL2_prob.min = min(clusterL2_prob))
  queryL2_clusterL2 <- left_join(queryL2_clusterL2,
                                 queryL2_clusterL2_prob,
                                 by = c(c1, "clusterL2"))  
  write.csv(queryL2_clusterL2,
            file = paste0(prefix, "_queryL2_Reference_Mapping.csv"),
            row.names = FALSE,
            quote = FALSE)
  
  queryL2.meta.data <- queryL2.meta.data %>%
    select(all_of(c3)) %>%
    distinct(.data[[c1]],
             .keep_all = TRUE)
  
  write.csv(queryL2.meta.data,
            file = paste0(prefix, "_queryL2_metadata_Reference_Mapping.csv"),
            row.names = FALSE,
            quote = FALSE)
    
  
  ## ----Extraction of Treg---------------------------------------------------------------------------
  query_Treg <- subset(queryL1,
                       clusterL1 %in% treg)
  
  query_Treg.meta.data <- query_Treg@meta.data
  
  query_Treg.meta.data <- query_Treg.meta.data %>%
    select(all_of(c3))
  
  query_Treg <- GetAssayData(object = query_Treg,
                               slot = "counts")
  
  
  ## ----saveRDS--------------------------------------------------------------------------------------
  saveRDS(query_Treg.meta.data,
          file = paste0(prefix, "_Treg_MetaData.rds"))
  
  saveRDS(query_Treg,
          file = paste0(prefix, "_Treg_AssayData.rds"))
  
}


