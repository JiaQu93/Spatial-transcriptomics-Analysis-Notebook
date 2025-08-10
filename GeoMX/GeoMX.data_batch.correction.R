
spe_lm <- geomxBatchCorrection(spe, 
                               batch = colData(spe)$sample_quality, method = "Limma",
                               design = model.matrix(~celltype, data = colData(spe)))


spe_filtered <- spe[, spe$celltype_quality %in% c("M1_1", "M2_1", "M2_2", "Tumor_1", "Tumor_2")]
spe_tumor <- spe[, spe$celltype_quality %in% c("Tumor_2")]
spe_M1 <- spe[, spe$celltype %in% c("M1")]
spe_M2 <- spe[, spe$celltype %in% c("M2")]



# Normalization and Batch correction  ================================================----

# 1. Predefined some functions for calculating silhoutte and Kbet score

#Silhoutte function: calculate the silhouette score, a measure of how well data points are clustered or separated according to batch labels.
cal_sil <- function(se.object = spe.list[[names(spe.list)[1]]],
                    #obj.name = names(spe.list)[1],
                    batch.name = unwantedVariable){
  data <- assay(se.object,i = 2)
  # find HVG
  dec <-scran:: modelGeneVar(data) #Estimates gene variance, identifying variability across samples.
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes,])
  
  # Perform PCA on HVGs and compute a distance matrix
  batch <- colData(se.object)[,batch.name]
  pca.data.f <- gmodels::fast.prcomp(data_hvg, center=TRUE)
  dd <- as.matrix(dist(pca.data.f$x[, 1:10]))
  
  # Compute silhouette score(measure how similar the distance among own batch and compare to other batch)
  batch.silhouette <- summary(cluster::silhouette(as.numeric(factor(batch,
                                                                    levels = sort(unique(batch)),
                                                                    labels = 1:length(sort(unique(batch))))), dd))$avg.width
  return(batch.silhouette) #compute the score for individual point then output the average silhouette width for entire dataset: closer to 1 means well separated cluster, closer to 0 or negative means poor separated or overlapping.
}

# kBET function (k nearest neighbor batch effect test)
cal_kbet <- function(se.object = spe.list[[names(spe.list)[1]]],
                     # obj.name = names(spe.list)[1],
                     batch.name = unwantedVariable){
  set.seed(123)
  data <- assay(se.object,i = 2)
  # find HVG
  dec <-scran:: modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes,])
  #
  batch <- colData(se.object)[,batch.name]
  k0=floor(mean(table(batch))) # neighbourhood size: mean batch size 
  knn <- get.knn(data_hvg, k=k0, algorithm = 'cover_tree') # search nearest neighbors
  batch.estimate <- kBET(data_hvg, batch, k = k0, knn = knn,plot = F) # test the nearest neighbors of each point are well distributed across the batches
  return(batch.estimate)
}


# 2. Define benchmarking parameter for RUV4 functions
unwantedVariable <- c("sample_quality") # define unwanted factor
WantedVariable <- c("Group") # define wanted factor
normalization_methods <- c( "CPM","upperquartile") # "CPM", "TMM", 
findNCG_topn <- seq(from = 2000, to = 3000, by = 500) # number of negative control gene
BatchCorrection_k <- seq(from = 3, to = 5, by = 1) # k of RUV4 (The number of unwanted factors to use. Can be 0. This is required for the RUV4 method.)
# show the number of combinations
print(paste0("Combination Parameter Number: ",length(BatchCorrection_k)*length(findNCG_topn)*length(normalization_methods)))


# 3. Run Benchmarking (run whole chunk): each parameter combinations + eachbatch -> score

# Total number of iterations/benchmarking combinations
total_iterations <- length(normalization_methods) * length(findNCG_topn) * length(BatchCorrection_k)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0

for (n_m in 1:length(normalization_methods)) { #outer loop: normalization methods
  tmp_n_m <- normalization_methods[n_m]
  tmp_spe <- geomxNorm(spe, method = tmp_n_m)
  if(tmp_n_m == "CPM"){tmp_spe <- scater::logNormCounts(spe)}
  
  for (topn in 1:length(findNCG_topn)) { #middle loop:negative control genes
    tmp_topn <- findNCG_topn[topn]
    #### here need to set up the batch factor based on meta data
    tmp_spe <- findNCGs(tmp_spe, batch_name = unwantedVariable, top_n = tmp_topn)
    
    for (k in 1:length(BatchCorrection_k)) { #inner loop: RUV4 batch correction
      tmp_k <- BatchCorrection_k[k]
      # [to do] change factors to biological variation remain based on meta data
      # here to setup wanted biology variable
      tmp_spe <- geomxBatchCorrection(tmp_spe, factors = WantedVariable,
                                      NCGs = metadata(tmp_spe)$NCGs, k = tmp_k)
      tmp_name <- paste0(tmp_n_m, "-nNCG_", tmp_topn, "-k_", tmp_k) #generate a unique name for each combination: "TMM-nNCG_100-k_3"
      spe.list <- c(spe.list, tmp_spe) #store batch corrected spe object
      names(spe.list)[length(spe.list)] <- tmp_name
      
      # Update progress bar
      iteration_count <- iteration_count + 1
      setTxtProgressBar(progress, iteration_count)
    }
  }
}

# Close progress bar
close(progress)


# 4. Evaluation using silouette and Kbet (Run whole chunk)

# Run for each of factor, including patientID (batch confounder), celltype (biological variation), and Disease (biological variation)
batch.name = c(unwantedVariable, WantedVariable)
sil_score_list <- list()
k_bet_list <- list()
# progress bar
total_iterations <- length(spe.list)*length(batch.name)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0
for (i in 1:length(spe.list)){ #outer loop for each batch corrected spe object
  se.object = spe.list[[names(spe.list)[i]]]
  obj.name = names(spe.list)[i]
  for (j in 1:length(batch.name)){ #inner loop for each batch variable calculate evaluating score
    sil_score <- cal_sil(se.object = se.object, batch.name = batch.name[j] )
    k_bet.obj <- cal_kbet(se.object = se.object, batch.name = batch.name[j] )
    sil_score_list <- c(sil_score_list, sil_score)
    k_bet_list <- c(k_bet_list, mean(k_bet.obj$stats$kBET.observed))
    names(sil_score_list)[length(sil_score_list)] <-  names(k_bet_list)[length(k_bet_list)] <- paste0(batch.name[j],"_",obj.name)
    # Update progress bar
    iteration_count <- iteration_count + 1
    setTxtProgressBar(progress, iteration_count)
  }
}#final output: sil_score and kbet score for all data sets and batch variables


# 5. Evaluation of all factors (evaluate the correction performance by analyzing sil_score and kbet score)
# 1) Minimize unwanted variance (evaluate the results for all combinations of unwanted variables(source of effect) and correction methods)
unwantedVariable_df <- matrix(rep(NA,length(unwantedVariable)*length(spe.list)*5),ncol =length(unwantedVariable)*5 )

unwantedVariable_combinations <- expand.grid(unwantedVariable,  c("_sil","_kbet","_Silrank","_kbetrank","_MeanRank"))

colnames(unwantedVariable_df) <- paste0(unwantedVariable_combinations$Var1, unwantedVariable_combinations$Var2)

rownames(unwantedVariable_df) <-  names(spe.list)
for (i in 1:length(unwantedVariable)){
  tmp_unwantedVariable <- unwantedVariable[i]
  tmp_kbet_name_vec <- grep(tmp_unwantedVariable,(names(k_bet_list)),value = F)
  tmp_kbet <-  unlist(k_bet_list[tmp_kbet_name_vec])
  tmp_sil_name_vec <- grep(tmp_unwantedVariable,(names(sil_score_list)),value = F)
  tmp_sil <-  unlist(sil_score_list[tmp_sil_name_vec])
  
  # rank scores
  names(tmp_sil) <- names(tmp_kbet) <- gsub(paste0("^",tmp_unwantedVariable,"_"),"",names(tmp_sil))
  tmp_unwantedVariable_kbet_rank <- rank(tmp_kbet)
  tmp_unwantedVariable_sil_rank <- rank(tmp_sil)
  # store metrics
  tmp_unwantedVariable_kbet_rank <- tmp_unwantedVariable_kbet_rank[rownames(unwantedVariable_df)]
  tmp_unwantedVariable_sil_rank <- tmp_unwantedVariable_sil_rank[rownames(unwantedVariable_df)]
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_sil"))] <- as.numeric(tmp_sil)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_kbet"))] <- as.numeric(tmp_kbet)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_Silrank"))] <- as.numeric(tmp_unwantedVariable_sil_rank)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_kbetrank"))] <- as.numeric(tmp_unwantedVariable_kbet_rank)
  # cal mean rank of sil and kbet
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_MeanRank"))] <- rowSums(unwantedVariable_df[,c(paste0(tmp_unwantedVariable,c("_Silrank")),
                                                                                                      paste0(tmp_unwantedVariable,c("_kbetrank")))])/2
} # the results assessing the rank across all parameter combinations to identify the best performance, which helps minimize variance

write.csv(unwantedVariable_df,file = "unwantedVariable_df.csv")


# 2) Maximize wanted variance(evaluate how well the meaningful biological variance preserved in batch correction)
WantedVariable_df <- matrix(rep(NA,length(WantedVariable)*length(spe.list)*5),ncol =length(WantedVariable)*5 )
WantedVariable_combinations <- expand.grid(WantedVariable, c("_sil","_kbet","_Silrank","_kbetrank","_MeanRank"))

colnames(WantedVariable_df) <- paste0(WantedVariable_combinations$Var1, WantedVariable_combinations$Var2)
rownames(WantedVariable_df) <-  names(spe.list)
for (i in 1:length(WantedVariable)){
  tmp_WantedVariable <- WantedVariable[i]
  tmp_kbet_name_vec <- grep(tmp_WantedVariable,(names(k_bet_list)),value = F)
  tmp_kbet <-  unlist(k_bet_list[tmp_kbet_name_vec])
  tmp_sil_name_vec <- grep(tmp_WantedVariable,(names(sil_score_list)),value = F)
  tmp_sil <-  unlist(sil_score_list[tmp_sil_name_vec])
  
  # scores ranked in descending order to maximize the wanted variance (higher score/lower rank means well separated/wanted variance)
  names(tmp_sil) <- names(tmp_kbet) <- gsub(paste0("^",tmp_WantedVariable,"_"),"",names(tmp_sil))
  tmp_WantedVariable_kbet_rank <- rank(-tmp_kbet)
  tmp_WantedVariable_sil_rank <- rank(-tmp_sil)
  # mean rank of batch
  tmp_WantedVariable_kbet_rank <- tmp_WantedVariable_kbet_rank[rownames(WantedVariable_df)]
  tmp_WantedVariable_sil_rank <- tmp_WantedVariable_sil_rank[rownames(WantedVariable_df)]
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_sil"))] <- as.numeric(tmp_sil)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_kbet"))] <- as.numeric(tmp_kbet)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_Silrank"))] <- as.numeric(tmp_WantedVariable_sil_rank)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_kbetrank"))] <- as.numeric(tmp_WantedVariable_kbet_rank)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_MeanRank"))] <- rowSums(WantedVariable_df[,c(paste0(tmp_WantedVariable,c("_Silrank")),
                                                                                                paste0(tmp_WantedVariable,c("_kbetrank")))])/2
}
write.csv(WantedVariable_df,file = "WantedVariable_df.csv")


# 6. Select your parameters (Action required for user): overall rank to identify the best performing parameter combination

# merge results of wanted variable and unwanted variable
merge_res <- cbind.data.frame(unwantedVariable_df,WantedVariable_df)
merge_res <- merge_res[,grep("_MeanRank",colnames(merge_res))]
merge_res$overall_rank <- rowSums(merge_res[,grep("_MeanRank",colnames(merge_res))])/length(merge_res[,grep("_MeanRank",colnames(merge_res))]) 
merge_res <- merge_res[order(merge_res$overall_rank),]
print(merge_res)
write.csv(merge_res,file = "merge_res.csv")





# RUV-4 batch correction (optimized parameter: CPM-nNCG_2500-k_4 )    ================================================================--------
unwantedVariable <- c("sample_quality") 
WantedVariable <- c("Group") 

tmp_spe <- geomxNorm(spe, method = "upperquartile")
# top_n = numebr of negative control gene
tmp_spe <- findNCGs(tmp_spe, batch_name = unwantedVariable, top_n = 3000)
# k is previous tuned parameter
tmp_spe <- geomxBatchCorrection(tmp_spe, factors = WantedVariable,    # Default methods is RUV4
                                NCGs = metadata(tmp_spe)$NCGs, k = 5)

#dec <-scran:: modelGeneVar(tmp_spe) # identify HVGs
#top_genes <- scran::getTopHVGs(dec, n = 1000)
# exp_mat is the post-correct expression matrix
ruv_logcounts <- assay(tmp_spe,i = 2)
# write the expression matrix
write.csv(ruv_logcounts,file = "/bmbl_data/qujia/GeoMX/CTCL/CTCL_data/CTCL_ruv_logcounts_upperquartile-nNCG_3000-k_5.csv",quote = F)

#save(tmp_spe, file = "Lymphnode_batchcorrected-TMM-nNCGH_3000-K_5.RData")




spe@colData$sample_quality <- as.character(spe@colData$sample_quality)
tmp_spe@colData$sample_quality <- as.character(tmp_spe@colData$sample_quality)


# Batch correction visualization  =========================================================================================
#1) PCA

#PCA for raw data
raw.PCA<- drawPCA(spe, assay = 2, col = Group) +
  ggtitle("PCA of raw data, group by Group")
print(raw.PCA)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/PCA/raw.PCA_Group.png", plot = raw.PCA, width = 8, height = 6, dpi = 300)

raw.PCA <- drawPCA(spe, assay = 2, col = patientID) +
  ggtitle("PCA of raw data, group by patientID")
print(raw.PCA)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/PCA/raw.PCA_patientID.png", plot = raw.PCA, width = 8, height = 6, dpi = 300)

#PCA for batch removed data
adj.PCA <- drawPCA(tmp_spe, assay = 2, col = Group) +
  ggtitle("PCA of batch corrected data, group by Group")
print(adj.PCA)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/PCA/adj.PCA_Group.png", plot = adj.PCA, width = 8, height = 6, dpi = 300)

adj.PCA <- drawPCA(tmp_spe, assay = 2, col = patientID) +
  ggtitle( "PCA of batch corrected data, group by patientID")
print(adj.PCA)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/PCA/adj.PCA_patientID.png", plot = adj.PCA, width = 8, height = 6, dpi = 300)


#2) UMAP: using Seurat to visualize the geomx sample relations (If PCA can not well separate the samples, we can use UMAP to visualize the sample relations)
setwd("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/UMAP/")
#raw data UMAP
raw.seurat <- CreateSeuratObject(spe@assays@data$logcounts)
raw.seurat <- NormalizeData(raw.seurat)
raw.seurat <- AddMetaData(raw.seurat, metadata = as.data.frame(spe@colData))
raw.seurat <- FindVariableFeatures(raw.seurat)
raw.seurat <- ScaleData(raw.seurat)
raw.seurat <- RunPCA(raw.seurat,npcs = 20)
raw.seurat <- RunUMAP(raw.seurat,dims = 1:5,reduction = "pca")
raw.seurat <- FindNeighbors(raw.seurat, dims = 1:5,reduction = "pca")  # Construct nearest-neighbor graph
#raw.seurat <- FindClusters(raw.seurat, resolution = 1)  # Assign clusters
raw.plot.df <- data.frame(raw.seurat@reductions$umap@cell.embeddings,as.data.frame(raw.seurat@meta.data))
#cluster_centroids <- raw.plot.df %>%
  #group_by(seurat_clusters) %>%  # Assuming clusters are stored in "seurat_clusters"
  #summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
# raw.plot.df[1:5,]

raw.UMAP<-ggplot(raw.plot.df, aes(x = UMAP_1, y = UMAP_2, color = patientID  )) + geom_point() +
  #geom_text_repel(data = cluster_centroids, aes(label = seurat_clusters), size = 5, color = "black") +  
  ggtitle("UMAP of raw data, group by patientID") 
print(raw.UMAP)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/UMAP/raw.UMAP_patientID.png", plot = raw.UMAP, width = 8, height = 6, dpi = 300)


raw.UMAP<-ggplot(raw.plot.df, aes(x = UMAP_1, y = UMAP_2, color = Group  )) + geom_point() +
  #geom_text_repel(data = cluster_centroids, aes(label = seurat_clusters), size = 5, color = "black") +  
  ggtitle("UMAP of raw data, group by Group") 
print(raw.UMAP)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/UMAP/raw.UMAP_Group.png", plot = raw.UMAP, width = 8, height = 6, dpi = 300)


# calculated ARI between seurat cluster and other labels
MixGHD::ARI(x=raw.plot.df$seurat_clusters, y=raw.plot.df$celltype)
MixGHD::ARI(x=raw.plot.df$seurat_clusters, y=raw.plot.df$sample_quality)
MixGHD::ARI(x=raw.plot.df$seurat_clusters, y=raw.plot.df$Best_Response)
MixGHD::ARI(x=raw.plot.df$seurat_clusters, y=raw.plot.df$Cancer_type)
MixGHD::ARI(x=raw.plot.df$seurat_clusters, y=raw.plot.df$Disease)
MixGHD::ARI(x=raw.plot.df$seurat_clusters, y=raw.plot.df$patientID)


#corrected(adj) data UMAP
adj.seurat <- CreateSeuratObject(tmp_spe@assays@data$logcounts)
adj.seurat <- AddMetaData(adj.seurat, metadata = as.data.frame(tmp_spe@colData))
adj.seurat <- FindVariableFeatures(adj.seurat)
adj.seurat <- ScaleData(adj.seurat)
adj.seurat <- RunPCA(adj.seurat,npcs = 20)
adj.seurat <- RunUMAP(adj.seurat,dims = 1:10)
adj.seurat <- FindNeighbors(adj.seurat, dims = 1:5)  # Construct nearest-neighbor graph
adj.seurat <- FindClusters(adj.seurat, resolution = 1)  # Assign clusters
adj.plot.df <- data.frame(adj.seurat@reductions$umap@cell.embeddings,as.data.frame(adj.seurat@meta.data))
cluster_centroids <- adj.plot.df %>%
  group_by(seurat_clusters) %>%  # Assuming clusters are stored in "seurat_clusters"
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

adj.UMAP<-ggplot(adj.plot.df, aes(x = UMAP_1, y = UMAP_2, color =  Group )) + geom_point() +
  ggtitle("UMAP of batch corrected data, group by Group") 
print(adj.UMAP)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/UMAP/adj.UMAP_Group.png", plot = adj.UMAP, width = 8, height = 6, dpi = 300)


adj.UMAP<-ggplot(adj.plot.df, aes(x = UMAP_1, y = UMAP_2, color =  sample_quality )) + geom_point() +
  ggtitle("UMAP of batch corrected data, group by sample_quality") 
print(adj.UMAP)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/UMAP/adj.UMAP_sample_quality.png", plot = adj.UMAP, width = 8, height = 6, dpi = 300)


adj.UMAP<-ggplot(adj.plot.df, aes(x = UMAP_1, y = UMAP_2, color =  celltype )) + geom_point() +
  ggtitle("UMAP of batch corrected data, group by celltype") 
print(adj.UMAP)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/UMAP/adj.UMAP_celltype.png", plot = adj.UMAP, width = 8, height = 6, dpi = 300)

adj.UMAP<-ggplot(adj.plot.df, aes(x = UMAP_1, y = UMAP_2, color =  celltype_response )) + geom_point() +
  ggtitle("UMAP of batch corrected data, group by celltype_response") 
print(adj.UMAP)


#3) RLE plot to check variantion
#check logCPM data: By using assay = 2 to run RLE on the logCPM data, we can see that most of the technical variations due to library size differences are removed.

raw.RLE <- plotRLExpr(spe, ordannots = "Group", assay = 2, color = Group)+
  ggtitle("RLE plot of raw data, group by Group")
print(raw.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/raw.RLE_Group.png", plot = raw.RLE, width = 10, height = 4, dpi = 300)

raw.RLE <- plotRLExpr(spe, ordannots = "patientID", assay = 2, color = patientID)+
  ggtitle("RLE plot of raw data, group by patientID")
print(raw.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/raw.RLE_patientID.png", plot = raw.RLE, width = 10, height = 4, dpi = 300)

raw.RLE <- plotRLExpr(spe, ordannots = "sample_quality", assay = 2, color = sample_quality)+
  ggtitle("RLE plot of raw data, group by sample_quality")
print(raw.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/raw.RLE_sample_quality.png", plot = raw.RLE, width = 10, height = 4, dpi = 300)



adj.RLE <- plotRLExpr(tmp_spe, ordannots = "Group", assay = 2, color = Group)+
  ggtitle("RLE plot of raw data, group by Group")
print(adj.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/adj.RLE_Group.png", plot = adj.RLE, width = 10, height = 4, dpi = 300)

adj.RLE <- plotRLExpr(tmp_spe, ordannots = "patientID", assay = 2, color = patientID)+
  ggtitle("RLE plot of raw data, group by patientID")
print(adj.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/adj.RLE_patientID.png", plot = adj.RLE, width = 10, height = 4, dpi = 300)

adj.RLE <- plotRLExpr(tmp_spe, ordannots = "sample_quality", assay = 2, color = sample_quality)+
  ggtitle("RLE plot of raw data, group by sample_quality")
print(adj.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/adj.RLE_sample_quality.png", plot = adj.RLE, width = 10, height = 4, dpi = 300)

adj.RLE <- plotRLExpr(tmp_spe, ordannots = "celltype_response", assay = 2, color = celltype_response)+
  ggtitle("RLE plot of raw data, group by celltype_response")
print(adj.RLE)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/adj.RLE_sample_quality.png", plot = adj.RLE, width = 10, height = 4, dpi = 300)



# Similarity statistics plot for evaluating batch correction ================================================================--------
#' Calculate statistics for evaluating batch correction
computeClusterEvalStats <- function(spe_object, foiColumn, precomputed = NULL,
                                    n_dimension = c(1, 2), assay = 2) {
  stopifnot(foiColumn %in% colnames(colData(spe_object)))
  
  # get PCA results
  
  # compute PCA
  if (is.null(precomputed)) {
    #calcPCA(assay(spe_object, assay), n_dimension)
    pca_object <- prcomp(t(assay(spe_object, assay)), center = TRUE, scale. = TRUE)
    pca_object <- pca_object$x[, n_dimension]
  } else {
    pca_object <- checkPrecomputedPCA(spe_object, precomputed)
  }
  
  ss <- getSilhouette(pca_object, spe_object, foiColumn)[[1]]
  
  # prepare for other stats to be computed
  k <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    dplyr::select(all_of(foiColumn))
  k <- k[, 1] |>
    unique() |>
    length()
  
  # clustering
  
  kc <- stats::kmeans(pca_object[, seq(2)], k)
  
  km_clusters <- kc$cluster
  
  types <- NULL
  
  # compute adjrand and jaccard
  df_out <- mclustcomp::mclustcomp(km_clusters, getSilhouette(pca_object, spe_object, foiColumn)[[2]]) |>
    filter(types %in% c("adjrand", "chisq", "jaccard", "overlap", "mirkin")) |>
    mutate(types = c(
      "Adjusted Rand Index", "Chi-Squared Coefficient", "Jaccard Index",
      "Overlap Coefficient", "Mirkin Distance"
    ))
  
  df_out[6, ] <- c("Silhouette Coefficient", ss) # remove 'mirkin' and 'overlap coefficient', their scores are not consistent with the other methods.
  
  return(df_out)
}



#' Compare and evaluate different batch corrected data with plotting.
plotClusterEvalStats <- function(spe_list, bio_feature_name, batch_feature_name,
                                 data_names, colors = NA) {
  #checkPackages("cluster")
  #checkPackages("mclustcomp")
  # Get stats for bio factors
  stat_bio <- lapply(spe_list, function(x) {
    computeClusterEvalStats(x, bio_feature_name)
  }) |>
    bind_rows() |>
    filter(types %in% c(
      "Adjusted Rand Index", 
      "Jaccard Index", 
      "Mirkin Distance",
      "Silhouette Coefficient"
    )) |>
    mutate(from = rep(data_names, each = 4)) # Adjust each to 4 stats
  
  # Get stats for batch factors
  stat_batch <- lapply(spe_list, function(x) {
    computeClusterEvalStats(x, batch_feature_name)
  }) |>
    bind_rows() |>
    filter(types %in% c(
      "Adjusted Rand Index", 
      "Jaccard Index", 
      "Mirkin Distance",
      "Silhouette Coefficient"
    )) |>
    mutate(from = rep(data_names, each = 4)) # Adjust each to 4 stats
  
  # Plot biology evaluation stats
  p_bio <- stat_bio |>
    mutate(from = factor(from, levels = data_names)) |>
    mutate(scores = as.numeric(scores)) |>
    ggplot(aes(from, scores, fill = types)) +
    geom_bar(stat = "identity", col = "black", width = .4) +
    facet_wrap(~types, scales = "free_y", ncol = 4) + # 4 columns for better space
    theme_bw() +
    theme(legend.position = "none") +
    xlab("Count data") +
    ylab("Scores") +
    ggtitle("Biology")
  
  # Plot batch correction evaluation stats
  p_batch <- stat_batch |>
    mutate(from = factor(from, levels = data_names)) |>
    mutate(scores = as.numeric(scores)) |>
    ggplot(aes(from, scores, fill = types)) +
    geom_bar(stat = "identity", col = "black", width = .4) +
    facet_wrap(~types, scales = "free_y", ncol = 4) + # 4 columns for better space
    theme_bw() +
    theme(legend.position = "none") +
    xlab("Count data") +
    ylab("Scores") +
    ggtitle("Batch")
  
  # Apply custom colors if provided
  if (!is.na(colors)) {
    p_bio <- p_bio + scale_fill_manual(values = colors)
    p_batch <- p_batch + scale_fill_manual(values = colors)
  }
  
  # Combine both plots
  p_bio + p_batch + patchwork::plot_layout(1, 2)
}

getSilhouette <- function(pca_object, spe, foiColumn) {
  # compute euclidean distance
  distm <- stats::dist(pca_object)
  
  # grouping
  label_by_factors <- colData(spe) |>
    as.data.frame(optional = TRUE) |>
    dplyr::select(all_of(foiColumn)) |>
    rownames_to_column() |>
    mutate(grp_id = as.numeric(factor(!!sym(foiColumn))))
  
  c <- label_by_factors$grp_id
  names(c) <- label_by_factors$rowname
  
  # calculate silhouette score
  si <- cluster::silhouette(c, distm)
  
  ss <- mean(si[, 3])
  
  return(list(ss, c))
}




set.seed(10)
spe_list <- list(spe,tmp_spe)
plotClusterEvalStats <- plotClusterEvalStats(spe_list = spe_list,
                                             bio_feature_name = "Group",
                                             batch_feature_name = "sample_quality",
                                             data_names = c("Raw","RUV4"))
print(plotClusterEvalStats)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/evaluation/plotClusterEvalStats.png", plot = plotClusterEvalStats, width = 12, height = 6, dpi = 300)


