# The ‘standR’ package can not be installed in OSC, so run the pipeline locally.
# This is a GeoMX data analysis pipeline: including load data, construct the GeoMX (Spatial Experiment) object, check data, QC, normalization & batch correction, differential expression analysis, and downstream analysis


# Loading packages   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("standR")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ExperimentHub")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("kBET")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msigdb")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva")

install_github('theislab/kBET')
install.packages('devtools')
install.packages('tidyverse')
install.packages('readxl')
install.packages('msigdb')
install.packages("ggvenn")

# Load environment   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
library(devtools)
library(standR)
library(SpatialExperiment)
library(limma)
library(edgeR)
library(ExperimentHub)
library(ggalluvial)
require(FNN)
library(kBET)
library(Seurat)
library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(msigdb)
library(devtools)
library(ggpubr)
library(sva)
library(scater)
library(enrichR)
library(ggvenn)
library(gridExtra)
library(grid)
library(msigdbr)
library(fgsea)

set.seed(100)

# Load data  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setwd("/bmbl_data/qujia/GeoMX/CTCL")


# Construct the GeoMX (Spatial Experiment) object   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sampleAnnoFile <- read_excel('/bmbl_data/qujia/GeoMX/CTCL/CTCL_data/ProbeQC.xlsx',sheet = "SegmentProperties")  %>% as.data.frame()
countFile <- read_excel('/bmbl_data/qujia/GeoMX/CTCL/CTCL_data/ProbeQC.xlsx',sheet = "TargetCountMatrix")  %>% as.data.frame()
#row.names(countFile) <- countFile[,1]
#countFile<-countFile[,-1]

featureAnnoFile <- read_excel(
  "/bmbl_data/qujia/GeoMX/CTCL/CTCL_data/ProbeQC.xlsx",
  sheet = "TargetProperties"
) %>% as.data.frame()


CTCL_metadata<-read.csv("/bmbl_data/qujia/GeoMX/CTCL/CTCL_data/CTCL_metadata.csv") # update metadata in sampleAnnoFile
sampleAnnoFile <- sampleAnnoFile %>%
  left_join(CTCL_metadata, by = c("SegmentDisplayName" = "sample.name"))

spe <- readGeoMx(countFile, sampleAnnoFile, featureAnnoFile,
                 colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),)
spe.list <- list()
spe@colData$patientID <- as.character(spe@colData$patientID)
spe@colData$Disease <- as.character(spe@colData$Disease)
#spe@colData$sample_quality <- as.character(spe@colData$sample_quality) (dont set as.character for RUV4 bacth correction)

#tmp_spe_tumor <- tmp_spe[, tmp_spe$celltype == "Tumor"]
#spe_M1 <- spe[, spe$celltype == "M1"]
#spe_M2 <- spe[, spe$celltype == "M2"]

# Note 1: By default, the readGeoMx function will look for the gene name column in both the countFile and featureAnnoFile 
# with the column name of "TargetName", and the sample name column in the sampleAnnoFile with the column name of "SegmentDisplayName", 
# these column names are given by the Nanostring in the default settings, if your data have been modified, you can indicate the 
# corresponding column names by specifying the parameter "colnames.as.rownames" in the readGeoMx function when loading the data.



# Check data    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# check count table
assayNames(spe)
assay(spe, "counts")[1:5, 1:5]
assay(spe, "logcounts")[1:5, 1:5]
# check sampel metadata
colData(spe)[1:5, 1:5]
# check gene meta date
rowData(spe)[1:5, 1:5]
# check negtive probes
metadata(spe)$NegProbes[, 1:5]
# check low quality tissue samples
colData(spe)$QCFlags
#row.names(countFile) <- countFile[,1]
#countFile<-countFile[,-1]


# QC    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#sample level QC
#plotSampleInfo(spe, column2plot = c("patientID","SegmentLabel"))

#gene level QC
# spe <- addPerROIQC(spe, rm_genes = TRUE) #sample_fraction = 0.9 # nolint
spe <-addPerROIQC(
  spe,
  sample_fraction = 0.95,
  rm_genes = TRUE,
  min_count = 5,
  design = NULL
)
dim(spe) # 17783   135
# plotGeneQC(spe, ordannots = "SegmentLabel", col = SegmentLabel, point_size = 2) 


#ROI level QC:  identify (if any) ROI(s) that have relatively low library size and/or low cell count because they are considered 
#as low quality samples due to insufficient sequencing depth or the lack of RNA in the selected region. 
# plotROIQC(spe, x_threshold = 10, color = patientID)
# 1) filter AOINucleiCount
qc <- colData(spe)$AOINucleiCount > 10 # remove 29 ROIs
table(qc)
spe <- spe[, qc]
# 2) filter sample_quality
spe<-spe[,spe$sample_quality == 1 | spe$sample_quality ==2]

dim(spe) # 17783   100



# Check technical variation (using relative log expression(RLE) distribution)
plotRLExpr(spe) #check the unfiltering count data
#check logCPM data: By using assay = 2 to run RLE on the logCPM data, we can see that most of the technical variations due to library size differences are removed.
RLEplot_patientID<-plotRLExpr(spe, ordannots = "patientID", assay = 2, color = patientID)
RLEplot_patientID
RLEplot<-plotRLExpr(spe, ordannots = "Group", assay = 2, color = Group)
RLEplot

ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/RLEplot_patientID.png", plot = RLEplot_patientID, width = 10, height = 4, dpi = 300)
ggsave("/bmbl_data/qujia/GeoMX/CTCL/batch_correction/RLEplot/RLEplot_Group.png", plot = RLEplot, width = 10, height = 4, dpi = 300)


# Check systemic variation (using dimension reduction: PCA and UMAP)
#pre-compute the PCA to keep results consistent
spe <- scater::runPCA(spe)
pca_results <- reducedDim(spe, "PCA")
pca_patientID<-drawPCA(spe, precomputed = pca_results, col = patientID)
print(pca_patientID)
pca_Group<-drawPCA(spe, precomputed = pca_results, col = Group)
print(pca_Group)



# Normalization and Batch correction (CPM.nNCG_2500.k_3)(Identify the optimized parameters for normalization by file "GeoMX.data_batch.correction.R")   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
unwantedVariable <- c("sample_quality") 
WantedVariable <- c("Best_Response") 

tmp_spe <- geomxNorm(spe, method = "upperquartile")
# tmp_spe <- scater::logNormCounts(spe)
# top_n = numebr of negative control gene
tmp_spe <- findNCGs(tmp_spe, batch_name = unwantedVariable, top_n = 3000)
# k is previous tuned parameter
tmp_spe <- geomxBatchCorrection(tmp_spe, factors = WantedVariable,
                                NCGs = metadata(tmp_spe)$NCGs, k = 5)

#dec <-scran:: modelGeneVar(tmp_spe) # identify HVGs
#top_genes <- scran::getTopHVGs(dec, n = 1000)
# exp_mat is the post-correct expression matrix
ruv_logcounts <- assay(tmp_spe,i = 2)
# write the expression matrix
write.csv(
  ruv_logcounts,
  file = "/bmbl_data/qujia/GeoMX/CTCL/CTCL_data/CTCL_ruv_logcounts_upper_5000_5_patientID_Group.csv", 
  quote = F
)

#save(tmp_spe, file = "Lymphnode_batchcorrected-TMM-nNCGH_3000-K_5.RData")


# Batch correction visualization  --------------------------------------------------------------------------------




# Differential expression analysis with limma-voom pipeline   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

setwd('/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results')
# check the weight matrix used in next linear model(the weight matrix generated from batch correction should be included as covariate)
colData(tmp_spe)[,seq(ncol(colData(tmp_spe))-1, ncol(colData(tmp_spe)))] |>
  head()

# Extract the count matrix and coldata
spe_tumor<-tmp_spe[,spe$celltype == "Tumor"] # set the spe_obj by user
counts_matrix <- assay(spe_obj,i = 2)
columndata <- colData(spe_obj)
group = spe_obj$Type
design_data <- data.frame(group = group, 
                          columndata[,grep("ruv_",colnames(columndata))],
                          CT_count = spe_obj$AOINucleiCount)
dge <- DGEList(counts = counts_matrix, group = group)


# Create DGEList object to incorperate the limma-voom ------------------RUV4_DEG
# establish a design matrix and contrast
design <- model.matrix(~0 + group +  ., data = design_data) 
colnames(design)

colnames(design)<-gsub('^group','',colnames(design)) # edit the colnames
colnames(design)<-gsub(' ','_',colnames(design))
colnames(design)

contr.matrix<-makeContrasts(BvT = B_cell_zone - T_cell_zone, levels= colnames(design)) # make comparison between B cell and T cell



# Differential Expression
#v<- voom(dge, design, plot=TRUE)
#fit <- lmFit(v)
fit <- lmFit(dge$counts, design) # dge$counts has been normalized and btch corrected, do not need voom transformed
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit_ruv <- eBayes(fit_contrast, robust = TRUE)
results_efit_ruv<- decideTests(efit_ruv, p.value = 0.05)
summary_efit_ruv <- summary(results_efit_ruv)
summary_efit_ruv  #   BvT: Down 1510 NotSig 16054 Up 1448

# DEG results: UP, DOWN, NOT DE
deg_result <- topTable(efit_ruv, coef = 1, sort.by = "P", n = Inf) %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val < 0.05, "UP", 
                     ifelse(logFC < 0 & adj.P.Val < 0.05, "DOWN", "NOT DE")))
write.csv(deg_result,file = "/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/deg_result_BvT.csv.csv")


# DEG list :sig DEGs; p.value = 0.05 (cutoff value for adjusted p-values)
deg_list <- topTable(efit_ruv, coef = 1, sort.by = "P", n = Inf, p.value = 0.01) %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val < 0.05, "UP", 
                     ifelse(logFC < 0 & adj.P.Val < 0.05, "DOWN", "NOT DE")))
write.csv(deg_list,file = "/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/BvT_DEG_sva.csv")

# make interactive table --------------------------------
library(DT) 
updn_cols <- c(RColorBrewer::brewer.pal(6, 'Greens')[2], RColorBrewer::brewer.pal(6, 'Purples')[2])
de_genes_top_BvT %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val", "DE")) %>%
  DT::datatable(caption = 'B cell zone vs. T cell zone in Lymph node (limma-voom)') %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)


# Load DEGs list
deg_result_BvT<-read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/deg_result_BvT.csv",row.names=1, check.names = FALSE)
deg_list_BvT<-read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/deg_list_BvT.csv",row.names=1, check.names = FALSE)



# Venn plot of DEGs   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
BvT_DEG_cb <- read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/BvT_DEG_cb.csv")
BvT_DEG_lm <- read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/BvT_DEG_lm.csv")
BvT_DEG_raw <- read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/BvT_DEG_raw.csv")
BvT_DEG_ruv <- read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/BvT_DEG_ruv.csv")
BvT_DEG_sva <- read.csv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/BvT_DEG_sva.csv")

deg_list <- list(
  "ComBat(387)"  = BvT_DEG_cb[,1],
  "Limma(2453)"  = BvT_DEG_lm[,1],
  "RAW(1682)" = BvT_DEG_raw[,1],
  "RUV-4(2958)" = BvT_DEG_ruv[,1],
  "SVA(2504)" = BvT_DEG_sva[,1]
)

overlap_results <- calculate.overlap(x = deg_list) #get the overlap and unique values of venn diagram
RUV_unique<-overlap_results$a4 #extract RUV-4 unique genes

RUV_unique<- as.data.frame(RUV_unique)
write.csv(RUV_unique,  file ="/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/RUV_unique_deg.csv")


venn_plot <- ggVennDiagram(deg_list[1:5], label_alpha = 0 , label = "count", 
                           color =  c("ComBat(387)" = "orange","Limma(2453)" ="steelblue",'RAW(1682)' = 'black', 'SVA(2504)' = 'red', "SVA" = "purple"),
                           set_color =  c("ComBat(387)" = "orange","Limma(2453)" ="steelblue",'RAW(1682)' = 'black', 'SVA(2504)' = 'red', "SVA" = "purple")) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  scale_color_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "right")


print(venn_plot)

ggsave("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/Venn_DEGs.png", 
       plot = venn_plot, width = 14, height = 8, dpi = 300)


# DEG visualization: MA plot
MA_plot_BvT<-deg_result_BvT %>%
  ggplot(aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = deg_result_BvT %>%
                    rownames_to_column()%>%
                    filter(abs(logFC) > 1 | adj.P.Val < 0.01), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("B cell zone vs. T cell zone in Lymph node (limma-voom)") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme(text = element_text(size = 15))
print(MA_plot_BvT)
ggsave(MA_plot_BvT , 
       filename = "/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/MA_plot_BvT.jpg",
       device = "jpg", 
       width = 6, 
       height = 10, 
       units = "in",
       dpi = 150)



# DEG visualization: volcano plot
#deg_result_BvT$genenames<-row.names(deg_result_BvT)
volcano_plot_BvT<- ggplot(deg_result_BvT,aes(x=logFC,y=-log10(adj.P.Val)))+
  geom_point(aes(color=DE), shape = 1, size = 1)+
  scale_color_manual(values=c("dodgerblue","grey","firebrick"))+
  geom_text_repel(data = deg_result_BvT %>%
                    rownames_to_column()%>%
                    filter(abs(logFC) > 2 | adj.P.Val < 0.01), aes(label = rowname),
                  box.padding = unit(2, "lines"),
                  point.padding = unit(0.5, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)+
  theme_bw()+
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) + #add vertical line
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)  #add horizontal line

print(volcano_plot_BvT)
ggsave(volcano_plot_BvT, 
       filename = "/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/volcano_plot_BvT.jpg",
       device = "jpg", 
       width = 6, 
       height = 10, 
       units = "in",
       dpi = 150)


# Gene Set Enrichment Analysis (using enrichR) --------------------------------------------------------------------------------
dbs <- c("Reactome_2022","KEGG_2021_Human")  #select pathway databases
# "GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023"

# UP regulated DEG with adj.p.value < 0.05 and log2foldchange > 1.5 for enrichment analysis
up_pathway <- enrichr(rownames(deg_list_BvT[which(deg_list_BvT$adj.P.Val < 0.05 & deg_list_BvT$logFC > 1.5),]),dbs)
up_pathway <- lapply(up_pathway, function(df) {
  df[df$Adjusted.P.value < 0.05, ]  # Filter significant pathways
})

names(up_pathway) # check enriched pathway for each database
up_KEGG <-up_pathway[["KEGG_2021_Human"]] # 35 enriched pathways
up_Reactome <-up_pathway[["Reactome_2022"]] # 43 enriched pathways

#print the top 20 enriched terms for GO_Biological_Process_2023
DT::datatable(head(up_enriched_combined$GO_Biological_Process_2023,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(pageLength = 5,scrollX = TRUE,scrollCollapse = TRUE,dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')
              ))

#write each table to a separate CSV file
for (db in names(up_enriched_combined)) {
  write.csv(up_enriched_combined[[db]], 
            file = file.path("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/UP_enrichment", 
                             paste0(db, "_enrichment_results_up.csv")), 
            row.names = FALSE)
}

# bar-plot for enrichment analysis
up_Reactome_2022 <- up_Reactome_2022 %>%arrange(desc(-log10(Adjusted.P.value))) %>%slice(1:20) # Select top 20 terms
up_Reactome_2022_barplot<-ggplot(up_Reactome_2022, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(x = "Reactome Pathways", y = "-log10(Adjusted.P.value)", title = "Top Reactome Enrichment Results") +
  theme_minimal()
ggsave(up_Reactome_2022_barplot , 
       filename = "/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/UP_enrichment/up_enrichment_plot/up_Reactome_2022_barplot.jpg",
       device = "jpg", width = 12, height = 6, units = "in",dpi = 150)

# DOWN regulated DEG with adj.p.value < 0.05 and log2foldchange < -1.5 for enrichment analysis
down_pathway <- enrichr(rownames(deg_list_BvT[which(deg_list_BvT$adj.P.Val < 0.05 & deg_list_BvT$logFC < -1.5),]),dbs)
down_pathway <- lapply(down_pathway, function(df) {
  df[df$Adjusted.P.value < 0.05, ]  # Filter significant pathways
})

names(down_pathway) # check enriched pathway for each database
down_KEGG <-down_pathway[["KEGG_2021_Human"]] # 20 enriched pathways
down_Reactome <-down_pathway[["Reactome_2022"]] # 26 enriched pathways





#print the top 20 enriched terms for GO_Biological_Process_2023
DT::datatable(head(down_enriched_combined$GO_Biological_Process_2023,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(pageLength = 5, scrollX = TRUE,scrollCollapse = TRUE,dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')
              ))

#write each table to a separate CSV file
for (db in names(down_enriched_combined)) {
  write.csv(down_enriched_combined[[db]], 
            file = file.path("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/DOWN_enrichment", 
                             paste0(db, "_enrichment_results_down.csv")), 
            row.names = FALSE)
}

#bar-plot for enrichment analysis
down_Reactome_2022 <- down_Reactome_2022 %>%arrange(desc(-log10(Adjusted.P.value))) %>%slice(1:20) # Select top 20 terms

down_Reactome_2022_barplot<-ggplot(down_Reactome_2022, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(x = "Reactome Pathways", y = "-log10(Adjusted.P.value)", title = "Top Reactome Enrichment Results") +
  theme_minimal()

ggsave(down_Reactome_2022_barplot , 
       filename = "/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/UP_enrichment/up_enrichment_plot/down_Reactome_2022_barplot.jpg",
       device = "jpg", 
       width = 12, 
       height = 6, 
       units = "in",
       dpi = 150)



# Combined analysis of up- and down-regulated genes with adj.p.value < 0.05 and abs(log2foldchange) > 1.5 for enrichment analysis
enriched_combined <- enrichr(rownames(de_results_BvT_with_DE[which(de_results_BvT_with_DE$adj.P.Val < 0.05 & abs(de_results_BvT_with_DE$logFC) > 1.5),]),dbs)

names(enriched_combined) # check enriched pathway for each database
down_GO_Biological_Process_2023 <-enriched_combined[["GO_Biological_Process_2023"]]
down_Reactome_2022 <-enriched_combined[["Reactome_2022"]]

#print the top 20 enriched terms for GO_Biological_Process_2023
DT::datatable(head(enriched_combined$GO_Biological_Process_2023,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(pageLength = 5, scrollX = TRUE,scrollCollapse = TRUE,dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel')
              ))


#write each table to a separate CSV file
for (db in names(enriched_combined)) {
  write.csv(enriched_combined[[db]], 
            file = file.path("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/DOWN_enrichment", 
                             paste0(db, "_enrichment_results_down.csv")), 
            row.names = FALSE)}


# Gene Set Enrichment Analysis (using GSEA)  --- All DEG:  pathway for all pathway database ==========================================================
# define pathway database 
#mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C5")
#mdb_kegg = mdb_c2 [grep("^KEGG",mdb_c2 $gs_name),]
#fgsea_sets<- mdb_kegg %>% split(x = .$gene_symbol, f = .$gs_name)


setwd("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA")
if(!file.exists("../all_gene_sets_human.qsave")) {
  all_gene_sets = msigdbr(species = "human")
  qs::qsave(all_gene_sets, "all_gene_sets_human.qsave")
  
} else {
  all_gene_sets <- qs::qread("../all_gene_sets_human.qsave")
}

fgsea_sets<- all_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

# DEG list
BvT_DEG_sva$genes = BvT_DEG_sva$X
deg_list<-BvT_DEG_sva
deg_genes<- deg_list %>% arrange(desc(logFC)) %>% dplyr::select(genes,logFC) #reoder DEGs by logfc
ranks<- deframe(deg_genes)

# run fgsea
fgsea<- fgsea(fgsea_sets, stats = ranks) # , minSize=15, maxSize=500, nperm=100000)
fgsea_list_sva<- fgsea[fgsea$padj < 0.05,]
library(data.table)

fwrite(fgsea_list_sva, file="/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/fgsea_list_sva.tsv", sep="\t", sep2=c("", " ", ""))
fgsea_list_raw<-readr::read_tsv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/fgsea_list_raw.tsv")
fgsea_list_cb<-readr::read_tsv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/fgsea_list_cb.tsv")
fgsea_list_lm<-readr::read_tsv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/fgsea_list_lm.tsv")
fgsea_list_ruv<-readr::read_tsv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/fgsea_list_ruv.tsv")
fgsea_list_sva<-readr::read_tsv("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/fgsea_list_sva.tsv")


# gsea plot
REACTOME_METALLOPROTEASE_DUBS<-plotEnrichment(fgsea_sets[["REACTOME_METALLOPROTEASE_DUBS"]],
                                              ranks) + labs(title="REACTOME_METALLOPROTEASE_DUBS")

REACTOME_METALLOPROTEASE_DUBS


# GSEA venn diagram 

 gsea_list <- list(
    "ComBat(76)"  = fgsea_list_cb[[1]], # extract characters form list
    "Limma(939)"  = fgsea_list_lm[[1]],
    "RAW(797)" = fgsea_list_raw[[1]],
    "RUV(1036)" = fgsea_list_ruv[[1]],
    "SVA(985)" = fgsea_list_sva[[1]]
  )

 overlap_gsea <- calculate.overlap(x = gsea_list) #get the overlap and unique values of venn diagram
 RUV_unique<-overlap_gsea$a4 #extract RUV-4 unique gsea
 
 RUV_unique<- as.data.frame(RUV_unique)
 write.csv(RUV_unique,"/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/RUV_unique_gsea.csv")
 
venn_plot_GSEA <- ggVennDiagram(gsea_list[1:5], label_alpha = 0 , label = "count", 
                           color =  c("ComBat(76)" = "orange","Limma(939)" ="steelblue",'RAW(797)' = 'black', 'RUV(1036)' = 'red', "SVA(985)" = "purple"),
                           set_color =  c("ComBat(76)" = "orange","Limma(939)" ="steelblue",'RAW(797)' = 'black', 'RUV(1036)' = 'red', "SVA(985)" = "purple")) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  scale_color_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "right")



print(venn_plot_GSEA)

ggsave("/bmbl_data/qujia/GeoMX/pipeline/lymphnode/DE.results/enrichment/GSEA/Venn_GSEA.png", 
       plot = venn_plot_GSEA, width = 14, height = 8, dpi = 300)






