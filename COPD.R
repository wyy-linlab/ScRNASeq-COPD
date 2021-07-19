# COPD code 

library(Matrix)
library(Seurat) 
library(stringr)
library(uwot) 
library(patchwork)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ggpubr)
library("MAST")  # 1.14.0
library(lme4)
library(edgeR)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggpubr)
library(grid)
library(cowplot)
library(apcluster) 
library(clustree) 
library(SingleR)
library(scHCL)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle)
library(DDRTree)
library(plyr)
library(dplyr)
library(flashClust)
library(slingshot)
library(mclust)
library("scales")
library(slingshot)
library(SingleCellExperiment)
library(destiny)



### data preparation --- 
data_dgcMatrix <- list(ctl.o1 = as(as.matrix(read.table("dgcMatrix/ctl.o1.txt",sep="\t")),"dgCMatrix"),    
                       ctl.o2 = as(as.matrix(read.table("dgcMatrix/ctl.o2.txt",sep="\t")),"dgCMatrix"),
                       copd.o1 = as(as.matrix(read.table("dgcMatrix/copd.o1.txt",sep="\t")),"dgCMatrix"), 
                       ctl.y1 = as(as.matrix(read.table("dgcMatrix/ctl.y1.txt",sep="\t")),"dgCMatrix"),
                       ctl.o3 = as(as.matrix(read.table("dgcMatrix/ctl.o3.txt",sep="\t")),"dgCMatrix"),   
                       copd.o2 = as(as.matrix(read.table("dgcMatrix/copd.o2.txt",sep="\t")),"dgCMatrix"),
                       copd.o3 = as(as.matrix(read.table("dgcMatrix/copd.o3.txt",sep="\t")),"dgCMatrix"), 
                       ctl.y2 = as(as.matrix(read.table("dgcMatrix/ctl.y2.txt",sep="\t")),"dgCMatrix"),
                       ctl.y3 = as(as.matrix(read.table("dgcMatrix/ctl.y3.txt",sep="\t")),"dgCMatrix"))  

function_CreateSeuratObject <- function(data){   
  data <- CreateSeuratObject(data,assay = "RNA")
  data$percent.mt <- PercentageFeatureSet(data, pattern = "^MT-")/100
  data <- data %>% subset(subset = nFeature_RNA > 300 & percent.mt < 0.3) %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000,verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000,verbose = F)
  return(data)
}

data_seurat <- lapply(data_dgcMatrix, function_CreateSeuratObject)  
mean(unlist(sapply(data_seurat, function(data){data$nFeature_RNA})))*1.5  

data_seurat.raw <- lapply(data_dgcMatrix, function_CreateSeuratObject_raw)

## merge: convert into a format that is convenient for subsequent analysiC

merged_seurat <- merge(x = data_seurat.raw$ctl.o1, 
                       y = data_seurat.raw[2:9], 
                       add.cell.id = names(data_seurat.raw))


merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

data_seurat <- subset(x=merged_seurat,
                          subset = nFeature_RNA > 300 & percent.mt < 0.3)



## integrate---
data_seurat.anchors <- FindIntegrationAnchors(data_seurat,dims = 1:30,anchor.features = 3000,verbose = F)
data_integrated <- IntegrateData(data_seurat.anchors,dims = 1:30,verbose = F)
data_integrated <- CellCycleScoring(data_integrated,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)


## add metedata---
data_integrated$individual <- as.character(data_integrated$orig.ident)
data_integrated$individual[data_integrated$individual=="copd.o1"] <- "COPD1"  
data_integrated$individual[data_integrated$individual=="copd.o2"] <- "COPD2"
data_integrated$individual[data_integrated$individual=="copd.o3"] <- "COPD3"

data_integrated$diagnosis <- factor(word(data_integrated$group,1,sep=fixed(".")),levels = c("ctl","COPD"))
data_integrated$age <-rep(c(73,71,63,28,75,73,50,35,27), times = sapply(data_seurat,ncol))
data_integrated$smoking = rep(c("smoker","smoker","smoker","non-smoker","non-smoker","non-smoker","smoker","non-smoker","non-smoker"),sapply(data_seurat,ncol))

data_integrated$group <- str_sub(data_integrated$orig.ident,1,-2)
data_integrated$group[data_integrated$group=="copd.o"] = "COPD"

data_integrated$cell_cycle <- as.character(data_integrated$Phase)  # cell cycle phase 
data_integrated$cell_cycle[data_integrated$cell_cycle=="G2M"]<-"dividing"
data_integrated$cell_cycle[data_integrated$cell_cycle!="dividing"]<-"non-dividing"


### scale data and FindCluster----
data_integrated <- ScaleData(data_integrated,verbose = F)
data_integrated_regress <- ScaleData(data_integrated,vars.to.regress = c("percent.mt","nCount_RNA"),verbose = F)

DefaultAssay(data_integrated) <- "integrated"

data_integrated_regress <- data_integrated_regress %>% RunPCA(npcs = 50,verbose = F) %>% 
  RunTSNE(reduction = "pca",dims = 1:14,verbose = F) %>% RunUMAP(reduction = "pca",dims = 1:14,verbose = F)  %>% 
  FindNeighbors(reduction = "pca", dims = 1:14,verbose = F) %>% FindClusters(resolution = c(0.1,0.2,0.3),verbose = F)  


## marker plot and define cell type---
cluster_cols=c("#F8766D","#E58700","#00BF7D","#00B0F6","#E76BF3",
               "#00BA38","#FD61D1","#00C0AF","#619CFF","#FF67A4",
               "#C99800","#00BCD8","#B983FF","#A3A500","#6BB100")
cluster_shapes <- c(rep(19,8),rep(17,7))
cluster_names <- c("macrophages","dendritic cells","monocytes","mast cells","neutrophils","natural killer cells","T cells","B cells","alveolar type 1 cells","alveolar type 2 cells","club cells","ciliated cells","stromal cells","endothelial cells","proliferating cells")


function_FeaturePlot <- function(data,features,reduction,ncol){
  DefaultAssay(data) <- "RNA"
  p <- FeaturePlot(data,features = features,reduction=reduction,cols = c("lightgrey","blue"),pt.size = 0.1,combine = F)
  p_add_theme <- lapply(p, function(x){
    return(x+theme_classic()+theme(legend.position = "none",
                                   axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
                                   plot.title = element_text(hjust = 0.5)))
  })
  wrap_plots(p_add_theme,ncol = ncol)
}
cluster_markers <- c("MARCO","FABP4","CLEC10A","CD1C","CD300E","FCN1","MS4A2","CPA3","FCGR3B","PROK2","KLRD1","NKG7","TRAC","CD79A","IGKC","AGER","PDPN","SFTPC","LAMP3","SCGB1A1","SCGB3A1","FOXJ1","PIFO","ACTA2","COL1A1","PDGFRA","CLDN5","VWF","MKI67","TOP2A")

function_FeaturePlot(data_integrated,cluster_markers[1:15],reduction="umap",ncol =5) 
function_FeaturePlot(data_integrated,cluster_markers[16:30],reduction="umap",ncol = 5) 

Idents(data_integrated_regress) <- data_integrated_regress$integrated_snn_res.0.3
data_integrated_regress <- RenameIdents(data_integrated_regress,"0"="macrophages_1","1"="macrophages_2","2"="alveolar type 2 cells","3"="T cells","4"="natural killer cells","5"="dendritic cells","6"="endothelial cells","7"="ciliated cells","8"="stromal cells","9"="mast cells","10"="club cells","11"="monocytes","12"="alveolar type 1 cells","13"="B cells_1","14"="neutrophils","15"="proliferating cells","16"="B cells_2")
data_integrated_regress$cell_type<-as.character(word(as.character(Idents(data_integrated_regress)),1,sep=fixed("_")))
data_integrated_regress$cell_type<-factor(data_integrated_regress$cell_type, levels=c("macrophages","dendritic cells","monocytes","mast cells","neutrophils","natural killer cells","T cells","B cells","alveolar type 1 cells","alveolar type 2 cells","club cells","ciliated cells","stromal cells","endothelial cells","proliferating cells"))

wrap_plots(list(DimPlot(data_integrated,group.by = "cell_type",pt.size = 0.1,label = T,cols = cluster_cols),
                DimPlot(data_integrated_regress,group.by = "cell_type",label = T,cols = cluster_cols)))
cluster_names <- levels(data_integrated$cell_type)


## Meta plot and percent ---
function_MetaPlot<-function(object,color,model){     # color for categrorical variable
  plot_data <- data.frame(umap_1=Embeddings(data_integrated,reduction = "umap")[,1],
                          umap_2=Embeddings(data_integrated,reduction = "umap")[,2],
                          plot_object=data_integrated[[object]][,1])
  if(model=="n")  # numeric
    p<-ggplot(plot_data,aes(umap_1,umap_2,color=plot_object))+
      geom_point(size=0.1)+ scale_color_gradient(low="grey",high="blue")+
      ggtitle(object) +
      theme_classic() + theme(legend.title = element_blank())
  else    # categorical
    p<-ggplot(plot_data,aes(umap_1,umap_2,color=plot_object)) +
      geom_point(size=0.1,alpha=0.1) + scale_color_manual(values =color) +
      ggtitle(object) +
      theme_classic()+theme(legend.title = element_blank(),legend.position = "none") +
      guides(colour = guide_legend(override.aes = list(size=3,alpha=0.9)))
   return(p)
}
function_MetaPlot("nCount_RNA", color = NULL, model = "n")
function_MetaPlot("nFeature_RNA", color = NULL, model = "n")
function_MetaPlot("percent.mt", color = NULL, model = "n")
function_MetaPlot("group", color = c("#FF6666","#87CEFF","#B2DF8A"), model = "c")   #
function_MetaPlot("smoking", color = c("#0000FF","#FFA07A"), model = "c") #
function_MetaPlot("cell_cycle", color = c("#FF6666","#87CEFF"), model = "c") #
function_MetaPlot("diagnosis", color =c("#B2DF8A","#1E90FF"), model = "c")
function_MetaPlot("individual", color =c("#F8766D","#E58700","#00BF7D","#00B0F6","#E76BF3","#00BA38","#FD61D1","#00C0AF","#619CFF","#FF6
7A4","#C99800","#00BCD8","#B983FF","#A3A500","#6BB100"),model = "c")


## metadata percent 
function_Meta_percent <- function(plot_object,color){
  data <- data_integrated
  data$cell_type <- factor(data$cell_type,levels = rev(levels(data$cell_type)))  
  p <- ggplot(data@meta.data,aes(x=cell_type,fill=data@meta.data[[plot_object]]))+
    geom_bar(position = "fill",alpha=0.7) + scale_fill_manual(values=color)+
    coord_flip()+
    ggtitle(plot_object) + xlab("cluster") + ylab("fraction of cells") +
    theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y = element_blank(),legend.title = element_blank())
    return(p)
}
function_Meta_percent(plot_object = "smoking",color = c("#0000FF","#FFA07A"))
function_Meta_percent(plot_object = "group",color = c("#FF6666","#B2DF8A","#87CEFF"))
function_Meta_percent(plot_object = "cell_cycle",color = c("#FF6666","#87CEFF"))
function_Meta_percent(plot_object = "diagnosis",color = c("#B2DF8A","#1E90FF"))
function_Meta_percent(plot_object = "individual",color = c("#F8766D","#E58700","#00BF7D","#00B0F6","#E76BF3","#00BA38","#FD61D1","#00C0AF","#619CFF","#FF67A4","#C99800","#00BCD8","#B983FF","#A3A500","#6BB100"))

## subClustering monocytes 
data_integrated_cluster <- SplitObject(data_integrated,split.by = "cell_type")
subcluster_cluster <- list()

subcluster_cluster$monocytes <- data_integrated_cluster$monocytes
DefaultAssay(subcluster_cluster$monocytes) <- "integrated"
subcluster_cluster$monocytes$umap1 <- subcluster_cluster$monocytes@reductions$umap@cell.embeddings[,1]
subcluster_cluster$monocytes$umap2 <- subcluster_cluster$monocytes@reductions$umap@cell.embeddings[,2]
subcluster_cluster$monocytes <- subset(subcluster_cluster$monocytes, subset=umap1< 2 & umap2 < 3.5)
wrap_plots(list(DimPlot(data_integrated_cluster$monocytes),
                DimPlot(subcluster_cluster$monocytes)))
subcluster_cluster$monocytes <- subcluster_cluster$monocytes %>%     # t-SNE
  RunTSNE(dims = 1:14,verbose = F) %>%
  FindNeighbors(reduction = "pca", dims = 1:14,verbose = F) %>% FindClusters(resolution = seq(0.1,0.3,by = 0.1),verbose = F)

DefaultAssay(subcluster_cluster$monocytes) <- "RNA"

subcluster_cluster$monocytes$subcl <- as.character(subcluster_cluster$monocytes$integrated_snn_res.0.2)    # subCluster type -
subcluster_cluster$monocytes$subcl[subcluster_cluster$monocytes$subcl %in% c(0,2)]="CD14+mono"
subcluster_cluster$monocytes$subcl[subcluster_cluster$monocytes$subcl == 1]="CD16+mono"
subcluster_cluster$monocytes$subcl <- factor(subcluster_cluster$monocytes$subcl,levels = c("CD14+mono","CD16+mono"))

DimPlot(subcluster_cluster$monocytes,reduction = "tsne",group.by = "subcl",cols = c("#F8766D","#00B0F6"),label = T)+xlab("t-SNE1")+ylab("t-SNE2")
ggsave("../result/Fig3/Fig3_monocytes_subcluster.pdf",device = "pdf",width = 4.9,height = 3.5)


### MAST ~  differential expression genes---
options(mc.cores=32)
function_MAST <- function(data,cluster,model){
  scData <- data[[cluster]]
  expressed <- rowSums(scData@assays$RNA@counts>0)/ncol(scData) > 0.1 
  dge <- DGEList(counts=scData@assays$RNA@counts[expressed,])
  dge <- calcNormFactors(dge)
  cpms <- edgeR::cpm(dge)
  logCPM1 <- log2(cpms+1)
  test_meta <- scData@meta.data   
  test_meta$nCount_RNA <- scale(test_meta$nCount_RNA)
  test_meta$age <- scale(test_meta$age)
  sca <- FromMatrix(exprsArray = logCPM1, cData = data.frame(wellKey = rownames(test_meta),test_meta))
  
  if(model == "random")   mast_zlm <- zlm( ~ diagnosis + age + (1|individual) + smoking + percent.mt + nCount_RNA , sca, method='glmer', ebayes=FALSE)
  if(model == "nr") mast_zlm <- zlm(~diagnosis + age + smoking + percent.mt + nCount_RNA , sca, ebayes=T)
  
  mast_res <- data.table()
  doLRT_object<-c("diagnosisCOPD","age","smokingsmoker")
  for(i in doLRT_object){
    mast_summary <- summary(mast_zlm , doLRT= i) 
    mast_summary_dt <- mast_summary$datatable
    fcHurdle <- merge(mast_summary_dt[contrast==i & component=='H',.(primerid, `Pr(>Chisq)`)], 
                      mast_summary_dt[contrast==i & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    colnames(fcHurdle) <- paste0(i ,"_",colnames(fcHurdle))
    mast_res <- cbind(mast_res,fcHurdle)
  }
  mast_table <- as.data.frame(mast_res)
  mast_table <- mast_table[,colnames(mast_table)[-grep("_ci.",colnames(mast_table))]] 
  if(model == "random"){  fdr = 0.05
  coef = 0}
  if(model=="nr"){  fdr =0.05
  coef =1}
  mast_table_diag <- na.omit(mast_table[,grep("diagnosis",colnames(mast_table))])
  mast_table_age <- na.omit(mast_table[,grep("age",colnames(mast_table))])
  mast_table_smoking <- na.omit(mast_table[,grep("smoking",colnames(mast_table))])
  
  deg <- list(COPD_up =    mast_table_diag[mast_table_diag$diagnosisCOPD_fdr < fdr & mast_table_diag$diagnosisCOPD_coef > coef, 1],
              COPD_down =  mast_table_diag[mast_table_diag$diagnosisCOPD_fdr < fdr & mast_table_diag$diagnosisCOPD_coef < -coef, 1],
              COPD_union = mast_table_diag[mast_table_diag$diagnosisCOPD_fdr < fdr & abs(mast_table_diag$diagnosisCOPD_coef) > coef, 1],
              age_up =    mast_table_age[mast_table_age$age_fdr < fdr & mast_table_age$age_coef > coef, 1],
              age_down =  mast_table_age[mast_table_age$age_fdr < fdr & mast_table_age$age_coef < -coef, 1],
              age_union = mast_table_age[mast_table_age$age_fdr < fdr & abs(mast_table_age$age_coef) > coef , 1],
              smoking_up =    mast_table_smoking[mast_table_smoking$smokingsmoker_fdr < fdr & mast_table_smoking$smokingsmoker_coef > coef, 1],
              smoking_down =  mast_table_smoking[mast_table_smoking$smokingsmoker_fdr < fdr & mast_table_smoking$smokingsmoker_coef < -coef, 1],
              smoking_union = mast_table_smoking[mast_table_smoking$smokingsmoker_fdr < fdr & abs(mast_table_smoking$smokingsmoker_coef) > coef, 1])
  print(paste0(cluster," MAST done !"))
  return(list(table=mast_table,
              deg=deg))
}
data_integrated_regress_cluster <- SplitObject(data_integrated_regress , "cell_type")
mast_result <- lapply(cluster_names, function_MAST, data = data_integrated_regress_cluster, model="random")
names(mast_result) <- cluster_names


## exact deg --
mast_deg <- lapply(mast_result, function(data){data$deg})
sapply(mast_deg$macrophages, length)
## MAST table --
mast_table <- do.call("rbind", lapply(mast_result,function(data){data$table}))
mast_table$cluster = factor(rep(names(mast_result),times=sapply(mast_result, function(data){nrow(data$table)})),levels = cluster_names)
mast_table_diagnosis <- na.omit(mast_table[,c(colnames(mast_table)[grep("diagnosis",colnames(mast_table))],"cluster")])
mast_table_age <- na.omit(mast_table[,c(colnames(mast_table)[grep("age",colnames(mast_table))],"cluster")])
mast_table_smoking <- na.omit(mast_table[,c(colnames(mast_table)[grep("smoking",colnames(mast_table))],"cluster")])

## mast table output --
sapply(c("diagnosis","age","smoking"), function(name){
  mast_table_sub <- na.omit(mast_table[,c(colnames(mast_table)[grep(name,colnames(mast_table))],"cluster")])
  colnames(mast_table_sub)[1:4] <- c("gene","pvalue","log2fc","fdr")
  write.table(mast_table_sub, file = paste0("../result/MAST/MAST_table_",name,".txt"),sep="\t",quote = F,row.names = F)
})


### slingshot ~ trajectory analysis  
seurat_integrated <- data_integrated
{DefaultAssay(seurat_integrated) <- "integrated"
  seurat_integrated <- seurat_integrated %>% RunPCA()
  
  seurat_sub <- subset(seurat_integrated_add_new_metadata, subset = cell.type.js %in% c("alveolar-type-1-cells","alveolar-type-2-cells","ciliated-cells","club-cells") )
  
  
  
  seurat_slingshot <- SingleCellExperiment(assays = List(counts = seurat_sub@assays$SCT@counts[seurat_sub@assays$integrated@var.features,]),
                                           colData=seurat_sub@meta.data)
  assays(seurat_slingshot)$norm <- seurat_sub@assays$integrated@data
  
  
  
  pca <- prcomp(t(assays(seurat_slingshot)$norm), scale. = T)
  rd1 <- pca$x[,1:2]
  
  plot(rd1, col = as.factor(seurat_sub$lv2), pch=16, asp = 1)
  dm <- DiffusionMap(t(as.matrix(seurat_integrated@assays$integrated@data)),n_eigs=5,k=10)
  
  #dm <- DiffusionMap((as.matrix(seurat_integrated@reductions$pca@cell.embeddings)))
  rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
  plot(rd2, col =as.factor(seurat_integrated$lv2) , pch=16, asp = 1)
  
  reducedDims(seurat_slingshot) <- SimpleList(pca = rd1,diff = rd2)
  
  colData(seurat_slingshot)$GMM <- seurat_integrated$integrated_snn_res.0.1
  cl2 <- kmeans(rd1, centers = 4)$cluster
  colData(seurat_slingshot)$kmeans <- cl2
  sim <- slingshot(seurat_slingshot, clusterLabels=seurat_integrated$integrated_snn_res.0.1,reducedDim = 'diff',start.clus=3,end.clus=0) #clusterLabels = 'GMM',
  sim_pca <- slingshot(seurat_slingshot, clusterLabels=seurat_integrated$integrated_snn_res.0.1,reducedDim = 'pca',start.clus=3,end.clus=0) #clusterLabels = 'GMM',
  
  summary(sim$slingPseudotime_1)
  plot(reducedDims(sim)$diff, col = seurat_integrated$integrated_snn_res.0.1, pch=16, asp = 1)
  lines(SlingshotDataSet(sim), lwd=2, col='black')
  lines(curves, lwd = 1, col = 'black')
  
  all_sim <- cbind(reducedDims(sim)$pca,seurat_integrated@meta.data)
  curve <- slingCurves(sim)$curve1$s
  
### save data---
save.image("COPD.RData")



