
### WGCNA----
library(WGCNA)
function_pickSoftThreshold <- function(data){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],type="n",xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  return(sft$powerEstimate)
}


### at2 gene expression correction
at2 <- data_integrated_cluster$`alveolar type 2 cells`
at2_expressed <- rowSums(as.matrix(at2@assays$RNA@counts>0))/ncol(at2) > 0.1 
at2_exp <- at2@assays$RNA@counts[at2_expressed,] %>% DGEList() %>% calcNormFactors() %>% edgeR::cpm()
at2_log_exp <- log2(at2_exp+1)
at2_variable <- at2@meta.data[,c(1,2,4,9,11,12)]
at2_corr <- empiricalBayesLM(t(at2_log_exp),at2_variable[,c(2,3,5,6)],at2_variable[,4]) 
at2_corr <-t(at2_corr$adjustedData)
at2_variable_filter <- filter(at2_variable,nCount_RNA>5000)  # 
at2_corr_filter <- at2_corr[,at2_variable$nCount_RNA>5000]  #
at2_power <- function_pickSoftThreshold(at2_corr_filter)

### at2 module establishment
TOM=TOMsimilarityFromExpr(t(at2_corr_filter),networkType = "signed", TOMType = "signed", power = 16,nThreads=40,corType="bicor")
colnames(TOM) =rownames(TOM) =rownames(at2_corr)
dissTOM=1-TOM
geneTree = flashClust(as.dist(dissTOM),method="average")
plot(geneTree, xlab="", sub="",cex=0.3)
minModuleSize = 50
function_cutreeDynamic <- function(cutHeight){
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,cutHeight=cutHeight,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
  return(labels2colors(dynamicMods))
}
dynamicColors_all <- data.frame(c95=function_cutreeDynamic(0.95),c96=function_cutreeDynamic(0.96),c97=function_cutreeDynamic(0.97),c98=function_cutreeDynamic(0.98),c99=function_cutreeDynamic(0.99))
plotDendroAndColors(geneTree, dynamicColors_all,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dynamicColors <- function_cutreeDynamic(0.96)
names(dynamicColors) <- rownames(at2_corr_filter)
dynamicColors_table <- data.frame(color = as.character(dynamicColors),
                                  gene=names(dynamicColors))
dynamicColors_table<-dynamicColors_table[order(dynamicColors_table$color),]
write.table(dynamicColors_table,"../result/WGCNA/module_gene.txt",sep = "\t",quote = F,row.names = F)

### at2 module expression, trait association
MEs0 = moduleEigengenes(t(at2_corr_filter), dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, as.numeric(at2_variable_filter$diagnosis)-1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 3391)
moduleTraitQvalue = p.adjust(moduleTraitPvalue,method = "fdr")

design <- model.matrix(~ diagnosis,at2_variable_filter)
fit <- lmFit(t(MEs), design)
fit <- eBayes(fit)
head(design)
MEs_pvalue <- topTable(fit, coef="diagnosisCOPD", number=Inf,sort.by="none")
module_plot <- cbind.data.frame(p=-log10(moduleTraitQvalue),size=c(780,838,1168,151,250,1645),fc=MEs_pvalue$logFC)

ggplot(module_plot,aes(x=fc,y=p))+
  geom_point(aes(size=size),colour=c("brown","blue","turquoise","green","yellow","grey"))+
  geom_vline(xintercept=0,lty="dashed")+
  xlab("module_expression_log_fold_change")+ylab("-log10(fdr)")+
  theme_classic()

blue_gene <- rownames(at2_corr_filter)[dynamicColors=="blue"]
write(blue_gene,"../result/AT2/blue_gene.txt")    

### hub gene 
dissTOM_blue <- dissTOM[dynamicColors=="blue",dynamicColors=="blue"]
blue_tree <- flashClust(as.dist(dissTOM_blue), method = "centroid")
plotTOM = (1-dissTOM_blue)^6
diag(plotTOM) = NA
TOMplot(plotTOM, blue_tree, dynamicColors[dynamicColors=="blue"], main = "")  # Fig 5D 
hub1 <- blue_tree$labels[blue_tree$order] # hub genes
geneModuleMembership = as.data.frame(cor(t(at2_corr_filter), MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 3391));
hub2 <- rownames(geneModuleMembership)[order(geneModuleMembership$MEblue,decreasing = T)] #hub genes

### heatmap 
blue_expr <- at2_corr_filter[dynamicColors=="blue",]
blue_expr_copd <- blue_expr[,at2_variable_filter$diagnosis=="COPD"]
blue_expr_ctrl <- blue_expr[,at2_variable_filter$diagnosis=="ctl"]
n=30

blue_expr_copd_mean <- collapseRows(t(blue_expr_copd),rowGroup = copd_group,rowID = rownames(t(blue_expr_copd)),method="Average")$datETcollapsed %>% t()
blue_expr_ctrl_mean <- collapseRows(t(blue_expr_ctrl),rowGroup = ctrl_group,rowID = rownames(t(blue_expr_ctrl)),method="Average")$datETcollapsed %>% t()
blue_expr_mean <- cbind(blue_expr_copd_mean,blue_expr_ctrl_mean)
anno=data.frame(group=c(rep("copd",n),rep("ctrl",n)   ))
rownames(anno)=colnames(blue_expr_mean)
pheatmap(blue_expr_mean,scale="row",show_rownames = F,show_colnames = F,cluster_rows=T,colorRampPalette(c("blue","white","red"))(16),
         annotation_col = anno,clustering_method = "average",clustering_distance_cols = "correlation")   # Fig 5B

load("COPD_SCT.Rdata")  
identical(as.numeric(seurat_integrated_wyy$nCount_RNA), as.numeric(data_integrated_regress$nCount_RNA))
seurat_integrated_wyy$cell_type <- data_integrated_regress$cell_type
seurat_integrated_wyy$smoking <- data_integrated_regress$smoking
seurat_integrated_wyy$group <- data_integrated_regress$group
seurat_integrated_wyy$individual <- data_integrated_regress$individual
seurat_integrated_wyy$diagnosis <- data_integrated_regress$diagnosis
seurat_integrated_wyy_cluster <- SplitObject(seurat_integrated_wyy,split.by = "cell_type")

AT2_RNA_expressed <- rowSums(data_integrated_regress_cluster$`alveolar type 2 cells`@assays$RNA@counts>0)/ncol(data_integrated_regress_cluster$`alveolar type 2 cells`) > 0.1 
AT2_RNA_logCPM <- log2(edgeR::cpm(calcNormFactors(DGEList(counts=data_integrated_regress_cluster$`alveolar type 2 cells`@assays$RNA@counts[AT2_RNA_expressed,])))+1)
AT2_RNA_CPM <- edgeR::cpm(calcNormFactors(DGEList(counts=data_integrated_regress_cluster$`alveolar type 2 cells`@assays$RNA@counts[AT2_RNA_expressed,])))

function_correct_vlnplot <- function(gene,assay,model,method,y1,y2){
  if(assay=="RNA"){
    if(model=="log"){  exp_data <- AT2_RNA_logCPM[gene,]}
    if(model=="unlog"){  exp_data <- AT2_RNA_CPM[gene,]}
    exp_data
  }
  if(assay=="SCT"){
    if(model=="unlog"){  exp_data <- exp(seurat_integrated_wyy_cluster$`alveolar type 2 cells`@assays$SCT@data[gene,])}
    if(model=="log"){ exp_data <- seurat_integrated_wyy_cluster$`alveolar type 2 cells`@assays$SCT@data[gene,]}
    exp_data
  }
  AT2_smoking <- (as.numeric(factor(seurat_integrated_wyy_cluster$`alveolar type 2 cells`$smoking))-1)
  AT2_age <- scale(seurat_integrated_wyy_cluster$`alveolar type 2 cells`$age.js)
  gene_coef = mast_result$`alveolar type 2 cells`$table[mast_result$`alveolar type 2 cells`$table$diagnosisCOPD_primerid==gene,]
  
  plot_data <- data.frame(row.names = names(exp_data),
                          raw_exp = exp_data,
                          exp=exp_data - AT2_smoking*(gene_coef$smokingsmoker_coef)- AT2_age*(gene_coef$age_coef),
                          diagnosis = seurat_integrated_wyy_cluster$`alveolar type 2 cells`$diagnosis,
                          group=seurat_integrated_wyy_cluster$`alveolar type 2 cells`$group)
  if(method == "filter"){  plot_data <- plot_data %>% filter(raw_exp > 0)  }
  
  p<-plot_data %>% ggplot(aes(x=diagnosis,y=exp,fill=diagnosis))+
    geom_violin()+
    scale_fill_manual(values = c("#6495ED","#EE3B3B"))+
    ggtitle(gene)+coord_cartesian(ylim = c(y1, y2))+
    theme_classic()+theme(axis.title.x = element_blank())
  print(p)
}








