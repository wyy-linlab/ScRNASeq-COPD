
### burden test
## count of deg 
mast_deg_count <- as.data.frame(sapply(mast_deg, function(x){sapply(x, length)}))
mast_deg_count$id <- row.names(mast_deg_count)
mast_deg_count <- mast_deg_count[grep("union",row.names(mast_deg_count)),] %>% reshape2::melt()
mast_deg_count$cell_count <- log10(rep(as.numeric(table(data_integrated$cell_type)),each=3))
mast_deg_count$id[mast_deg_count$id=="COPD_union"] ="COPD"
mast_deg_count$id[mast_deg_count$id=="age_union"] ="age"
mast_deg_count$id[mast_deg_count$id=="smoking_union"] ="smoking"
colnames(mast_deg_count)[c(1,2)]=c("condition","cluster")

mast_deg_count %>% filter(cluster!= "proliferating cells") %>%
  ggplot(aes(x=cell_count,y=log10(value),shape=condition,color=cluster))+
  geom_point(size=5,stroke = 0.8)+
  scale_color_manual(values = cluster_cols)+ scale_shape_manual(values = c(21,8,24))+
  xlab("log10 n_cells ")+ylab("log10 count of DEGs ")+
  theme_classic()
ggsave("/data/wyy/data_yy/hqq_copd/results/fig2_update/count_DEGs_v3.pdf",device = "pdf",width = 6.5,height = 5.5)

options(mc.cores=32)
function_burden_test <- function(cluster,downsampling_n,times){
  function_burden <- function(cluster,downsampling_n){
    seurat_cluster <- data_integrated_regress_cluster[[cluster]]
    if(ncol(seurat_cluster) >= downsampling_n) cells_downsampling <- sample(colnames(seurat_cluster), downsampling_n)   # samplesize >=
    dwonsmple
    if(ncol(seurat_cluster) < downsampling_n)  cells_downsampling <- sample(colnames(seurat_cluster), downsampling_n, replace = T)   # samplesize < dwonsmple , need replace
    
    count_downsampling <- seurat_cluster@assays$RNA@counts[,cells_downsampling]   #count
    count_downsampling_expressed <- count_downsampling[((Matrix::rowSums(count_downsampling>0)/ncol(count_downsampling))>0.1),]  # exp pct >0.1
    meta_downsampling <- seurat_cluster@meta.data[cells_downsampling,]   # metadata
    meta_downsampling$nCount_RNA <- scale(meta_downsampling$nCount_RNA)
    meta_downsampling$age <- scale(meta_downsampling$age)
    colnames(count_downsampling_expressed) <- row.names(meta_downsampling)
    
    dge <- DGEList(counts=count_downsampling_expressed)      # mast --
    dge <- calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    logCPM1 <- log2(cpms+1)
    sca <- FromMatrix(exprsArray = logCPM1, cData = data.frame(wellKey = rownames(meta_downsampling),meta_downsampling))
    mast_zlm <- zlm(~diagnosis + age + (1|orig.ident) + smoking + percent.mt + nCount_RNA , sca, method='glmer', ebayes=FALSE)
    
    doLRT_object<-c("diagnosisCOPD","age","smokingsmoker")  # doLRT --
    mast_res<-data.table()
    
    for(i in doLRT_object){
      mast_summary <- summary(mast_zlm,doLRT= i)
      mast_summary_dt <- mast_summary$datatable
      fcHurdle <- merge(mast_summary_dt[contrast==i & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                        mast_summary_dt[contrast==i & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
      fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
      colnames(fcHurdle)<-paste0(i ,"_",colnames(fcHurdle))
      mast_res<-cbind(mast_res,fcHurdle)
    }
    mast_table <- as.data.frame(mast_res)
    mast_table <- mast_table[,colnames(mast_table)[-grep("_ci.",colnames(mast_table))]]
    mast_table_diag <- na.omit(mast_table[,grep("diagnosis",colnames(mast_table))])
    mast_table_age <- na.omit(mast_table[,grep("age",colnames(mast_table))])
    mast_table_smoking <- na.omit(mast_table[,grep("smoking",colnames(mast_table))])
    
    downsampling_res <- list()       # output
    downsampling_res$count_pvalue_0.05 <- c(sum(mast_table_diag$`diagnosisCOPD_Pr(>Chisq)` < 0.05),sum(mast_table_age$`age_Pr(>Chisq)` < 0.05),sum(mast_table_smoking$`smokingsmoker_Pr(>Chisq)` < 0.05))
    downsampling_res$count_fdr_0.05 <- c(sum(mast_table_diag$diagnosisCOPD_fdr < 0.05),sum(mast_table_age$age_fdr < 0.05),sum(mast_table_smoking$smokingsmoker_fdr < 0.05))
    names(downsampling_res$count_pvalue_0.05) <- c("diagnosis","age","smoking")
    names(downsampling_res$count_fdr_0.05) <- c("diagnosis","age","smoking")
    
    downsampling_res$mast_res_all <- mast_table    # mast table -
    downsampling_res$mast_res_diagnosis <- mast_table_diag
    downsampling_res$mast_res_age <- mast_table_age
    downsampling_res$mast_res_smoking <- mast_table_smoking
    
    downsampling_res$count <- count_downsampling
    downsampling_res$count_expressed <- count_downsampling_expressed
    downsampling_res$meta <- meta_downsampling
    
    return(downsampling_res)   
  }
  
  downsampling_times_res <- list()
  for(time in 1:times){
    downsampling_times_res[[time]] <- function_burden(cluster,downsampling_n)
    print(paste(cluster,"time:",time,"----- done"))
  }
  names(downsampling_times_res) <- paste0(cluster,"_downsampling_",seq(1,times,by = 1))
  return(downsampling_times_res)
}

load("COPD_downsampling.RData")

function_burden_boxplot <- function(downsampling_res,ncells,pq,object){
  downsampling <- sapply(cluster_names[1:14],function(cluster){sapply(downsampling_res[[cluster]],
                                                                      function(data){data[[pq]][object]})})
  downsampling <- as.data.frame(downsampling)
  downsampling$times <- paste0("downsampling_",seq(1:10))
  downsampling %>% reshape2::melt(id="times") %>%
    ggplot(aes(x=variable,y=value,col=variable))+
    ggtitle(paste0("burden analysis : ",object))+
    geom_boxplot()+
    scale_color_manual(values = c("#F8766D","#E58700","#00BF7D","#00B0F6","#E76BF3",
                                  "#00BA38","#FD61D1","#00C0AF","#619CFF","#FF67A4",
                                  "#C99800","#00BCD8","#B983FF","#A3A500","#6BB100"))+
    ylab("count of significant genes")+
    theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1),
                          axis.title.x = element_blank(),
                          legend.position = "none")
}

function_burden_boxplot(downsampling_res = burden_result$n700,ncells = 700,pq = "count_fdr_0.05",object = "diagnosis")
ggsave("/data/js/project/COPD/result/burden/burden_plot_diagnosis.pdf",device = "pdf",width = 3.2,height = 3)
function_burden_boxplot(downsampling_res = burden_result$n700,ncells = 700,pq = "count_fdr_0.05",object = "age")
ggsave("/data/js/project/COPD/result/burden/burden_plot_age.pdf",device = "pdf",width = 3.2,height = 3)
function_burden_boxplot(downsampling_res = burden_result$n700,ncells = 700,pq = "count_fdr_0.05",object = "smoking")
ggsave("/data/js/project/COPD/result/burden/burden_plot_smoking.pdf",device = "pdf",width = 3.2,height = 3)


