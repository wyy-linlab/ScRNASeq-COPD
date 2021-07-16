
### enrichment analysis
## GWAS risk enrichment
gwas_risk <- as.character(read.table("gwas/GWAS.tsv",sep="\t",header = T)$REPORTED.GENE.S.)
sapply(cluster_names, function(cluster){intersect(mast_deg[[cluster]]$COPD_union,gwas_risk)})
sapply(cluster_names, function(cluster){intersect(mast_deg[[cluster]]$COPD_up,gwas_risk)})

function_fisher_GWAS_diagnosis <- function(cluster,GWAS_risk){   # fisher exact test --
  deg <- mast_deg[[cluster]]$COPD_union
  background <- mast_result[[cluster]]$table$diagnosisCOPD_primerid
  risk <- intersect(GWAS_risk,background)  # gwas_risk must in background (gene expressed > 0.1)
  its <- length(intersect(deg,risk))
  matrix_test <- matrix(c(its , length(deg) - its , length(risk) - its , length(background) -length(deg) - (length(risk) - its)),ncol=2)
  fisher_result <- fisher.test(matrix_test,alternative = "greater")
  return(fisher_result)
}
fisher_GWAS_diagnosis <-lapply(cluster_names[1:14], function_fisher_GWAS_diagnosis, GWAS_risk = gwas_risk)
names(fisher_GWAS_diagnosis) <- cluster_names[1:14]
t(sapply(fisher_GWAS_diagnosis, function(data){c(p.value=data$p.value,data$estimate)}))


# fisher plot --
fisher_GWAS_diagnosis_plot <- data.frame(cluster = names(fisher_GWAS_diagnosis),
                                         gwas = rep("gwas",length(names(fisher_GWAS_diagnosis))),
                                         OR = sapply(fisher_GWAS_diagnosis, function(data){data$estimate})) %>% filter(OR >0)
fisher_GWAS_diagnosis_plot$cluster <- factor(fisher_GWAS_diagnosis_plot$cluster,levels=rev(fisher_GWAS_diagnosis_plot[order(fisher_GWAS_diagnosis_plot$OR,decreasing = T),"cluster"]))

ggplot(fisher_GWAS_diagnosis_plot,aes(x=cluster,y=gwas)) +
  geom_tile(alpha=0.5,aes(fill= OR)) + scale_fill_gradient(low="#FFF68F",high="#E31A1C")+
  geom_text(aes(label= round(OR,2)),size=3)+
  coord_flip()+
  theme_classic()+theme(axis.text.x = element_text(hjust = 1),axis.title = element_blank())
ggsave("../result/MAST/gwas_OR.pdf",device = "pdf",height = 4,width = 3)


## Go enrichment
function_clusterprolifer_result <- function(list,name,width,height){
  compare=compareCluster(geneCluster = list, fun = "enrichGO",keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH")
  compare_sim <- clusterProfiler::simplify(compare, cutoff= 0.7, by="p.adjust", select_fun=min)
  compare_sim@compareClusterResult <- na.omit(compare_sim@compareClusterResult)
  dotplot(compare_sim,font.size=8,showCategory=20)
  ggsave(paste0("../result/MAST/clusterprolifer_",name,".pdf"),device = "pdf",width = width,height = height)
  return(compare_sim)
}
clusterProfiler_result <- list(copd_union = function_clusterprolifer_result(list = sapply(cluster_names[1:14], function(cluster){mast_deg[[cluster]]$COPD_union}),name = "diagnosis_union_top20",width = 18,height = 33),
                               age_union = function_clusterprolifer_result(list = sapply(cluster_names[1:14], function(cluster){mast_deg[[cluster]]$age_union}),name = "age_union_top20",width = 18,height = 49),
                               smoking_union = function_clusterprolifer_result(list = sapply(cluster_names[1:14], function(cluster){mast_deg[[cluster]]$smoking_union}),name = "smoking_union_top20",width = 18,height = 33))
names(clusterProfiler_result) <- c("diagnosis","age","smoking")

write.table(clusterProfiler_result$diagnosis@compareClusterResult,"../data/clusterProlifer/clusterProfiler_diagnosis.txt",sep="\t",row.names = F,quote = F)
write.table(clusterProfiler_result$age@compareClusterResult,"../data/clusterProlifer/clusterProfiler_age.txt",sep="\t",row.names = F,quote = F)
write.table(clusterProfiler_result$smoking@compareClusterResult,"../data/clusterProlifer/clusterProfiler_smoking.txt",sep="\t",row.names = F,quote = F)


# monocytes COPD up network and table(Fig3A) ----
function_creat_network_form_ClusterProlifer_or_DAVID <- function(file,collapse){
  network <- matrix(c(rep(file$pathways[1],length(unique(unlist(str_split(str_c(file[file$pathways==file$pathways[1],]$gene,collapse = collapse),pattern = collapse))))),
                      unique(unlist(str_split(str_c(file[file$pathways==file$pathways[1],]$gene,collapse = collapse),pattern = collapse)))),ncol=2)  # pathway 1
  for( i in unique(file$pathways)[-1]){   # pathway 2-end
    gene_list <- unique(unlist(str_split(str_c(file[file$pathways==i,]$gene,collapse = collapse),pattern = collapse)))
    table <- matrix(c(rep(i,length(gene_list)),gene_list),ncol = 2)
    network  <- rbind(network ,table)
  }
  network  <- as.data.frame(network )
  colnames(network ) <- c("pathway","gene")
  return(network)
}

function_creat_table_from_network <- function(network){
  table <- data.frame(row.names = unique(network$gene),
                      gene=unique(network$gene))
  for(i in unique(network$pathway)){
    table[gsub(" ","_",i)] <- as.numeric(table$gene %in% network[network$pathway==i,]$gene)
  }
  table$sum <- rep("piechart",nrow(table))
  return(table)
}

monocytes_clusterProlifer_COPD_up_end <- as.data.frame(read_excel("../result/monocytes_network/monocytes_clusterProlifer_COPD_up_end.xls
x",sheet = 1))
colnames(monocytes_clusterProlifer_COPD_up_end)[c(1,9)] <- c("pathways","gene")
table(monocytes_clusterProlifer_COPD_up_end$pathways)

monocytes_DAVID_COPD_up_end <- as.data.frame(read_excel("../result/monocytes_network/monocytes_DAVID_COPD_up_end.xlsx",sheet = 1))
colnames(monocytes_DAVID_COPD_up_end)[c(1,6)] <- c("pathways","gene")

monocytes_clusterProlifer_COPD_up_end_network <- function_creat_network_form_ClusterProlifer_or_DAVID(monocytes_clusterProlifer_COPD_up_end,"/")
monocytes_DAVID_COPD_up_end_network <- function_creat_network_form_ClusterProlifer_or_DAVID(monocytes_DAVID_COPD_up_end,", ")
table(monocytes_clusterProlifer_COPD_up_end_network$pathway)

monocytes_clusterProlifer_COPD_up_end_network <- monocytes_clusterProlifer_COPD_up_end_network[monocytes_clusterProlifer_COPD_up_end_network$gene %in% mast_result$monocytes$deg$COPD_up,]










