

### slingshot

library(slingshot)
library(SingleCellExperiment)
library(destiny)

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


ggplot(all_sim,aes(x=DC1,y=DC2))+
  geom_point(aes(col=as.factor(lv2)))+
  geom_line(data=as.data.frame(curve),mapping = aes(x=DC1,y=DC2))+
  theme_classic()+
  facet_wrap(~ group)
}

seurat_sub1 <- subset(seurat_integrated_add_new_metadata, subset = cell.type.js %in% c("alveolar-type-1-cells","alveolar-type-2-cells") )
counts <- seurat_sub1@assays$SCT@counts[seurat_sub1@assays$integrated@var.features,]

s <- apply(counts,1,sd)
counts <- counts[s!=0,]

seurat_slingshot1 <- SingleCellExperiment(assays = List(counts = counts),
                                         colData=seurat_sub1@meta.data)
assays(seurat_slingshot1)$norm <- counts

pca <- prcomp(t(assays(seurat_slingshot1)$norm), scale. = T)
rd1 <- pca$x[,1:2]
reducedDims(seurat_slingshot1) <- SimpleList(pca = rd1)

colData(seurat_slingshot1)$GMM <- seurat_sub1$cell.type.js
cl2 <- kmeans(rd1, centers = 2)$cluster
colData(seurat_slingshot1)$kmeans <- cl2
sim_at2_at1 <- slingshot(seurat_slingshot1, clusterLabels=seurat_sub1$cell.type.js,reducedDim = 'pca',start.clus="alveolar-type-2-cells",
                         end.clus="alveolar-type-1-cells") #clusterLabels = 'GMM',,start.clus="club-cells"
summary(sim_at2_at1$slingPseudotime_1)
plot(reducedDims(sim_at2_at1)$pca, col = as.factor(seurat_sub1$cell.type.js), pch=16, asp = 1,cex=0.5,colour=c("#619CFF","#FF67A4"))
lines(SlingshotDataSet(sim_at2_at1), lwd=1, col='red')

all_sim1 <- cbind(reducedDims(sim_at2_at1)$pca,seurat_sub1@meta.data,pseudotime=sim_at2_at1$slingPseudotime_1)

ggplot(all_sim1,aes(PC1,PC2))+
  geom_point(aes(col=as.factor(cell.type.js)))+
  geom_line(data=as.data.frame(slingCurves(sim_at2_at1)$curve1$s),mapping = aes(x=PC1,y=PC2))+
  scale_colour_manual(values=c("#619CFF","#FF67A4"))+
  theme_classic()


ggplot(all_sim1)+
  geom_density(data=filter(all_sim1,cell.type.js %in% c("alveolar-type-1-cells","alveolar-type-2-cells")),aes(x=pseudotime,col=cell.type.js,linetype=diagnosis.js),adjust=4 )+
  scale_colour_manual(values=c("#619CFF","#FF67A4"))+
  theme_classic()

ggplot(all_sim1)+
  geom_boxplot(aes(x=cell.type.js,y=pseudotime,fill=cell.type.js) )+
  scale_fill_manual(values=c("#619CFF","#FF67A4"))+
  theme_classic()

deg <- read.table("MAST_random_diagnosis.txt",sep='\t',header = T)

at2_deg <- deg[deg$cluster=="alveolar type 2 cells" & deg$fdr<0.05,]$gene

sim_at2_at1_fit <- fitGAM(as.matrix(counts),sds=SlingshotDataSet(sim_at2_at1))
ATres_at <- associationTest(sim_at2_at1_fit)

topgenes <- rownames(ATres_at[order(ATres_at$waldStat,decreasing = T), ])[1:100]
pst.ord <- order(sim_at2_at1$slingPseudotime_1, na.last = NA)
heatdata <- assays(sim_at2_at1)$counts[topgenes, pst.ord]
heatclus <- sim_at2_at1$GMM[pst.ord]

heatmap(log1p(as.matrix(heatdata[intersect(at2_deg,topgenes),])), Colv = NA,ColSideColors = c("#619CFF","#FF67A4")[heatclus],scale="row")

plot(heatdata["KLF4",])
KLF4 <- cbind.data.frame(exp=heatdata["HES1",],index=c(1:8945)) 
ggplot(KLF4,aes(index,exp))+
  geom_smooth()
gene2 <- cbind.data.frame(exp=assays(sim_at2_at1)$counts["KLF4",],pseudotime=sim_at2_at1$slingPseudotime_1) 
ggplot(gene2,aes(pseudotime,exp))+labs(title="KLF4")+
  geom_smooth(method="loess",span = 1)+theme_classic()


seurat_sub <- subset(seurat_integrated_add_new_metadata, subset = cell.type.js %in% c("alveolar-type-1-cells","alveolar-type-2-cells","club-cells","ciliated-cells") )
counts <- seurat_sub@assays$SCT@counts[seurat_sub@assays$integrated@var.features,]
s <- apply(counts,1,sd)
counts <- counts[s!=0,]
seurat_slingshot <- SingleCellExperiment(assays = List(counts = counts),
                                         colData=seurat_sub@meta.data)
assays(seurat_slingshot)$norm <- seurat_sub@assays$integrated@data
assays(seurat_slingshot)$norm <- counts

pca <- prcomp(t(assays(seurat_slingshot)$norm), scale. = T)
rd1 <- pca$x[,1:2]
reducedDims(seurat_slingshot) <- SimpleList(pca = rd1)

colData(seurat_slingshot)$GMM <- seurat_sub$cell.type.js
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(seurat_slingshot)$kmeans <- cl2
sim_at2 <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="alveolar-type-2-cells") #clusterLabels = 'GMM',,start.clus="club-cells"
sim_at1 <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="alveolar-type-1-cells")
sim_cili <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="ciliated-cells")
sim_loess <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',smoother="loess") #clusterLabels = 'GMM',,start.clus="club-cells"
sim_orig <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca')
sim_club <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="club-cells")
sim_club_loess <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="club-cells",smoother="loess")
sim_kmean <- slingshot(seurat_slingshot, clusterLabels="kmeans",reducedDim = 'pca',start.clus="1")
sim_at2_at1 <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="alveolar-type-2-cells",
                         end.clus="alveolar-type-1-cells") #clusterLabels = 'GMM',,start.clus="club-cells"
sim_club_noextend <- slingshot(seurat_slingshot, clusterLabels=seurat_sub$cell.type.js,reducedDim = 'pca',start.clus="club-cells",extend="n")


summary(sim_club$slingPseudotime_3)
col <- as.character(seurat_sub$cell.type.js)
col[col=="alveolar-type-1-cells"]="#619CFF";col[col=="alveolar-type-2-cells"]="#FF67A4";col[col=="ciliated-cells"]="#00BCD8";col[col=="club-cells"]="#C99800"


plot(reducedDims(sim_club)$pca, col = col, pch=16, asp = 1,cex=0.5,xlim=c(-5,60))
lines(SlingshotDataSet(sim_club), lwd=1, col='black',lty="dashed")

all_sim <- cbind(reducedDims(sim)$pca,seurat_sub@meta.data,pseudotime=sim$slingPseudotime_1)
all_sim <- cbind(reducedDims(sim_club)$pca,seurat_sub@meta.data,pseudotime1=sim_club$slingPseudotime_1,
                 pseudotime2=sim_club$slingPseudotime_2,pseudotime3=sim_club$slingPseudotime_3)

curve <- slingCurves(sim)$curve1$s
curve2 <- slingCurves(sim)$curve2$s

ggplot(all_sim,aes(PC1,PC2))+
  geom_point(aes(col=as.factor(cell.type.js)))+
  geom_point(data=as.data.frame(slingCurves(sim_club)$curve1$s),mapping = aes(x=PC1,y=PC2),size=0.5)+
  geom_point(data=as.data.frame(slingCurves(sim_club)$curve2$s),mapping = aes(x=PC1,y=PC2),size=0.5)+
  geom_point(data=as.data.frame(slingCurves(sim_club)$curve3$s),mapping = aes(x=PC1,y=PC2),size=0.5)+
  xlim(-15,60)+
  scale_colour_manual(values=c("#619CFF","#FF67A4","#00BCD8","#C99800"))+
  theme_classic()

age_group <- all_sim$age.js
age_group[age_group%in%c(27,28,35)]="young";age_group[age_group%in%c(50,63,71,73,75)]="old";
all_sim$age_group <- age_group


ggplot(all_sim)+
  geom_density(data=filter(all_sim,cell.type.js %in% c("alveolar-type-1-cells","club-cells")),aes(x=pseudotime2,col=cell.type.js,linetype=as.factor(age_group)),adjust=4 )+
  scale_colour_manual(values=c("#619CFF","#C99800"))+
  theme_classic()

ggplot(all_sim)+
  geom_boxplot(aes(x=as.factor(age.js),y=pseudotime3) )+
  #scale_fill_manual(values=c("#619CFF","#FF67A4"))+
  theme_classic()


ggplot(all_sim,aes(x=PC1,y=PC2))+
  geom_point(aes(col=as.factor(orig.ident)))+
  #geom_line(data=as.data.frame(curve),mapping = aes(x=PC1,y=PC2))+
  #geom_line(data=as.data.frame(curve2),mapping = aes(x=PC1,y=PC2))+
  theme_classic()

ggplot(all_sim,aes(x=cell.type.js,y=pseudotime3))+geom_boxplot(aes(fill=cell.type.js))+
  theme(axis.text.x = element_text(angle = 60,vjust=0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot(all_sim)+
  geom_jitter(data=filter(all_sim,cell.type.js %in% c("alveolar-type-2-cells","club-cells")),aes(x=diagnosis.js,y=pseudotime1,col=cell.type.js) )+
  theme_classic()

ggplot(all_sim)+
  geom_density(data=filter(all_sim,cell.type.js %in% c("alveolar-type-1-cells")),aes(x=pseudotime4,col=diagnosis.js) )+
  theme_classic()

ggplot(all_sim)+
  geom_violin(data=filter(all_sim,cell.type.js %in% c("club-cells")),aes(x=diagnosis.js,y=pseudotime3) )+
  theme_classic()

ggplot(all_sim)+
  geom_violin(data=filter(all_sim,cell.type.js %in% c("alveolar-type-1-cells")),aes(x=diagnosis.js,y=pseudotime2) )+
  theme_classic()

sim_club_fit <- fitGAM(as.matrix(counts),sds=SlingshotDataSet(sim_club))
ATres_club_lineage <- associationTest(sim_club_fit,lineages = T)
topgenes <- rownames(ATres_club_lineage[order(ATres_club_lineage$waldStat,decreasing = T), ])[1:200]
pst.ord <- order(sim_club$slingPseudotime_3, na.last = NA)
heatdata <- assays(sim_club)$counts[c("SOX2","TPPP3"), pst.ord]
heatclus <- sim_club$GMM[pst.ord]

heatmap(log1p(as.matrix(heatdata)), Colv = NA,ColSideColors = brewer.pal(9,"Set1")[heatclus],scale="row")

topgenes_club_at2 <- rownames(ATres_club_lineage[order(ATres_club_lineage$waldStat_1,decreasing = T), ])[1:400]
topgenes_club_at1 <- rownames(ATres_club_lineage[order(ATres_club_lineage$waldStat_2,decreasing = T), ])[1:400]
topgenes_club_cili <- rownames(ATres_club_lineage[order(ATres_club_lineage$waldStat_3,decreasing = T), ])[1:400]

KLF4 <- cbind.data.frame(exp=heatdata["SOX2",],index=c(1:8945)) 
ggplot(KLF4,aes(index,exp))+
  geom_smooth()
gene2 <- cbind.data.frame(exp=assays(sim_club)$counts["CTSS",],pseudotime1=sim_club$slingPseudotime_1,pseudotime2=sim_club$slingPseudotime_2,pseudotime3=sim_club$slingPseudotime_3) 
ggplot(gene2)+labs(title="CTSS")+
  geom_smooth(aes(pseudotime1,exp),method="loess",span = 1,col="#FF67A4")+
  geom_smooth(aes(pseudotime2,exp),method="loess",span = 1,col="#619CFF")+
  geom_smooth(aes(pseudotime3,exp),method="loess",span = 1,col="#00BCD8")+
  theme_classic()


deg_age <- read.table("MAST_random_age.txt",sep='\t',header = T)

club_deg <- deg_age[deg_age$cluster=="club cells" & deg_age$fdr<0.05,]$gene

deg_age[deg_age$cluster=="club cells" & deg_age$gene%in%intersect(topgenes,club_deg),]

deg_age[deg_age$cluster=="club cells" & deg_age$gene%in%intersect(intersect(intersect(topgenes_club_at2,topgenes_club_at1),topgenes_club_cili),club_deg),]


heatmap(log1p(as.matrix(heatdata[intersect(topgenes,club_deg),])), Colv = NA,ColSideColors = brewer.pal(9,"Set1")[heatclus],scale="row")

write(intersect(topgenes,club_deg),"testgene.txt",sep="\n")

plot(heatdata["KLF4",])
KLF4 <- cbind.data.frame(exp=heatdata["HES1",],index=c(1:8945)) 
ggplot(KLF4,aes(index,exp))+
  geom_smooth()
gene2 <- cbind.data.frame(exp=assays(sim_at2_at1)$counts["KLF4",],pseudotime=sim_at2_at1$slingPseudotime_1) 
ggplot(gene2,aes(pseudotime,exp))+labs(title="KLF4")+
  geom_smooth(method="loess",span = 1)+theme_classic()


seurat_sub2 <- subset(seurat_integrated_add_new_metadata, subset = cell.type.js %in% c("alveolar-type-2-cells","club-cells") )
counts <- seurat_sub2@assays$SCT@counts[seurat_sub2@assays$integrated@var.features,]
s <- apply(counts,1,sd)
counts <- counts[s!=0,]
seurat_slingshot2 <- SingleCellExperiment(assays = List(counts = counts),
                                       colData=seurat_sub2@meta.data)
assays(seurat_slingshot2)$norm <- counts
pca <- prcomp(t(assays(seurat_slingshot2)$norm), scale. = T)
rd1 <- pca$x[,1:2]
reducedDims(seurat_slingshot2) <- SimpleList(pca = rd1)
colData(seurat_slingshot2)$GMM <- seurat_sub2$cell.type.js
cl2 <- kmeans(rd1, centers = 2)$cluster
colData(seurat_slingshot2)$kmeans <- cl2
sim_clue_at2 <- slingshot(seurat_slingshot2, clusterLabels=seurat_sub2$cell.type.js,reducedDim = 'pca',start.clus="club-cells",
                          end.clus="alveolar-type-2-cells") #clusterLabels = 'GMM',,start.clus="club-cells"
summary(sim_clue_at2$slingPseudotime_1)
plot(reducedDims(sim_clue_at2)$pca, col = as.factor(seurat_sub2$cell.type.js), pch=16, asp = 1,cex=0.5,colour=c("#619CFF","#C99800"))
lines(SlingshotDataSet(sim_clue_at2), lwd=1, col='red')
all_sim2 <- cbind(reducedDims(sim_clue_at2)$pca,seurat_sub2@meta.data,pseudotime=sim_clue_at2$slingPseudotime_1)
age_group <- all_sim2$age.js
age_group[age_group%in%c(27,28,35)]="young";age_group[age_group%in%c(50,63,71,73,75)]="old";
all_sim2$age_group <- age_group
ggplot(all_sim2)+
  geom_density(aes(x=pseudotime,col=cell.type.js,linetype=as.factor(age_group)),adjust=6 )+
  scale_colour_manual(values=c("#FF67A4","#C99800"))+
  theme_classic()
ggplot(all_sim2)+
  geom_boxplot(data=filter(all_sim2,cell.type.js=="club-cells"),aes(x=age_group,y=pseudotime,fill=age_group),outlier.colour="white" )+
  ylim(60,110)+
  #scale_colour_manual(values=c("#00BCD8","#C99800"))+
  theme_classic()

all_sim2_1 <- filter(all_sim2,cell.type.js=="club-cells")

ggplot(all_sim2_1,aes(x=pseudotime,y=age_group,fill=..x..)) +
  geom_density_ridges_gradient(bandwidth = 4,scale=3) +
  theme_ridges() + 
  scale_fill_viridis_c() #+
  #theme(legend.position = "none")

q1=ggplot(data=filter(all_sim2,cell.type.js=="club-cells"),aes(y=age_group,x=pseudotime,fill=age_group)) +
geom_density_ridges(bandwidth=5,quantile_lines = TRUE, quantiles = 2) +
theme_ridges(grid = FALSE, center_axis_labels = T) +
scale_x_continuous(expand = c(0, 0)) +
scale_y_discrete(expand = c(0, 0)) +
coord_cartesian(clip = "off") +
labs(title = 'club-cells to AT2')


seurat_sub3 <- subset(seurat_integrated_add_new_metadata, subset = cell.type.js %in% c("alveolar-type-1-cells","club-cells") )
counts <- seurat_sub3@assays$SCT@counts[seurat_sub3@assays$integrated@var.features,]
s <- apply(counts,1,sd)
counts <- counts[s!=0,]
seurat_slingshot3 <- SingleCellExperiment(assays = List(counts = counts),
                                          colData=seurat_sub3@meta.data)
assays(seurat_slingshot3)$norm <- counts
pca <- prcomp(t(assays(seurat_slingshot3)$norm), scale. = T)
rd1 <- pca$x[,1:2]
reducedDims(seurat_slingshot3) <- SimpleList(pca = rd1)
colData(seurat_slingshot3)$GMM <- seurat_sub3$cell.type.js
cl2 <- kmeans(rd1, centers = 2)$cluster
colData(seurat_slingshot3)$kmeans <- cl2
sim_clue_at1 <- slingshot(seurat_slingshot3, clusterLabels=seurat_sub3$cell.type.js,reducedDim = 'pca',start.clus="club-cells",
                          end.clus="alveolar-type-1-cells") #clusterLabels = 'GMM',,start.clus="club-cells"
all_sim3 <- cbind(reducedDims(sim_clue_at1)$pca,seurat_sub3@meta.data,pseudotime=sim_clue_at1$slingPseudotime_1)
age_group <- all_sim3$age.js
age_group[age_group%in%c(27,28,35)]="young";age_group[age_group%in%c(50,63,71,73,75)]="old";
all_sim3$age_group <- age_group
ggplot(all_sim3)+
  geom_density(aes(x=pseudotime,col=cell.type.js,linetype=as.factor(age_group)),adjust=6 )+
  scale_colour_manual(values=c("#619CFF","#C99800"))+
  theme_classic()
ggplot(all_sim3)+
  geom_boxplot(data=filter(all_sim3,cell.type.js=="club-cells"),aes(x=age_group,y=pseudotime,fill=age_group),outlier.colour="white" )+
  ylim(35,60)+
  #scale_colour_manual(values=c("#00BCD8","#C99800"))+
  theme_classic()

ggplot(data=filter(all_sim3,cell.type.js=="club-cells"),aes(y=age_group,x=pseudotime,fill=age_group)) +
  geom_density_ridges(bandwidth=3,quantile_lines = TRUE, quantiles = 2) +
  theme_ridges(grid = FALSE, center_axis_labels = T) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = 'club-cells to AT1')

seurat_sub4 <- subset(seurat_integrated_add_new_metadata, subset = cell.type.js %in% c("ciliated-cells","club-cells") )
counts <- seurat_sub4@assays$SCT@counts[seurat_sub4@assays$integrated@var.features,]
s <- apply(counts,1,sd)
counts <- counts[s!=0,]
seurat_slingshot4 <- SingleCellExperiment(assays = List(counts = counts),
                                          colData=seurat_sub4@meta.data)
assays(seurat_slingshot4)$norm <- counts
pca <- prcomp(t(assays(seurat_slingshot4)$norm), scale. = T)
rd1 <- pca$x[,1:2]
reducedDims(seurat_slingshot4) <- SimpleList(pca = rd1)
colData(seurat_slingshot4)$GMM <- seurat_sub4$cell.type.js
cl2 <- kmeans(rd1, centers = 2)$cluster
colData(seurat_slingshot4)$kmeans <- cl2
sim_clue_cili <- slingshot(seurat_slingshot4, clusterLabels=seurat_sub4$cell.type.js,reducedDim = 'pca',start.clus="club-cells",
                          end.clus="ciliated-cells") #clusterLabels = 'GMM',,start.clus="club-cells"
all_sim4 <- cbind(reducedDims(sim_clue_cili)$pca,seurat_sub4@meta.data,pseudotime=sim_clue_cili$slingPseudotime_1)
age_group <- all_sim4$age.js
age_group[age_group%in%c(27,28,35)]="young";age_group[age_group%in%c(50,63,71,73,75)]="old";
all_sim4$age_group <- age_group
ggplot(all_sim4)+
  geom_density(aes(x=pseudotime,col=cell.type.js,linetype=as.factor(age_group)),adjust=6 )+
  scale_colour_manual(values=c("#00BCD8","#C99800"))+
  theme_classic()
ggplot(all_sim4)+
  geom_boxplot(data=filter(all_sim4,cell.type.js=="club-cells"),aes(x=age_group,y=pseudotime,fill=age_group),outlier.colour="white" )+
  ylim(5,25)+
  #scale_colour_manual(values=c("#00BCD8","#C99800"))+
  theme_classic()

toppDEG <- list(at2=topgenes,club=topgenes_club,club_at2=topgenes_club_at2,club_at1=topgenes_club_at1,club_cili=topgenes_club_cili)
DEG_compare=compareCluster(geneCluster = toppDEG, 
                           fun = "enrichGO",
                           keyType = "SYMBOL",
                           OrgDb    = org.Hs.eg.db,
                           ont      = "BP",
                           pAdjustMethod = "BH")
DEG_compare_sim <- simplify(DEG_compare, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(DEG_compare_sim,font.size=8,showCategory=10,title="Top 10 GO terms")




library("scales")
hex_codes1 <- hue_pal()(8)
show_col(hex_codes1)  

seurat_slingshot.label <- list(seurat_integrated$integrated_snn_res.0.1,seurat_integrated$lv3)

function_getLineages <- function(label,start.cluster,reduction,end.clus){
  if(reduction == "umap")
    reducedDim <- reducedDims(seurat_slingshot)$umap
  if(reduction == "pca")
    reducedDim <- reducedDims(seurat_slingshot)$pca
  lineages <- getLineages( reducedDim, clusterLabels = label,start.clus=start.cluster,end.clus=end.clus)
  #plot( reducedDim, col = data_col[label], asp = 1, pch = 16,cex=0.1)
  plot( reducedDim, col = label, asp = 1, pch = 16,cex=0.1)
  lines(lineages, lwd = 1, col = 'black', show.constraints = TRUE)
  lineages
}
function_getCurves <- function(label,lineages,reduction){
  curves <- getCurves(lineages,approx_points = 100,extend = 'n')
  plot(reduction, col = label, asp = 1, pch = 16,cex=0.1)
  lines(curves, lwd = 1, col = 'black')
  curves
}  


par(mfrow=c(1,1))
slingshot_lineages<- lapply(seurat_slingshot.label[1],function_getLineages,start.cluster="3",reduction="pca",end.clus = "0")
slingshot_curves <- function_getCurves(label=seurat_slingshot.label[1],lineages=slingshot_lineages[[1]],reduction=reducedDims(seurat_slingshot)$pca)
plot(reducedDims(seurat_slingshot)$pca, col = seurat_slingshot.label[[1]], asp = 1, pch = 16,cex=0.1)
lines(slingshot_curves, lwd = 1, col = 'black')








