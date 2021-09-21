library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(limma)
library(patchwork)

data<-Read10X(data.dir ="Sample_2812-CL-1_GCACTGAG-TTCACGCA/filtered_feature_bc_matrix")
Control <- CreateSeuratObject(counts = data, project = "Vaccine", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="Control"))
data<-Read10X(data.dir ="Sample_2450-CL-2_AAGATTGG-AAATCCCG/outs/filtered_feature_bc_matrix")
IONP <- CreateSeuratObject(counts = data, project = "Vaccine", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="IONP + PD-1"))
data<-Read10X(data.dir ="Sample_2450-CL-3_TGTAGTCA-TACGATCA/outs/filtered_feature_bc_matrix")
VSM <- CreateSeuratObject(counts = data, project = "Vaccine", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="VSM + PD-1"))
data<-Read10X(data.dir ="Sample_2450-CL-5_GTGGATCA-CAGGGTTG/outs/filtered_feature_bc_matrix")
VSML <- CreateSeuratObject(counts = data, project = "Vaccine", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="VSM&LIGHT + PD-1"))

BC.list<-list(Control=Control,IONP=IONP,VSM=VSM,VSML=VSML)
for (i in 1:length(BC.list)) {
  BC.list[[i]] <- NormalizeData(BC.list[[i]], verbose = FALSE)
  BC.list[[i]] <- FindVariableFeatures(BC.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}
BC.anchors <- FindIntegrationAnchors(object.list = BC.list, dims = 1:15)
BC.integrated <- IntegrateData(anchorset = BC.anchors, dims = 1:15)
DefaultAssay(BC.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
BC.integrated <- ScaleData(BC.integrated, verbose = FALSE)
BC.integrated <- RunPCA(BC.integrated, npcs = 15, verbose = FALSE)

#Cell Cycle Regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
BC.integrated <- CellCycleScoring(BC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
BC.integrated <- ScaleData(BC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(BC.integrated))
##BC.integrated <- RunPCA(BC.integrated, features = c(s.genes, g2m.genes))
BC.integrated <- RunPCA(BC.integrated, npcs = 15, verbose = FALSE)

# t-SNE and Clustering
BC.integrated <- RunUMAP(BC.integrated, reduction = "pca", dims = 1:15)
BC.integrated <- FindNeighbors(BC.integrated, reduction = "pca", dims = 1:15)
BC.integrated <- FindClusters(BC.integrated, resolution = 0.4)

#save data
saveRDS(BC.integrated,"BC.integrated.rds")

#cell type
DefaultAssay(BC.integrated) <- "integrated"
BC.integrated <- RenameIdents(BC.integrated, `0` = "Other T Cells", `1` = "Other T Cells", `2` = "CD4+ T Cells", 
                                `3` = "B Cells", `4` = "CD8+ T Cells", `5` = "B Cells", `6` = "Regulatory T Cells", `7` = "Innate Immune Cells", `8` = "Regulatory T Cells", `9` = "Macrophages", 
                                `10` = "Other T Cells", `11` = "Monocytes", `12` = "Natural Killer Cells", `13` = "Dendritic Cells")
cell_order=c("B Cells","CD4+ T Cells","CD8+ T Cells","Regulatory T Cells","Other T Cells","Macrophages","Dendritic Cells","Natural Killer Cells","Monocytes","Innate Immune Cells")
levels(BC.integrated) <- cell_order
pdf(file = "BC_total.pdf",width = 8,height = 6)
DimPlot(BC.integrated, label = F)
dev.off()

#Graphs
raw.reads=c("Control"=167167,"VSM&LIGHT + PD-1"=31892,"VSM + PD-1"=35358,"IONP + PD-1"=175634,"HER2 + PD-1"=110588)
data=data<-as.matrix(BC.integrated@assays$RNA@data)
raw.read.scale=max(raw.reads)/raw.reads
raw.scale.vector=as.vector(raw.read.scale[as.vector(BC.integrated@meta.data$sample)])
n.data=t((t(data))*raw.scale.vector)
IG=c('Vegfc', 'Ccl19', 'Ccl21a', 'Ccl17', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl12', 'Ccl22', 'Ccl5', 'Ccr7', 'Cxcl10', 'Cxcl13', 'Cxcl12', 'Cxcl9', 'Cxcl16', 'Icam1', 'Icam2', 'Sell', 'Vcam1', 'Ccl8', 'Cxcr4', 'Ifng', 'Tnf', 'Bcl2', 'Tgfb1', 'Tgfb3', 'Lgals1', 'Cd70', 'Cd27', 'Bcl6', 'Cd34', 'Cd40', 'Glycam1', 'MAdcam1', 'Il2', 'Il6', 'Il7', 'Il15', 'Il17', 'Il21', 'Il27', 'Lamp3')
IG=intersect(rownames(data),IG)
for(gene in IG){
  pdf(file = paste0("All_",gene,"_box.pdf"),width = 8,height = 2)
  gene.data<-n.data[gene,]
  meta.data<-BC.integrated@meta.data$sample
  n.plot.data<-data.frame(Expression=as.vector(gene.data),Gene=gene,Sample=meta.data)
  n.plot.data<-n.plot.data[which(n.plot.data$Expression>0),]
  write.csv(n.plot.data,paste0("All_",gene,"_box.csv"),row.names = F)
  print(ggplot(n.plot.data, aes(x=Sample, y=Expression, color=Sample)) +
          geom_boxplot()+ facet_wrap(~ Gene,ncol=1)+ theme(legend.position="none"))
  dev.off()
  
}

sub.data<-subset(BC.integrated,idents="B Cells")
data=data<-as.matrix(sub.data@assays$RNA@data)
raw.read.scale=max(raw.reads)/raw.reads
raw.scale.vector=as.vector(raw.read.scale[as.vector(sub.data@meta.data$sample)])
n.data=t((t(data))*raw.scale.vector)
IG=c('Cxcr4', 'Icam2', 'Sell', 'Vcam1', 'Ccl4', 'Ccr7', 'Cxcl10', 'Ccr5', 'Cxcr3', 'Il2ra', 'Cd38', 'Cd40', 'Ki67', 'Aicda', 'Bcl6', 'Il21', 'Ms4a1', 'Itgal', 'Itga4', 'Cxcr5', 'Il2', 'Il6', 'Il27')
IG=intersect(rownames(data),IG)
for(gene in IG){
  gene.data<-n.data[gene,]
  meta.data<-sub.data@meta.data$sample
  n.plot.data<-data.frame(Expression=as.vector(gene.data),Gene=gene,Sample=meta.data)
  n.plot.data<-n.plot.data[which(n.plot.data$Expression>0),]
  if(nrow(n.plot.data)>0){
  pdf(file = paste0("Bcell_",gene,"_box.pdf"),width = 8,height = 2)
  write.csv(n.plot.data,paste0("Bcell_",gene,"_box.csv"),row.names = F)
  print(ggplot(n.plot.data, aes(x=Sample, y=Expression, color=Sample)) +
          geom_boxplot()+ facet_wrap(~ Gene,ncol=1)+ theme(legend.position="none"))
  dev.off()
  }
}

sub.data<-subset(BC.integrated,seurat_clusters %in% c(0,1,2,4,6,8,10))
data=data<-as.matrix(sub.data@assays$RNA@data)
raw.read.scale=max(raw.reads)/raw.reads
raw.scale.vector=as.vector(raw.read.scale[as.vector(sub.data@meta.data$sample)])
n.data=t((t(data))*raw.scale.vector)
IG=c('Cxcr4', 'Cxr7', 'Ccl4', 'Ccl5', 'Cxcl9', 'Icam2', 'Sell', 'Vcam1', 'Cxcl13', 'Icos', 'Pdcd1', 'Tigit', 'Sgpp2', 'Sh2d1a', 'Fbln7', 'Itgal', 'Itga4', 'Cd40l', 'Cxcr3', 'Fbln7', 'Il21', 'Il2', 'Il6', 'Il27', 'Ccr7', 'Ccr5', 'Gzmb', 'Ifng', 'Tgfb1', 'Tgfb3',"Bcl6","Cd40lg")
IG=intersect(rownames(data),IG)
for(gene in IG){
  gene.data<-n.data[gene,]
  meta.data<-sub.data@meta.data$sample
  n.plot.data<-data.frame(Expression=as.vector(gene.data),Gene=gene,Sample=meta.data)
  n.plot.data<-n.plot.data[which(n.plot.data$Expression>0),]
  if(nrow(n.plot.data)>0){
    pdf(file = paste0("Tcell_",gene,"_box.pdf"),width = 8,height = 2)
    write.csv(n.plot.data,paste0("Tcell_",gene,"_box.csv"),row.names = F)
    print(ggplot(n.plot.data, aes(x=Sample, y=Expression, color=Sample)) +
            geom_boxplot()+ facet_wrap(~ Gene,ncol=1)+ theme(legend.position="none"))
    dev.off()
  }
}

B.gene=c('Ebf1', 'Bank1', 'Mef2c', 'Cr2', 'Fus', 'Pou2f2', 'Dock10', 'Igkc', 'Cd74', 'Mzb1', 'H2−Ab1', 'H2−Eb1', 'Bcl6', 'Fcer1g', 'Ccr2', 'Tyrobp', 'Aif1', 'Lgals3', 'Lgals1', 'Il7r', 'Cd28', 'Cd27', 'Il6ra', 'Spn', 'Cd6')
T4.gene=c('Isg15', 'Ifit1', 'Ifit3', 'Ifit47', 'Ifit206', 'Ifit208', 'Ifi214', 'Oas3', 'Rtp4', 'BST-2', 'Irf7', 'Iigp1', 'Igtp', 'Mndal', 'Stat1', 'Stat2', 'Ly6a', 'Cxcr3', 'Ccl4', 'Pdcd1', 'Ccl5', 'Cxcr6', 'Ifngr1', 'Tnfsf18', 'Itga4')
T8.gene=c('Ifit1', 'Ifit3', 'Igtp', 'Irf7', 'Mndal', 'Bst2', 'Oas3', 'Cxcr3', 'Eomes', 'Bcl2', 'Ccr9', 'S1pr1', 'Klf13', 'Sell', 'Tcf7', 'Ccr7', 'Itga4', 'Lag3', 'Pdcd1', 'Ifng', 'Klrc1', 'Gzmb', 'Gzmk', 'Cxcr6', 'Tnfrsf4', 'Egr3', 'Egr1', 'Cd69', 'Cd82', 'Batf')
B.gene.F=c('Igkc',  'Nme2',  'Actb',  'Ftl1',  'H2−Ab1',  'H2−Eb1',  'Bcl6',  'Fus',  'Pou2f2',  'Stap1',  'Dock10',  'Fcer1g',  'Ccr2',  'Tyrobp',  'Aif1',  'Ifitm3',  'Mif',  'Lgals3',  'Lgals1',  'Il7r',  'Cd28',  'Cd27',  'Il6ra',  'Grap2',  'Spn',  'Cd6',  'Rictor',  'Tcf7',  'Mybl1',  'S1pr2',  'Neil1',  'Hmces',  'Siglecg',  'Cst3',  'Jchain',  'Ebf1',  'Igkc',  'Basp1',  'Cr2',  'Fth1',  'Ighg1',  'Bank1',  'Cd74',  'Eaf2',  'Ciita',  'Cd79a',  'Mef2b',  'Lrrk2',  'Plac8',  'Mzb1',  'Bcl2',  'H2-Ab1',  'Aicda',  'Ikzf3',  'H2-Eb1',  'Cd83',  'Bcl11a',  'Iglc2',  'Cd72',  'Mef2c',  'Iglc1',  'Iglv1',  'Cd37',  'Cd52',  'Xbp1',  'Ms4a4b',  'Fcmr',  'Slpi',  'Lck',  'Iglc3',  'Mcm3',  'Bcl11b',  'H2-DMa',  'Mcm6',  'Malt1',  'Btf3',  'Pdcd1',  'Wdfy4',  'Apoe',  'Mki67',  'Ifng',  'Cd79b',  'Dock8',  'Bcl-2',  'Pfn1',  'Bach2',  'Irf4',  'Ly6d',  'Igll1',  'Irf8',  'Vpreb3',  'Igll5',  'Myc',  'Ctss',  'Ccr6',  'Mcl1',  'Gnai2',  'Tcf3',  'H2-Oa', 'Spib')
T4.gene.F=c('Isg15',  'Ifit47',  'Ifit44',  'Ifit206',  'Ifit203',  'Ifit208',  'Oas3',  'Rtp4',  'BST-2',  'Igtp',  'Stat2',  'Socs3',  'Cxcr3',  'Ccl4',  'Ccl5',  'Tnfsf18',  'Itga4',  'Pim1',  'Igfbp4',  'Ifi27l2a',  'Tnfsf8',  'Pdlim4',  'Ifit1',  'Bcl2a1b',  'Lef1',  'Iigp1',  'Cd82',  'Ccr7',  'Bst2',  'Pdcd1',  'Emp3',  'Ifit3',  'Cxcr6',  'Ikzf2',  'Irf7',  'Cd44',  'Xist',  'Stat1',  'Bcl2a1d',  'Maf',  'Ifi44',  'Ccr2',  'Ahnak',  'Ifi203',  'Cd200',  'Bcl2',  'Ifi47',  'Cd28',  'Tsix',  'Mndal',  'Ctla2a',  'Malt1',  'Ifit3b',  'Cd74',  'Icos',  'Ifi214',  'H2-Aa',  'Wnk1',  'Ifi206',  'H2-Eb1',  'Ctla4',  'Ifi208',  'Igkc',  'Itgae',  'Izumo1r',  'H2-Ab1',  'Tnfrsf4',  'Macf1',  'Bcl11b',  'Il2ra',  'Foxp1',  'Cd37',  'Tnfrsf18',  'Cd55',  'Cxcr5',  'Ly6a',  'Sell',  'Cd40lg',  'Tnfrsf9',  'Cd84',  'Il21',  'Ikzf4',  'Il7r',  'Il4',  'Socs1',  'Il4ra',  'Mki67',  'Socs2',  'Ccr9',  'Tbx21',  'Ifng',  'Ifngr1',  'Havcr2',  'Bcl-2',  'Cxcl10',  'Cd70',  'Cxcl9',  'Cxcl16',  'Cd27',  'Tnfsf14',  'Il15',  'Eomes',  'Lag3',  'Tcf7')
T8.gene.F=c('Ifit1',  'Ifit3',  'Isg15',  'Rtp4',  'Stat1',  'Ly6a',  'Ly6c2',  'Igtp',  'Irf7',  'Mndal',  'Bst2',  'Oas3',  'Iigp1',  'Plac8',  'Irgm1',  'Ifi203',  'Ifi214',  'Ifi208',  'Ifi47',  'Ifi209',  'Il7r',  'Cxcr3',  'Eomes',  'Bcl2',  'Ccr9',  'S1pr1',  'Bcl11b',  'Klf13',  'Satb1',  'Myh9',  'Est1',  'Itgb7',  'Sell',  'Tcf7',  'Ccr7',  'Itga4',  'Il17ra',  'Lag3',  'Pdcd1',  'Ifng',  'Klrc1',  'Fgl2',  'Klrc2',  'Lgals3',  'Ccl4',  'Ccl3',  'Cxcr6',  'Gzmb',  'Klrk1',  'Nkg7',  'Lgals1',  'Ccl5',  'Klrd1',  'Ctla2a',  'Gzmk',  'Tnfrsf18',  'Itgb1',  'Ifngr1',  'Ptma',  'Tnfrsf4',  'Irf4',  'Egr3',  'Tnfrsf9',  'Cd69',  'Cd82',  'Nfkb1',  'Rel',  'Egr1',  'Batf',  'Mif',  'Isg15',  'Ccna2',  'Il4ra',  'Miki67',  'Cd8a',  'Ifit1',  'Cxcr5',  'Dapl1',  'Ifit3',  'Klk8',  'Cd40lg',  'Ms4a4c',  'Cd7',  'Il21',  'Ly6c1',  'Ccr10',  'Il4',  'Maf',  'Ctsw',  'Ctla4',  'Ly6e',  'Tbx21',  'Xcl1',  'Havcr2',  'Mki67',  'Bcl-2',  'Nfkbia',  'Cd3g',  'Cxcl9',  'Nr4a1',  'Ifi44',  'Cxcl10',  'Tagap',  'Gzma',  'Ifi35',  'Cxcl16',  'Nr4a3',  'Il15',  'Nfkbid',  'Il2rb',  'Stat2',  'Cd70',  'Icam1',  'Jaml',  'Birc5',  'Cd27',  'Ccr5',  'Nusap1')  
#Heatmap.genes=read.csv("Heatmap_genes.csv",header = T)
for(celltype in c("B Cells","CD4+ T Cells","CD8+ T Cells")){
  sub.data<-subset(BC.integrated,idents=celltype)
  sample_order=c("Control","IONP + PD-1","VSM + PD-1","VSM&LIGHT + PD-1")
  sub.data@meta.data$sample <- factor(x = sub.data@meta.data$sample, levels = sample_order)
  DefaultAssay(sub.data) <- "integrated"
  sub.data <- ScaleData(sub.data, verbose = FALSE)
  sub.data <- RunPCA(sub.data, npcs = 15, verbose = FALSE)
  sub.data <- RunUMAP(sub.data, reduction = "pca", dims = 1:15)
  sub.data <- FindNeighbors(sub.data, reduction = "pca", dims = 1:15)
  sub.data <- FindClusters(sub.data, resolution = 0.4)
  if(celltype == "B Cells"){
    sub.data <- RenameIdents(sub.data, `0` = "0", `1` = "1", `2` = "0",
                                  `3` = "0", `4` = "2",`5` = "3", `6` = "0")
    temp=as.vector(sub.data@meta.data$seurat_clusters)
    temp.map=data.frame(row.names = 0:6,raw=0:6,new=c(0,1,0,0,2,3,0))
    new.temp=temp.map[temp,]
    sub.data@meta.data$seurat_clusters=new.temp$new
  }else if(celltype == "CD4+ T Cells"){
    sub.data <- RenameIdents(sub.data, `0` = "0", `1` = "0", `2` = "0",
                             `3` = "1", `4` = "0")
    temp=as.vector(sub.data@meta.data$seurat_clusters)
    temp.map=data.frame(row.names = 0:4,raw=0:4,new=c(0,0,0,1,0))
    new.temp=temp.map[temp,]
    sub.data@meta.data$seurat_clusters=new.temp$new    
  }else if(celltype == "CD8+ T Cells"){
    sub.data <- RenameIdents(sub.data, `0` = "0", `1` = "1", `2` = "2",
                             `3` = "3", `4` = "0",`5` = "4", `6` = "0")
    temp=as.vector(sub.data@meta.data$seurat_clusters)
    temp.map=data.frame(row.names = 0:6,raw=0:6,new=c(0,1,2,3,0,4,0))
    new.temp=temp.map[temp,]
    sub.data@meta.data$seurat_clusters=new.temp$new    
  }
  final.table<-data.frame()
  seurat_clusters=unique(sub.data@meta.data$seurat_clusters)
  for(s in unique(sub.data@meta.data$sample)){
    sub.data.temp<-subset(sub.data@meta.data, sample==s)
    celltype.freq<-round(100*table(sub.data.temp$seurat_clusters)/nrow(sub.data.temp),2)
    celltype.freq<-celltype.freq[1:length(seurat_clusters)]
    final.table<-rbind(final.table,celltype.freq)
    
  }
  colnames(final.table)<-paste0("Cluster",0:(length(seurat_clusters)-1))
  row.names(final.table)<-unique(sub.data@meta.data$sample)
  write.csv(final.table,file = paste0(celltype,"_Cluster_proportion.csv"))
  
  pdf(file = paste0(celltype,"_umap.pdf"),width = 8,height = 6)
  print(DimPlot(sub.data, label = F))
  dev.off()
  
  DefaultAssay(sub.data) <- "RNA"
  #markers <- FindAllMarkers(sub.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #top10 <- markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)
  HG=top10$gene
  if(celltype == "B Cells"){HG=intersect(B.gene,row.names(sub.data@assays$RNA@data))}
  if(celltype == "CD4+ T Cells"){HG=intersect(T4.gene,row.names(sub.data@assays$RNA@data))}
  if(celltype == "CD8+ T Cells"){HG=intersect(T8.gene,row.names(sub.data@assays$RNA@data))}
  
  sub.data <- ScaleData(sub.data, verbose = FALSE)
  pdf(file = paste0(celltype,"_heatmap.pdf"),width = 8,height = 0.13*length(HG))
  print(DoHeatmap(sub.data, features = HG) + NoLegend())
  dev.off()
  
  for(gene in HG){
    pdf(file = paste0(celltype,"_",gene,"_umap_all.pdf"),width = 6,height = 6)
    print(FeaturePlot(sub.data, features =gene))
    dev.off()
    # pdf(file = paste0(celltype,"_",gene,"_umap_split.pdf"),width = 30,height = 6)
    # print(FeaturePlot(sub.data, features =gene,split.by = "sample"))
    # dev.off()
  }
  
  DefaultAssay(sub.data) <- "RNA"
  data=data<-as.matrix(sub.data@assays$RNA@data)
  raw.read.scale=max(raw.reads)/raw.reads
  raw.scale.vector=as.vector(raw.read.scale[as.vector(sub.data@meta.data$sample)])
  n.data=t((t(data))*raw.scale.vector)

  if(celltype == "B Cells"){IG=intersect(rownames(data),B.gene.F)}
  if(celltype == "CD4+ T Cells"){IG=intersect(rownames(data),T4.gene.F)}
  if(celltype == "CD8+ T Cells"){IG=intersect(rownames(data),T8.gene.F)}

  for(gene in IG){
    pdf(file = paste0(celltype,"_",gene,"_umap.pdf"),width = 6,height = 6)
    print(FeaturePlot(sub.data, features =gene))
    dev.off()
  }
  

}

BC.integrated<-readRDS("BC.integrated.rds")
sample_order=c("Control","IONP + PD-1","VSM + PD-1","HER2&LIGHT + PD-1","VSM&LIGHT + PD-1")
BC.integrated@meta.data$sample <- factor(x = BC.integrated@meta.data$sample, levels = sample_order)
DefaultAssay(BC.integrated) <- "integrated"
IG=c('Ifng', 'Tnf', 'Bcl-2', 'Cxcl12', 'Tgfb1', 'Tgfb3', 'Lgals1', 'Cxcl9', 'Cxcl10', 'Cxcl16', 'Ccl12', 'Il15', 'Cd70', 'Cd27', 'Ccl2')
data=data<-as.matrix(BC.integrated@assays$RNA@data)
raw.read.scale=max(raw.reads)/raw.reads
raw.scale.vector=as.vector(raw.read.scale[as.vector(BC.integrated@meta.data$sample)])
n.data=t((t(data))*raw.scale.vector)
IG=intersect(rownames(data),IG)
for(gene in IG){
  pdf(file = paste0("All_",gene,"_box.pdf"),width = 8,height = 2)
  gene.data<-n.data[gene,]
  meta.data<-BC.integrated@meta.data$sample
  n.plot.data<-data.frame(Expression=as.vector(gene.data),Gene=gene,Sample=meta.data)
  n.plot.data<-n.plot.data[which(n.plot.data$Expression>0),]
  write.csv(n.plot.data,paste0("All_",gene,"_box.csv"),row.names = F)
  print(ggplot(n.plot.data, aes(x=Sample, y=Expression, color=Sample)) +
          geom_boxplot()+ facet_wrap(~ Gene,ncol=1)+ theme(legend.position="none"))
  dev.off()

}

BC.integrated<-readRDS("BC.integrated.rds")
BC.integrated=subset(BC.integrated,seurat_clusters %in% setdiff(0:17,c(0,1,9,12)))
sample_order=c("Control","IONP + PD-1","VSM + PD-1","HER2&LIGHT + PD-1","VSM&LIGHT + PD-1")
BC.integrated@meta.data$sample <- factor(x = BC.integrated@meta.data$sample, levels = sample_order)
DefaultAssay(BC.integrated) <- "integrated"
IG=c('Ifng', 'Tnf', 'Bcl-2', 'Cxcl12', 'Tgfb1', 'Tgfb3', 'Lgals1', 'Cxcl9', 'Cxcl10', 'Cxcl16', 'Ccl12', 'Il15', 'Cd70', 'Cd27', 'Ccl2')
data=data<-as.matrix(BC.integrated@assays$RNA@data)
raw.read.scale=max(raw.reads)/raw.reads
raw.scale.vector=as.vector(raw.read.scale[as.vector(BC.integrated@meta.data$sample)])
n.data=t((t(data))*raw.scale.vector)
IG=intersect(rownames(data),IG)
for(gene in IG){
  pdf(file = paste0("All_immune_",gene,"_box.pdf"),width = 8,height = 2)
  gene.data<-n.data[gene,]
  meta.data<-BC.integrated@meta.data$sample
  n.plot.data<-data.frame(Expression=as.vector(gene.data),Gene=gene,Sample=meta.data)
  n.plot.data<-n.plot.data[which(n.plot.data$Expression>0),]
  write.csv(n.plot.data,paste0("All_immune_",gene,"_box.csv"),row.names = F)
  print(ggplot(n.plot.data, aes(x=Sample, y=Expression, color=Sample)) +
          geom_boxplot()+ facet_wrap(~ Gene,ncol=1)+ theme(legend.position="none"))
  dev.off()
  
}

